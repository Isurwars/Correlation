/**
 * @file RayPicker.cpp
 * @brief Implementation of ray-sphere atom picking and camera HUD projection utilities.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/RayPicker.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <limits>
#include <numbers>

namespace correlation::render3d {

namespace {

constexpr float DEG_TO_RAD = std::numbers::pi_v<float> / 180.0F;
constexpr float TARGET_BAR_PX = 84.0F;
constexpr std::array<float, 12> NICE_ANGSTROM_STEPS{
    0.5F, 1.0F, 2.0F, 5.0F, 10.0F, 20.0F, 50.0F, 100.0F, 200.0F, 500.0F, 1000.0F, 2000.0F};

struct CameraBasis {
  math::Vector3<float> eye;
  math::Vector3<float> forward;
  math::Vector3<float> right;
  math::Vector3<float> up;
  float half_height{1.0F};
  float half_width{1.0F};
};

[[nodiscard]] math::Vector3<float> safeNormalize(const math::Vector3<float> &vec) noexcept {
  const float length = math::norm(vec);
  if (length < 1e-6F) {
    return {0.0F, 0.0F, 1.0F};
  }
  return vec / length;
}

[[nodiscard]] CameraBasis buildCameraBasis(const Camera3D &camera) noexcept {
  const math::Vector3<float> eye = camera.eyePosition();
  const math::Vector3<float> forward = safeNormalize(camera.target() - eye);
  const math::Vector3<float> world_up{0.0F, 1.0F, 0.0F};
  const math::Vector3<float> right = safeNormalize(math::cross(forward, world_up));
  const math::Vector3<float> cam_up = safeNormalize(math::cross(right, forward));

  const float fov_rad = camera.fov() * DEG_TO_RAD;
  const float half_height = std::max(camera.distance() * std::tan(fov_rad * 0.5F), 1e-4F);
  const float half_width = half_height * std::max(camera.aspectRatio(), 1e-4F);

  return CameraBasis{
      .eye = eye,
      .forward = forward,
      .right = right,
      .up = cam_up,
      .half_height = half_height,
      .half_width = half_width,
  };
}

[[nodiscard]] std::optional<float> intersectRaySphere(const math::Vector3<float> &origin,
                                                      const math::Vector3<float> &direction,
                                                      const math::Vector3<float> &center,
                                                      float radius) noexcept {
  const math::Vector3<float> offset = origin - center;
  const float half_b = math::dot(offset, direction);
  const float c_term = math::dot(offset, offset) - (radius * radius);
  const float discriminant = (half_b * half_b) - c_term;
  if (discriminant < 0.0F) {
    return std::nullopt;
  }
  const float sqrt_disc = std::sqrt(discriminant);
  const float near_t = -half_b - sqrt_disc;
  if (near_t > 0.0F) {
    return near_t;
  }
  const float far_t = -half_b + sqrt_disc;
  if (far_t > 0.0F) {
    return far_t;
  }
  return std::nullopt;
}

} // namespace

std::optional<std::size_t> pickAtom(const Camera3D &camera, const Scene3D &scene, float norm_x,
                                    float norm_y) noexcept {
  const auto &atoms = scene.atoms();
  if (atoms.empty() || norm_x < 0.0F || norm_x > 1.0F || norm_y < 0.0F || norm_y > 1.0F) {
    return std::nullopt;
  }

  const CameraBasis basis = buildCameraBasis(camera);
  const float ndc_x = (2.0F * norm_x) - 1.0F;
  const float ndc_y = 1.0F - (2.0F * norm_y);

  math::Vector3<float> ray_origin = basis.eye;
  math::Vector3<float> ray_dir = basis.forward;

  if (camera.isOrthographic()) {
    ray_origin = basis.eye + (basis.right * (ndc_x * basis.half_width)) +
                 (basis.up * (ndc_y * basis.half_height));
  } else {
    const math::Vector3<float> target_pt = (basis.forward * camera.distance()) +
                                           (basis.right * (ndc_x * basis.half_width)) +
                                           (basis.up * (ndc_y * basis.half_height));
    ray_dir = safeNormalize(target_pt);
  }

  std::optional<std::size_t> best_index;
  float best_dist = std::numeric_limits<float>::max();

  for (std::size_t idx = 0; idx < atoms.size(); ++idx) {
    const auto &atom = atoms[idx];
    if (atom.dim_factor <= 0.0F && atom.halo_color.a == 0) {
      continue;
    }
    const float pick_radius = std::max(atom.radius * 1.15F, 0.25F);
    const auto hit_t = intersectRaySphere(ray_origin, ray_dir, atom.position, pick_radius);
    if (hit_t.has_value() && *hit_t < best_dist) {
      best_dist = *hit_t;
      best_index = idx;
    }
  }

  return best_index;
}

CameraHudMetrics computeCameraHudMetrics(const Camera3D &camera, float viewport_height_px) {
  const CameraBasis basis = buildCameraBasis(camera);
  CameraHudMetrics metrics{};

  metrics.axis_x = {basis.right.x(), -basis.up.x()};
  metrics.axis_y = {basis.right.y(), -basis.up.y()};
  metrics.axis_z = {basis.right.z(), -basis.up.z()};

  const float safe_height = std::max(viewport_height_px, 100.0F);
  const float angstroms_per_px = (2.0F * basis.half_height) / safe_height;
  const float target_angstroms = angstroms_per_px * TARGET_BAR_PX;

  float chosen_step = NICE_ANGSTROM_STEPS.back();
  for (const float step : NICE_ANGSTROM_STEPS) {
    if (step >= target_angstroms * 0.7F) {
      chosen_step = step;
      break;
    }
  }

  metrics.scale_bar_px = std::clamp(chosen_step / std::max(angstroms_per_px, 1e-5F), 32.0F, 140.0F);
  if (chosen_step < 1.0F) {
    metrics.scale_bar_label = std::format("{:.1f} Å", chosen_step);
  } else {
    metrics.scale_bar_label = std::format("{:.0f} Å", chosen_step);
  }

  return metrics;
}

} // namespace correlation::render3d
