/**
 * @file Camera3D.cpp
 * @brief Implementation of orbit camera for 3D molecular visualization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/Camera3D.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>

namespace correlation::render3d {

namespace {

constexpr float K_DEG_TO_RAD = std::numbers::pi_v<float> / 180.0F;
constexpr float K_MIN_ELEVATION = -89.0F;
constexpr float K_MAX_ELEVATION = 89.0F;
constexpr float K_MIN_DISTANCE = 0.5F;
constexpr float K_MAX_DISTANCE = 5000.0F;

[[nodiscard]] math::Vector3<float> safeNormalize(const math::Vector3<float> &vec) noexcept {
  const float len = math::norm(vec);
  if (len < 1e-6F) {
    return {0.0F, 0.0F, 1.0F};
  }
  return vec / len;
}

} // namespace

Camera3D::Camera3D() noexcept = default;

void Camera3D::orbit(float delta_azimuth_deg, float delta_elevation_deg) noexcept {
  azimuth_deg_ += delta_azimuth_deg;
  while (azimuth_deg_ > 360.0F) {
    azimuth_deg_ -= 360.0F;
  }
  while (azimuth_deg_ < 0.0F) {
    azimuth_deg_ += 360.0F;
  }

  elevation_deg_ =
      std::clamp(elevation_deg_ + delta_elevation_deg, K_MIN_ELEVATION, K_MAX_ELEVATION);
}

void Camera3D::pan(float delta_x, float delta_y) noexcept {
  const math::Vector3<float> eye = eyePosition();
  const math::Vector3<float> forward = safeNormalize(target_ - eye);
  const math::Vector3<float> cam_right = safeNormalize(math::cross(forward, world_up_));
  const math::Vector3<float> cam_up = safeNormalize(math::cross(cam_right, forward));

  const float factor = distance_ * 0.0015F;
  target_ = target_ - cam_right * (delta_x * factor) + cam_up * (delta_y * factor);
}

void Camera3D::zoom(float factor) noexcept {
  distance_ = std::clamp(distance_ * (1.0F - factor), K_MIN_DISTANCE, K_MAX_DISTANCE);
}

void Camera3D::fitToBox(const math::Vector3<float> &min_pt,
                        const math::Vector3<float> &max_pt) noexcept {
  target_ = (min_pt + max_pt) * 0.5F;
  const math::Vector3<float> diag = max_pt - min_pt;
  const float radius = std::max(math::norm(diag) * 0.5F, 1.0F);

  const float fov_rad = fov_y_deg_ * K_DEG_TO_RAD;
  const float half_fov = std::max(fov_rad * 0.5F, 0.01F);
  distance_ = (radius / std::sin(half_fov)) * 1.25F;

  z_near_ = std::max(0.1F, distance_ - radius * 2.5F);
  z_far_ = distance_ + radius * 3.5F;
}

void Camera3D::setViewport(float width, float height) noexcept {
  if (height > 0.0F) {
    aspect_ratio_ = width / height;
  }
}

void Camera3D::setFov(float fov_y_deg) noexcept {
  fov_y_deg_ = std::clamp(fov_y_deg, 10.0F, 120.0F);
}

void Camera3D::setTarget(const math::Vector3<float> &target) noexcept { target_ = target; }

math::Vector3<float> Camera3D::eyePosition() const noexcept {
  const float az_rad = azimuth_deg_ * K_DEG_TO_RAD;
  const float el_rad = elevation_deg_ * K_DEG_TO_RAD;

  const float cos_el = std::cos(el_rad);
  const float sin_el = std::sin(el_rad);
  const float sin_az = std::sin(az_rad);
  const float cos_az = std::cos(az_rad);

  const float offset_x = distance_ * cos_el * sin_az;
  const float offset_y = distance_ * sin_el;
  const float offset_z = distance_ * cos_el * cos_az;

  return target_ + math::Vector3<float>(offset_x, offset_y, offset_z);
}

Matrix4<float> Camera3D::viewMatrix() const {
  return Matrix4<float>::lookAt(eyePosition(), target_, world_up_);
}

Matrix4<float> Camera3D::projectionMatrix() const {
  const float fov_rad = fov_y_deg_ * K_DEG_TO_RAD;
  return Matrix4<float>::perspective(fov_rad, aspect_ratio_, z_near_, z_far_);
}

Matrix4<float> Camera3D::viewProjectionMatrix() const { return projectionMatrix() * viewMatrix(); }

} // namespace correlation::render3d
