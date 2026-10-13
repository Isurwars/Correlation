/**
 * @file SoftwareRasterizer.cpp
 * @brief Implementation of CPU software rasterizer for 3D molecular visualization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/cpu/SoftwareRasterizer.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>

namespace correlation::render3d {

namespace {

constexpr float K_DEG_TO_RAD = std::numbers::pi_v<float> / 180.0F;

[[nodiscard]] math::Vector3<float> safeNormalize(const math::Vector3<float> &vec) noexcept {
  const float len = math::norm(vec);
  if (len < 1e-6F) {
    return {0.0F, 0.0F, 1.0F};
  }
  return vec / len;
}

[[nodiscard]] RgbaColor shadeColor(RgbaColor base_color, float normal_x, float normal_y,
                                   float normal_z, const RasterizerOptions &options) noexcept {
  const math::Vector3<float> norm_vec(normal_x, normal_y, normal_z);
  const math::Vector3<float> light_dir = safeNormalize(options.light_dir);
  const math::Vector3<float> view_dir(0.0F, 0.0F, 1.0F);
  const math::Vector3<float> halfway = safeNormalize(light_dir + view_dir);

  const float n_dot_l = std::max(0.0F, math::dot(norm_vec, light_dir));
  const float n_dot_h = std::max(0.0F, math::dot(norm_vec, halfway));

  const float diffuse_term = options.diffuse * n_dot_l;
  const float specular_term = options.specular * std::pow(n_dot_h, options.shininess);
  const float intensity = options.ambient + diffuse_term;

  const auto red =
      static_cast<int>(static_cast<float>(base_color.r) * intensity + 255.0F * specular_term);
  const auto green =
      static_cast<int>(static_cast<float>(base_color.g) * intensity + 255.0F * specular_term);
  const auto blue =
      static_cast<int>(static_cast<float>(base_color.b) * intensity + 255.0F * specular_term);

  return RgbaColor{
      .r = static_cast<uint8_t>(std::clamp(red, 0, 255)),
      .g = static_cast<uint8_t>(std::clamp(green, 0, 255)),
      .b = static_cast<uint8_t>(std::clamp(blue, 0, 255)),
      .a = base_color.a,
  };
}

[[nodiscard]] bool projectPoint(const math::Vector3<float> &world_pt,
                                const Matrix4<float> &view_proj, uint32_t width, uint32_t height,
                                math::Vector3<float> &out_screen) noexcept {
  const Vector4<float> clip = view_proj * Vector4<float>(world_pt, 1.0F);
  if (clip.w() <= 1e-4F) {
    return false;
  }

  const float inv_w = 1.0F / clip.w();
  const float ndc_x = clip.x() * inv_w;
  const float ndc_y = clip.y() * inv_w;
  const float ndc_z = clip.z() * inv_w;

  const float screen_x = (ndc_x * 0.5F + 0.5F) * static_cast<float>(width);
  const float screen_y = (1.0F - (ndc_y * 0.5F + 0.5F)) * static_cast<float>(height);

  out_screen = math::Vector3<float>(screen_x, screen_y, ndc_z);
  return true;
}

void drawLineSegment(math::Vector3<float> start_pt, math::Vector3<float> end_pt, RgbaColor color,
                     SoftwareFramebuffer &framebuffer) noexcept {
  const auto x0_coord = static_cast<int>(std::round(start_pt.x()));
  const auto y0_coord = static_cast<int>(std::round(start_pt.y()));
  const auto x1_coord = static_cast<int>(std::round(end_pt.x()));
  const auto y1_coord = static_cast<int>(std::round(end_pt.y()));

  const int delta_x = std::abs(x1_coord - x0_coord);
  const int delta_y = std::abs(y1_coord - y0_coord);
  const int steps = std::max(delta_x, delta_y);

  if (steps == 0) {
    if (x0_coord >= 0 && y0_coord >= 0) {
      framebuffer.testAndSetPixel(static_cast<uint32_t>(x0_coord), static_cast<uint32_t>(y0_coord),
                                  start_pt.z(), color);
    }
    return;
  }

  const float x_inc = static_cast<float>(x1_coord - x0_coord) / static_cast<float>(steps);
  const float y_inc = static_cast<float>(y1_coord - y0_coord) / static_cast<float>(steps);
  const float z_inc = (end_pt.z() - start_pt.z()) / static_cast<float>(steps);

  auto current_x = static_cast<float>(x0_coord);
  auto current_y = static_cast<float>(y0_coord);
  float current_z = start_pt.z();

  for (int step = 0; step <= steps; ++step) {
    const auto px_coord = static_cast<int>(std::round(current_x));
    const auto py_coord = static_cast<int>(std::round(current_y));
    if (px_coord >= 0 && py_coord >= 0) {
      framebuffer.testAndSetPixel(static_cast<uint32_t>(px_coord), static_cast<uint32_t>(py_coord),
                                  current_z, color);
    }
    current_x += x_inc;
    current_y += y_inc;
    current_z += z_inc;
  }
}

void rasterizeSphere(const AtomInstance &atom, const Matrix4<float> &view,
                     const Matrix4<float> &view_proj, float fov_rad,
                     const RasterizerOptions &options, SoftwareFramebuffer &framebuffer) noexcept {
  const math::Vector3<float> p_view = view.transformPoint(atom.position);
  const float view_dist = -p_view.z();
  if (view_dist <= 0.1F) {
    return;
  }

  math::Vector3<float> screen_center;
  if (!projectPoint(atom.position, view_proj, framebuffer.width(), framebuffer.height(),
                    screen_center)) {
    return;
  }

  const float tan_half_fov = std::tan(fov_rad * 0.5F);
  const float screen_r =
      (atom.radius * static_cast<float>(framebuffer.height())) / (2.0F * tan_half_fov * view_dist);
  if (screen_r < 0.5F) {
    return;
  }

  const auto width_bound = static_cast<int>(framebuffer.width());
  const auto height_bound = static_cast<int>(framebuffer.height());

  const auto min_x =
      std::clamp(static_cast<int>(std::floor(screen_center.x() - screen_r)), 0, width_bound - 1);
  const auto max_x =
      std::clamp(static_cast<int>(std::ceil(screen_center.x() + screen_r)), 0, width_bound - 1);
  const auto min_y =
      std::clamp(static_cast<int>(std::floor(screen_center.y() - screen_r)), 0, height_bound - 1);
  const auto max_y =
      std::clamp(static_cast<int>(std::ceil(screen_center.y() + screen_r)), 0, height_bound - 1);

  const float r_squared = screen_r * screen_r;
  const float inv_r = 1.0F / screen_r;

  for (int pixel_y = min_y; pixel_y <= max_y; ++pixel_y) {
    const float diff_y = static_cast<float>(pixel_y) - screen_center.y();
    const float diff_y_sq = diff_y * diff_y;

    for (int pixel_x = min_x; pixel_x <= max_x; ++pixel_x) {
      const float diff_x = static_cast<float>(pixel_x) - screen_center.x();
      const float dist_sq = diff_x * diff_x + diff_y_sq;

      if (dist_sq <= r_squared) {
        const float norm_x = diff_x * inv_r;
        const float norm_y = -diff_y * inv_r;
        const float norm_z = std::sqrt(std::max(0.0F, 1.0F - norm_x * norm_x - norm_y * norm_y));

        const float z_offset = norm_z * (atom.radius / view_dist) * 0.5F;
        const float depth = screen_center.z() - z_offset;

        const RgbaColor lit_color = shadeColor(atom.color, norm_x, norm_y, norm_z, options);
        framebuffer.testAndSetPixel(static_cast<uint32_t>(pixel_x), static_cast<uint32_t>(pixel_y),
                                    depth, lit_color);
      }
    }
  }
}

} // namespace

SoftwareRasterizer::SoftwareRasterizer(RasterizerOptions options) noexcept : options_(options) {}

void SoftwareRasterizer::render(const Scene3D &scene, const Camera3D &camera,
                                SoftwareFramebuffer &framebuffer) const {
  if (framebuffer.width() == 0 || framebuffer.height() == 0) {
    return;
  }

  framebuffer.clear(options_.background_color);

  const Matrix4<float> view = camera.viewMatrix();
  const Matrix4<float> view_proj = camera.viewProjectionMatrix();
  const float fov_rad = camera.fov() * K_DEG_TO_RAD;

  // 1. Render unit cell wireframe lines
  if (options_.render_wireframe) {
    for (const auto &line : scene.unitCellLines()) {
      math::Vector3<float> p0_screen;
      math::Vector3<float> p1_screen;
      if (projectPoint(line.start, view_proj, framebuffer.width(), framebuffer.height(),
                       p0_screen) &&
          projectPoint(line.end, view_proj, framebuffer.width(), framebuffer.height(), p1_screen)) {
        drawLineSegment(p0_screen, p1_screen, line.color, framebuffer);
      }
    }
  }

  // 2. Render covalent bonds
  if (options_.render_bonds) {
    for (const auto &bond : scene.bonds()) {
      math::Vector3<float> p0_screen;
      math::Vector3<float> p1_screen;
      if (projectPoint(bond.start, view_proj, framebuffer.width(), framebuffer.height(),
                       p0_screen) &&
          projectPoint(bond.end, view_proj, framebuffer.width(), framebuffer.height(), p1_screen)) {
        drawLineSegment(p0_screen, p1_screen, bond.color, framebuffer);
      }
    }
  }

  // 3. Render spherical atoms
  if (options_.render_atoms) {
    for (const auto &atom : scene.atoms()) {
      rasterizeSphere(atom, view, view_proj, fov_rad, options_, framebuffer);
    }
  }
}

} // namespace correlation::render3d
