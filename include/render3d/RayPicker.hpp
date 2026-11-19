/**
 * @file RayPicker.hpp
 * @brief Ray-sphere picking and camera HUD projection utilities for 3D molecular viewports.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "render3d/Camera3D.hpp"
#include "render3d/Scene3D.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <string>

namespace correlation::render3d {

/**
 * @struct CameraHudMetrics
 * @brief Projected 2D axis endpoints and calibrated scale bar metrics derived from Camera3D.
 */
struct CameraHudMetrics {
  std::array<float, 2> axis_x{1.0F, 0.0F};  ///< Screen-space X axis direction in [-1, 1].
  std::array<float, 2> axis_y{0.0F, -1.0F}; ///< Screen-space Y axis direction in [-1, 1].
  std::array<float, 2> axis_z{0.0F, 0.0F};  ///< Screen-space Z axis direction in [-1, 1].
  float scale_bar_px{64.0F};                ///< Pixel length of the calibrated scale bar.
  std::string scale_bar_label{"5 Å"};       ///< Human-readable scale label.
};

/**
 * @brief Picks the closest visible atom sphere intersected by a viewport click ray.
 * @param[in] camera Active 3D orbit camera.
 * @param[in] scene Active 3D molecular scene.
 * @param[in] norm_x Normalized horizontal viewport coordinate in [0, 1] (left to right).
 * @param[in] norm_y Normalized vertical viewport coordinate in [0, 1] (top to bottom).
 * @return Zero-based index into scene.atoms() of the closest hit atom, or std::nullopt.
 */
[[nodiscard]] std::optional<std::size_t> pickAtom(const Camera3D &camera, const Scene3D &scene,
                                                  float norm_x, float norm_y) noexcept;

/**
 * @brief Computes projected axis triad vectors and a nice-number Angstrom scale bar.
 * @param[in] camera Active 3D orbit camera.
 * @param[in] viewport_height_px Height of the 3D viewport in pixels.
 * @return Populated CameraHudMetrics struct.
 */
[[nodiscard]] CameraHudMetrics computeCameraHudMetrics(const Camera3D &camera,
                                                       float viewport_height_px);

} // namespace correlation::render3d
