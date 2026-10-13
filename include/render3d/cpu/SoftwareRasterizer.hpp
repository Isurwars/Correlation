/**
 * @file SoftwareRasterizer.hpp
 * @brief High-performance CPU software rasterizer with analytical ray-cast spheres and depth
 * testing.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "render3d/Camera3D.hpp"
#include "render3d/Scene3D.hpp"
#include "render3d/cpu/SoftwareFramebuffer.hpp"

namespace correlation::render3d {

/**
 * @struct RasterizerOptions
 * @brief Configuration parameters for software rasterizer lighting and shading.
 */
struct RasterizerOptions {
  RgbaColor background_color{.r = 24, .g = 24, .b = 29, .a = 255}; ///< Dark slate default.
  math::Vector3<float> light_dir{0.35F, 0.55F, 0.75F};             ///< Directional light vector.
  float ambient{0.30F};                                            ///< Ambient lighting factor.
  float diffuse{0.65F};                                            ///< Diffuse lighting factor.
  float specular{0.35F};                                           ///< Specular highlight factor.
  float shininess{24.0F};                                          ///< Specular exponent.
  bool render_atoms{true};                                         ///< Flag to draw atom spheres.
  bool render_bonds{true};                                         ///< Flag to draw bonds.
  bool render_wireframe{true};                                     ///< Flag to draw simulation box.
};

/**
 * @class SoftwareRasterizer
 * @brief Executes software-rendered 3D visualization using depth-tested billboards and lines.
 */
class SoftwareRasterizer {
public:
  /**
   * @brief Constructs software rasterizer with default lighting options.
   * @param options Lighting and rendering options.
   */
  explicit SoftwareRasterizer(RasterizerOptions options = RasterizerOptions{}) noexcept;

  /**
   * @brief Renders the 3D scene from the camera viewpoint into the framebuffer.
   * @param scene Molecular geometry scene graph.
   * @param camera View and projection camera.
   * @param framebuffer Output color and depth buffer.
   */
  void render(const Scene3D &scene, const Camera3D &camera, SoftwareFramebuffer &framebuffer) const;

  [[nodiscard]] const RasterizerOptions &options() const noexcept { return options_; }
  void setOptions(const RasterizerOptions &options) noexcept { options_ = options; }

private:
  RasterizerOptions options_;
};

} // namespace correlation::render3d
