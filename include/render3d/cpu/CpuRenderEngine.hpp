/**
 * @file CpuRenderEngine.hpp
 * @brief CPU software rasterization implementation of IRenderEngine.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "render3d/IRenderEngine.hpp"
#include "render3d/cpu/SoftwareFramebuffer.hpp"
#include "render3d/cpu/SoftwareRasterizer.hpp"

namespace correlation::render3d {

/**
 * @class CpuRenderEngine
 * @brief Offscreen CPU software rasterizer implementing IRenderEngine.
 */
class CpuRenderEngine final : public IRenderEngine {
public:
  /**
   * @brief Constructs CPU render engine with optional rasterizer configuration.
   * @param options Rasterizer lighting and feature flags.
   */
  explicit CpuRenderEngine(RasterizerOptions options = RasterizerOptions{}) noexcept;

  [[nodiscard]] bool initialize(const ViewportDimensions &dimensions) override;
  void resize(const ViewportDimensions &dimensions) override;
  void setScene(std::shared_ptr<const Scene3D> scene) override;
  void setCamera(const Camera3D &camera) override;
  void setFeatureVisibility(bool render_bonds, bool render_wireframe) override;
  [[nodiscard]] std::span<const uint8_t> renderFrame() override;
  [[nodiscard]] ViewportDimensions dimensions() const noexcept override;
  [[nodiscard]] std::string_view backendName() const noexcept override;
  [[nodiscard]] bool isHardwareAccelerated() const noexcept override;

private:
  SoftwareRasterizer rasterizer_;
  SoftwareFramebuffer framebuffer_;
  std::shared_ptr<const Scene3D> scene_;
  Camera3D camera_;
  ViewportDimensions dimensions_{800, 600};
};

} // namespace correlation::render3d
