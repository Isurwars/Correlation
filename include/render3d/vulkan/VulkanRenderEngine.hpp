/**
 * @file VulkanRenderEngine.hpp
 * @brief Headless Vulkan offscreen pipeline implementation of IRenderEngine.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "render3d/IRenderEngine.hpp"
#include <vector>

namespace correlation::render3d {

/**
 * @class VulkanRenderEngine
 * @brief Headless Vulkan offscreen rendering pipeline producing RGBA pixel buffers.
 */
class VulkanRenderEngine final : public IRenderEngine {
public:
  VulkanRenderEngine() = default;
  ~VulkanRenderEngine() override;

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
  ViewportDimensions dimensions_{.width = 800, .height = 600};
  std::shared_ptr<const Scene3D> scene_;
  Camera3D camera_;
  std::vector<uint8_t> frame_pixels_;
  bool initialized_{false};
};

} // namespace correlation::render3d
