/**
 * @file VulkanRenderEngine.hpp
 * @brief Headless Vulkan offscreen pipeline implementation of IRenderEngine.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "render3d/Camera3D.hpp"
#include "render3d/IRenderEngine.hpp"
#include "render3d/Scene3D.hpp"
#include <memory>
#include <vector>

#ifdef CORRELATION_USE_VULKAN
#include <vulkan/vulkan.h>
#endif

namespace correlation::render3d {

namespace vulkan {
class VulkanDevice;
class VulkanOffscreenTarget;
class VulkanPipelineManager;
class VulkanBufferManager;
} // namespace vulkan

/**
 * @class VulkanRenderEngine
 * @brief Headless Vulkan offscreen rendering pipeline producing RGBA pixel buffers.
 */
class VulkanRenderEngine final : public IRenderEngine {
public:
  VulkanRenderEngine();
  ~VulkanRenderEngine() override;

  VulkanRenderEngine(const VulkanRenderEngine &) = delete;
  VulkanRenderEngine &operator=(const VulkanRenderEngine &) = delete;
  VulkanRenderEngine(VulkanRenderEngine &&) noexcept = default;
  VulkanRenderEngine &operator=(VulkanRenderEngine &&) noexcept = default;

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
  void cleanup() noexcept;
  void recordRenderCommands();
  void bindViewportAndDescriptors();
  void drawWireframe();
  void drawBonds();
  void drawAtoms();

  ViewportDimensions dimensions_{.width = 800, .height = 600};
  std::shared_ptr<const Scene3D> scene_;
  Camera3D camera_;
  bool render_bonds_{true};
  bool render_wireframe_{true};
  bool initialized_{false};
  uint32_t active_staging_index_{0};

#ifdef CORRELATION_USE_VULKAN
  std::unique_ptr<vulkan::VulkanDevice> device_;
  std::unique_ptr<vulkan::VulkanOffscreenTarget> target_;
  std::unique_ptr<vulkan::VulkanPipelineManager> pipeline_manager_;
  std::unique_ptr<vulkan::VulkanBufferManager> buffer_manager_;
  VkCommandBuffer render_command_buffer_{VK_NULL_HANDLE};
  VkFence render_fence_{VK_NULL_HANDLE};
#endif

  std::vector<uint8_t> fallback_pixels_;
};

} // namespace correlation::render3d
