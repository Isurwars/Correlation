/**
 * @file VulkanPipelineManager.hpp
 * @brief Graphics pipeline and descriptor management for Vulkan 3D molecular impostors.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#ifdef CORRELATION_USE_VULKAN
#include <vulkan/vulkan.h>
#endif

#include <span>

namespace correlation::render3d::vulkan {

#ifdef CORRELATION_USE_VULKAN

/**
 * @class VulkanPipelineManager
 * @brief Manages graphics pipelines, layouts, and shader modules for sphere/cylinder impostors and
 * wireframe lines.
 */
class VulkanPipelineManager {
public:
  VulkanPipelineManager() noexcept = default;
  ~VulkanPipelineManager() noexcept { destroy(); }

  VulkanPipelineManager(const VulkanPipelineManager &) = delete;
  VulkanPipelineManager &operator=(const VulkanPipelineManager &) = delete;

  VulkanPipelineManager(VulkanPipelineManager &&other) noexcept { moveFrom(other); }
  VulkanPipelineManager &operator=(VulkanPipelineManager &&other) noexcept {
    if (this != &other) {
      destroy();
      moveFrom(other);
    }
    return *this;
  }

  /**
   * @brief Initializes descriptor set layouts, pipeline layouts, and graphics pipelines.
   * @param device Vulkan logical device handle.
   * @param render_pass Compatible render pass handle.
   * @return True if all pipelines compiled successfully.
   */
  [[nodiscard]] bool initialize(VkDevice device, VkRenderPass render_pass);

  /**
   * @brief Destroys all pipelines, pipeline layouts, and descriptor set layouts.
   */
  void destroy() noexcept;

  [[nodiscard]] VkPipeline spherePipeline() const noexcept { return sphere_pipeline_; }
  [[nodiscard]] VkPipeline cylinderPipeline() const noexcept { return cylinder_pipeline_; }
  [[nodiscard]] VkPipeline wireframePipeline() const noexcept { return wireframe_pipeline_; }
  [[nodiscard]] VkPipelineLayout pipelineLayout() const noexcept { return pipeline_layout_; }
  [[nodiscard]] VkDescriptorSetLayout descriptorSetLayout() const noexcept {
    return descriptor_set_layout_;
  }
  [[nodiscard]] bool isValid() const noexcept { return sphere_pipeline_ != VK_NULL_HANDLE; }

private:
  void moveFrom(VulkanPipelineManager &other) noexcept;
  [[nodiscard]] VkShaderModule createShaderModule(std::span<const uint32_t> spv_code) const;
  [[nodiscard]] bool createDescriptorSetLayout();
  [[nodiscard]] bool createPipelineLayout();
  [[nodiscard]] bool createSpherePipeline(VkRenderPass render_pass);
  [[nodiscard]] bool createCylinderPipeline(VkRenderPass render_pass);
  [[nodiscard]] bool createWireframePipeline(VkRenderPass render_pass);

  VkDevice device_{VK_NULL_HANDLE};
  VkDescriptorSetLayout descriptor_set_layout_{VK_NULL_HANDLE};
  VkPipelineLayout pipeline_layout_{VK_NULL_HANDLE};
  VkPipeline sphere_pipeline_{VK_NULL_HANDLE};
  VkPipeline cylinder_pipeline_{VK_NULL_HANDLE};
  VkPipeline wireframe_pipeline_{VK_NULL_HANDLE};
};

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d::vulkan
