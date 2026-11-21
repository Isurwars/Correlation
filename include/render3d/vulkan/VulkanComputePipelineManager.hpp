/**
 * @file VulkanComputePipelineManager.hpp
 * @brief Manages Vulkan compute pipelines for GPU order parameter and halo evaluation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#ifdef CORRELATION_USE_VULKAN
#include <vulkan/vulkan.h>
#endif

#include <array>
#include <cstdint>
#include <span>

namespace correlation::render3d::vulkan {

#ifdef CORRELATION_USE_VULKAN

enum class ComputeKernelType : uint8_t {
  NeighborGrid = 0,
  SteinhardtOrder = 1,
  TetrahedralOrder = 2,
  LocalEntropy = 3,
  TenWoldeFrenkel = 4,
  HaloApply = 5,
  HistogramReduce = 6,
  CellListBuild = 7
};

class VulkanComputePipelineManager {
public:
  static constexpr size_t KERNEL_COUNT = 8;

  VulkanComputePipelineManager() = default;
  ~VulkanComputePipelineManager() noexcept { cleanup(); }

  VulkanComputePipelineManager(const VulkanComputePipelineManager &) = delete;
  VulkanComputePipelineManager &operator=(const VulkanComputePipelineManager &) = delete;
  VulkanComputePipelineManager(VulkanComputePipelineManager &&) noexcept = default;
  VulkanComputePipelineManager &operator=(VulkanComputePipelineManager &&) noexcept = default;

  [[nodiscard]] bool initialize(VkDevice device);
  void cleanup() noexcept;

  [[nodiscard]] VkPipeline pipeline(ComputeKernelType type) const noexcept;
  [[nodiscard]] VkPipelineLayout pipelineLayout(ComputeKernelType type) const noexcept;
  [[nodiscard]] VkDescriptorSetLayout descriptorSetLayout(ComputeKernelType type) const noexcept;

private:
  [[nodiscard]] bool createComputePipeline(ComputeKernelType type, std::span<const uint32_t> spirv,
                                           uint32_t push_constant_size,
                                           uint32_t num_storage_buffers);

  VkDevice device_{VK_NULL_HANDLE};
  std::array<VkDescriptorSetLayout, KERNEL_COUNT> descriptor_layouts_{};
  std::array<VkPipelineLayout, KERNEL_COUNT> pipeline_layouts_{};
  std::array<VkPipeline, KERNEL_COUNT> pipelines_{};
  bool initialized_{false};
};

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d::vulkan
