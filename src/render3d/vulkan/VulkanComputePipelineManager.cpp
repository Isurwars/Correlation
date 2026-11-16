/**
 * @file VulkanComputePipelineManager.cpp
 * @brief Implementation of Vulkan compute pipelines for GPU order parameter evaluation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/vulkan/VulkanComputePipelineManager.hpp"
#include "render3d/vulkan/VulkanCommon.hpp"
#include "render3d/vulkan/shaders/HaloApplyCompSpv.hpp"
#include "render3d/vulkan/shaders/HistogramReduceCompSpv.hpp"
#include "render3d/vulkan/shaders/LocalEntropyCompSpv.hpp"
#include "render3d/vulkan/shaders/NeighborGridCompSpv.hpp"
#include "render3d/vulkan/shaders/SteinhardtOrderCompSpv.hpp"
#include "render3d/vulkan/shaders/TenWoldeFrenkelCompSpv.hpp"
#include "render3d/vulkan/shaders/TetrahedralOrderCompSpv.hpp"

#include <vector>

namespace correlation::render3d::vulkan {

#ifdef CORRELATION_USE_VULKAN

void VulkanComputePipelineManager::cleanup() noexcept {
  if (device_ != VK_NULL_HANDLE) {
    for (auto &pipe : pipelines_) {
      if (pipe != VK_NULL_HANDLE) {
        vkDestroyPipeline(device_, pipe, nullptr);
        pipe = VK_NULL_HANDLE;
      }
    }
    for (auto &layout : pipeline_layouts_) {
      if (layout != VK_NULL_HANDLE) {
        vkDestroyPipelineLayout(device_, layout, nullptr);
        layout = VK_NULL_HANDLE;
      }
    }
    for (auto &desc_layout : descriptor_layouts_) {
      if (desc_layout != VK_NULL_HANDLE) {
        vkDestroyDescriptorSetLayout(device_, desc_layout, nullptr);
        desc_layout = VK_NULL_HANDLE;
      }
    }
    device_ = VK_NULL_HANDLE;
  }
  initialized_ = false;
}

bool VulkanComputePipelineManager::createComputePipeline(
    ComputeKernelType type, std::span<const uint32_t> spirv,
    uint32_t push_constant_size, uint32_t num_storage_buffers) {
  const auto idx = static_cast<size_t>(type);

  std::vector<VkDescriptorSetLayoutBinding> bindings(num_storage_buffers);
  for (uint32_t binding_idx = 0; binding_idx < num_storage_buffers; ++binding_idx) {
    bindings.at(binding_idx) = VkDescriptorSetLayoutBinding{
        .binding = binding_idx,
        .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
        .descriptorCount = 1,
        .stageFlags = VK_SHADER_STAGE_COMPUTE_BIT,
        .pImmutableSamplers = nullptr,
    };
  }

  const VkDescriptorSetLayoutCreateInfo layout_info{
      .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .bindingCount = static_cast<uint32_t>(bindings.size()),
      .pBindings = bindings.data(),
  };

  if (vkCreateDescriptorSetLayout(device_, &layout_info, nullptr, &descriptor_layouts_.at(idx)) !=
      VK_SUCCESS) {
    return false;
  }

  const VkPushConstantRange push_range{
      .stageFlags = VK_SHADER_STAGE_COMPUTE_BIT,
      .offset = 0,
      .size = push_constant_size,
  };

  const VkPipelineLayoutCreateInfo pipe_layout_info{
      .sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .setLayoutCount = 1,
      .pSetLayouts = &descriptor_layouts_.at(idx),
      .pushConstantRangeCount = push_constant_size > 0 ? 1U : 0U,
      .pPushConstantRanges = push_constant_size > 0 ? &push_range : nullptr,
  };

  if (vkCreatePipelineLayout(device_, &pipe_layout_info, nullptr, &pipeline_layouts_.at(idx)) !=
      VK_SUCCESS) {
    return false;
  }

  const VkShaderModuleCreateInfo module_info{
      .sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .codeSize = spirv.size() * sizeof(uint32_t),
      .pCode = spirv.data(),
  };

  VkShaderModule shader_module = VK_NULL_HANDLE;
  if (vkCreateShaderModule(device_, &module_info, nullptr, &shader_module) != VK_SUCCESS) {
    return false;
  }

  const VkComputePipelineCreateInfo compute_info{
      .sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO,
      .pNext = nullptr,
      .flags = 0,
      .stage =
          {
              .sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
              .pNext = nullptr,
              .flags = 0,
              .stage = VK_SHADER_STAGE_COMPUTE_BIT,
              .module = shader_module,
              .pName = "main",
              .pSpecializationInfo = nullptr,
          },
      .layout = pipeline_layouts_.at(idx),
      .basePipelineHandle = VK_NULL_HANDLE,
      .basePipelineIndex = -1,
  };

  const VkResult res = vkCreateComputePipelines(device_, VK_NULL_HANDLE, 1, &compute_info, nullptr,
                                                &pipelines_.at(idx));
  vkDestroyShaderModule(device_, shader_module, nullptr);
  return res == VK_SUCCESS;
}

bool VulkanComputePipelineManager::initialize(VkDevice device) {
  cleanup();
  device_ = device;

  // 1. NeighborGrid: 3 buffers, 48-byte push constants
  if (!createComputePipeline(ComputeKernelType::NeighborGrid, shaders::NeighborGridCompSpv, 48, 3)) {
    cleanup();
    return false;
  }

  // 2. SteinhardtOrder: 5 buffers, 48-byte push constants
  if (!createComputePipeline(ComputeKernelType::SteinhardtOrder, shaders::SteinhardtOrderCompSpv, 48, 5)) {
    cleanup();
    return false;
  }

  // 3. TetrahedralOrder: 4 buffers, 40-byte push constants
  if (!createComputePipeline(ComputeKernelType::TetrahedralOrder, shaders::TetrahedralOrderCompSpv, 40, 4)) {
    cleanup();
    return false;
  }

  // 4. LocalEntropy: 4 buffers, 48-byte push constants
  if (!createComputePipeline(ComputeKernelType::LocalEntropy, shaders::LocalEntropyCompSpv, 48, 4)) {
    cleanup();
    return false;
  }

  // 5. TenWoldeFrenkel: 4 buffers, 16-byte push constants
  if (!createComputePipeline(ComputeKernelType::TenWoldeFrenkel, shaders::TenWoldeFrenkelCompSpv, 16, 4)) {
    cleanup();
    return false;
  }

  // 6. HaloApply: 2 buffers, 48-byte push constants
  if (!createComputePipeline(ComputeKernelType::HaloApply, shaders::HaloApplyCompSpv, 48, 2)) {
    cleanup();
    return false;
  }

  // 7. HistogramReduce: 2 buffers, 20-byte push constants
  if (!createComputePipeline(ComputeKernelType::HistogramReduce, shaders::HistogramReduceCompSpv, 20, 2)) {
    cleanup();
    return false;
  }

  initialized_ = true;
  return true;
}

VkPipeline VulkanComputePipelineManager::pipeline(ComputeKernelType type) const noexcept {
  return pipelines_.at(static_cast<size_t>(type));
}

VkPipelineLayout
VulkanComputePipelineManager::pipelineLayout(ComputeKernelType type) const noexcept {
  return pipeline_layouts_.at(static_cast<size_t>(type));
}

VkDescriptorSetLayout
VulkanComputePipelineManager::descriptorSetLayout(ComputeKernelType type) const noexcept {
  return descriptor_layouts_.at(static_cast<size_t>(type));
}

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d::vulkan
