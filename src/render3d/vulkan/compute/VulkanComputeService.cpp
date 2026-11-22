/**
 * @file VulkanComputeService.cpp
 * @brief Implementation of VulkanComputeService GPU-only order parameter evaluation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/vulkan/compute/VulkanComputeService.hpp"
#include "render3d/vulkan/VulkanCommon.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace correlation::render3d::vulkan {

#ifdef CORRELATION_USE_VULKAN

constexpr uint32_t MAX_NEIGHBORS = 64;
constexpr int32_t MAX_GRID_DIM = 64;
constexpr uint32_t EMPTY_CELL_SENTINEL = 0xFFFFFFFFU;

namespace {

struct PushConstantsCellListBuild {
  std::array<float, 4> box_min;
  std::array<float, 4> box_max;
  std::array<int32_t, 4> grid_dims;
  uint32_t num_atoms;
  uint32_t pad0;
  uint32_t pad1;
  uint32_t pad2;
};
static_assert(sizeof(PushConstantsCellListBuild) == 64);

struct PushConstantsNeighborGrid {
  std::array<float, 4> box_min;
  std::array<float, 4> box_max;
  std::array<int32_t, 4> grid_dims;
  float cutoff;
  uint32_t max_neighbors;
  uint32_t num_atoms;
  uint32_t pad;
};
static_assert(sizeof(PushConstantsNeighborGrid) == 64);

[[nodiscard]] std::array<int32_t, 4> computeGridDims(const std::array<float, 4> &box_min,
                                                     const std::array<float, 4> &box_max,
                                                     float cutoff) noexcept {
  const float len_x = box_max[0] - box_min[0];
  const float len_y = box_max[1] - box_min[1];
  const float len_z = box_max[2] - box_min[2];
  if (len_x <= 0.0F || len_y <= 0.0F || len_z <= 0.0F || cutoff <= 0.0F) {
    return {1, 1, 1, 0};
  }
  const auto dim = [cutoff](float length) noexcept -> int32_t {
    const auto count = static_cast<int32_t>(std::floor(length / cutoff));
    return std::clamp(count, 1, MAX_GRID_DIM);
  };
  return {dim(len_x), dim(len_y), dim(len_z), 1};
}

void emitComputeBarrier(VkCommandBuffer cmd_buffer) noexcept {
  const VkMemoryBarrier barrier{
      .sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER,
      .pNext = nullptr,
      .srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT,
      .dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
  };
  vkCmdPipelineBarrier(cmd_buffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                       VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &barrier, 0, nullptr, 0,
                       nullptr);
}

struct PushConstantsSteinhardt {
  std::array<float, 4> box_min;
  std::array<float, 4> box_max;
  uint32_t max_neighbors;
  uint32_t num_atoms;
  float cutoff;
  uint32_t pad;
};

struct PushConstantsHaloApply {
  std::array<float, 4> active_halo_color;
  float min_threshold;
  float max_threshold;
  float halo_radius_ratio;
  float halo_intensity;
  uint32_t dim_non_matching;
  uint32_t hide_non_matching;
  uint32_t num_atoms;
  uint32_t pad;
};

struct PushConstantsHistogram {
  float min_val;
  float max_val;
  float threshold_min;
  float threshold_max;
  uint32_t num_atoms;
};

} // namespace

void VulkanComputeService::cleanup() noexcept {
  if (device_ != VK_NULL_HANDLE) {
    if (descriptor_pool_ != VK_NULL_HANDLE) {
      vkDestroyDescriptorPool(device_, descriptor_pool_, nullptr);
      descriptor_pool_ = VK_NULL_HANDLE;
    }
    if (command_pool_ != VK_NULL_HANDLE) {
      vkDestroyCommandPool(device_, command_pool_, nullptr);
      command_pool_ = VK_NULL_HANDLE;
    }
    pipelines_.cleanup();
    cell_head_.destroy();
    atom_next_.destroy();
    neighbor_counts_.destroy();
    neighbor_indices_.destroy();
    descriptor_values_.destroy();
    complex_harmonics_.destroy();
    histogram_ssbo_.destroy();
    histogram_readback_.destroy();
    descriptor_readback_.destroy();
    host_positions_buf_.destroy();
    scratch_instances_buf_.destroy();
    device_ = VK_NULL_HANDLE;
  }
  is_initialized_ = false;
  allocated_atom_capacity_ = 0;
  allocated_cell_capacity_ = 0;
}

bool VulkanComputeService::initialize(VkDevice device, VkPhysicalDevice physical_device,
                                      VkQueue compute_queue, uint32_t queue_family_index) {
  cleanup();
  device_ = device;
  physical_device_ = physical_device;
  queue_ = compute_queue;
  queue_family_ = queue_family_index;

  if (!pipelines_.initialize(device_)) {
    cleanup();
    return false;
  }

  const VkCommandPoolCreateInfo pool_info{
      .sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO,
      .pNext = nullptr,
      .flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT,
      .queueFamilyIndex = queue_family_,
  };

  if (vkCreateCommandPool(device_, &pool_info, nullptr, &command_pool_) != VK_SUCCESS) {
    cleanup();
    return false;
  }

  const std::array<VkDescriptorPoolSize, 1> pool_sizes{
      VkDescriptorPoolSize{
          .type = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
          .descriptorCount = 128,
      },
  };

  const VkDescriptorPoolCreateInfo desc_pool_info{
      .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO,
      .pNext = nullptr,
      .flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT,
      .maxSets = 32,
      .poolSizeCount = static_cast<uint32_t>(pool_sizes.size()),
      .pPoolSizes = pool_sizes.data(),
  };

  if (vkCreateDescriptorPool(device_, &desc_pool_info, nullptr, &descriptor_pool_) != VK_SUCCESS) {
    cleanup();
    return false;
  }

  is_initialized_ = true;
  return true;
}

void VulkanComputeService::ensureBuffers(uint32_t num_atoms) {
  if (num_atoms <= allocated_atom_capacity_) {
    return;
  }
  allocated_atom_capacity_ = std::max(num_atoms * 2U, 4096U);

  const BufferAllocationDesc device_desc{
      .size = 0,
      .usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
      .properties = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
  };

  auto desc = device_desc;
  desc.size = allocated_atom_capacity_ * sizeof(uint32_t);
  atom_next_.ensureCapacity(device_, physical_device_, desc);
  neighbor_counts_.ensureCapacity(device_, physical_device_, desc);

  desc.size =
      static_cast<VkDeviceSize>(allocated_atom_capacity_) * MAX_NEIGHBORS * sizeof(uint32_t);
  neighbor_indices_.ensureCapacity(device_, physical_device_, desc);

  desc.size = static_cast<VkDeviceSize>(allocated_atom_capacity_) * 4U * sizeof(float);
  descriptor_values_.ensureCapacity(device_, physical_device_, desc);

  desc.size = static_cast<VkDeviceSize>(allocated_atom_capacity_) * 13U * 2U * sizeof(float);
  complex_harmonics_.ensureCapacity(device_, physical_device_, desc);

  desc.size = (100U + 2U) * sizeof(uint32_t);
  desc.usage = static_cast<VkBufferUsageFlags>(VK_BUFFER_USAGE_STORAGE_BUFFER_BIT) |
               static_cast<VkBufferUsageFlags>(VK_BUFFER_USAGE_TRANSFER_DST_BIT) |
               static_cast<VkBufferUsageFlags>(VK_BUFFER_USAGE_TRANSFER_SRC_BIT);
  histogram_ssbo_.ensureCapacity(device_, physical_device_, desc);

  desc.usage = VK_BUFFER_USAGE_TRANSFER_DST_BIT;
  desc.properties = static_cast<VkMemoryPropertyFlags>(VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT) |
                    static_cast<VkMemoryPropertyFlags>(VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
  histogram_readback_.ensureCapacity(device_, physical_device_, desc);
}

void VulkanComputeService::ensureCellGridBuffers(uint32_t total_cells) {
  if (total_cells <= allocated_cell_capacity_) {
    return;
  }
  allocated_cell_capacity_ = std::max(total_cells, 4096U);
  const BufferAllocationDesc cell_desc{
      .size = static_cast<VkDeviceSize>(allocated_cell_capacity_) * sizeof(uint32_t),
      .usage = static_cast<VkBufferUsageFlags>(VK_BUFFER_USAGE_STORAGE_BUFFER_BIT) |
               static_cast<VkBufferUsageFlags>(VK_BUFFER_USAGE_TRANSFER_DST_BIT),
      .properties = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
  };
  cell_head_.ensureCapacity(device_, physical_device_, cell_desc);
}

void VulkanComputeService::recordCellListBuildPass(
    VkCommandBuffer cmd_buffer, VkBuffer atom_positions_buffer, uint32_t num_atoms,
    const std::array<float, 4> &box_min, const std::array<float, 4> &box_max,
    const std::array<int32_t, 4> &grid_dims, uint32_t total_cells) {
  VkDescriptorSetLayout layout = pipelines_.descriptorSetLayout(ComputeKernelType::CellListBuild);
  const VkDescriptorSetAllocateInfo ds_alloc{
      .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
      .pNext = nullptr,
      .descriptorPool = descriptor_pool_,
      .descriptorSetCount = 1,
      .pSetLayouts = &layout,
  };
  VkDescriptorSet descriptor_set = VK_NULL_HANDLE;
  vkAllocateDescriptorSets(device_, &ds_alloc, &descriptor_set);

  const VkDescriptorBufferInfo desc_buf0{
      .buffer = atom_positions_buffer,
      .offset = 0,
      .range = static_cast<VkDeviceSize>(num_atoms) * sizeof(float) * 4U,
  };
  const VkDescriptorBufferInfo desc_buf1{
      .buffer = cell_head_.buffer(),
      .offset = 0,
      .range = static_cast<VkDeviceSize>(total_cells) * sizeof(uint32_t),
  };
  const VkDescriptorBufferInfo desc_buf2{
      .buffer = atom_next_.buffer(),
      .offset = 0,
      .range = static_cast<VkDeviceSize>(num_atoms) * sizeof(uint32_t),
  };

  const std::array<VkWriteDescriptorSet, 3> writes{
      VkWriteDescriptorSet{
          .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
          .pNext = nullptr,
          .dstSet = descriptor_set,
          .dstBinding = 0,
          .dstArrayElement = 0,
          .descriptorCount = 1,
          .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
          .pImageInfo = nullptr,
          .pBufferInfo = &desc_buf0,
          .pTexelBufferView = nullptr,
      },
      VkWriteDescriptorSet{
          .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
          .pNext = nullptr,
          .dstSet = descriptor_set,
          .dstBinding = 1,
          .dstArrayElement = 0,
          .descriptorCount = 1,
          .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
          .pImageInfo = nullptr,
          .pBufferInfo = &desc_buf1,
          .pTexelBufferView = nullptr,
      },
      VkWriteDescriptorSet{
          .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
          .pNext = nullptr,
          .dstSet = descriptor_set,
          .dstBinding = 2,
          .dstArrayElement = 0,
          .descriptorCount = 1,
          .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
          .pImageInfo = nullptr,
          .pBufferInfo = &desc_buf2,
          .pTexelBufferView = nullptr,
      },
  };
  vkUpdateDescriptorSets(device_, static_cast<uint32_t>(writes.size()), writes.data(), 0, nullptr);

  vkCmdBindPipeline(cmd_buffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                    pipelines_.pipeline(ComputeKernelType::CellListBuild));
  vkCmdBindDescriptorSets(cmd_buffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                          pipelines_.pipelineLayout(ComputeKernelType::CellListBuild), 0, 1,
                          &descriptor_set, 0, nullptr);

  const PushConstantsCellListBuild push_constants{
      .box_min = box_min,
      .box_max = box_max,
      .grid_dims = grid_dims,
      .num_atoms = num_atoms,
      .pad0 = 0,
      .pad1 = 0,
      .pad2 = 0,
  };
  vkCmdPushConstants(cmd_buffer, pipelines_.pipelineLayout(ComputeKernelType::CellListBuild),
                     VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(push_constants), &push_constants);

  const uint32_t groups = (num_atoms + 255U) / 256U;
  vkCmdDispatch(cmd_buffer, groups, 1, 1);
}

void VulkanComputeService::recordNeighborGridPass(
    VkCommandBuffer cmd_buffer, VkBuffer atom_positions_buffer, uint32_t num_atoms,
    const std::array<float, 4> &box_min, const std::array<float, 4> &box_max,
    const std::array<int32_t, 4> &grid_dims, float cutoff, uint32_t total_cells) {
  VkDescriptorSetLayout layout = pipelines_.descriptorSetLayout(ComputeKernelType::NeighborGrid);
  const VkDescriptorSetAllocateInfo ds_alloc{
      .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
      .pNext = nullptr,
      .descriptorPool = descriptor_pool_,
      .descriptorSetCount = 1,
      .pSetLayouts = &layout,
  };
  VkDescriptorSet descriptor_set = VK_NULL_HANDLE;
  vkAllocateDescriptorSets(device_, &ds_alloc, &descriptor_set);

  const VkDescriptorBufferInfo desc_buf0{
      .buffer = atom_positions_buffer,
      .offset = 0,
      .range = static_cast<VkDeviceSize>(num_atoms) * sizeof(float) * 4U,
  };
  const VkDescriptorBufferInfo desc_buf1{
      .buffer = neighbor_counts_.buffer(),
      .offset = 0,
      .range = static_cast<VkDeviceSize>(num_atoms) * sizeof(uint32_t),
  };
  const VkDescriptorBufferInfo desc_buf2{
      .buffer = neighbor_indices_.buffer(),
      .offset = 0,
      .range = static_cast<VkDeviceSize>(num_atoms) * MAX_NEIGHBORS * sizeof(uint32_t),
  };
  const VkDescriptorBufferInfo desc_buf3{
      .buffer = cell_head_.buffer(),
      .offset = 0,
      .range = static_cast<VkDeviceSize>(total_cells) * sizeof(uint32_t),
  };
  const VkDescriptorBufferInfo desc_buf4{
      .buffer = atom_next_.buffer(),
      .offset = 0,
      .range = static_cast<VkDeviceSize>(num_atoms) * sizeof(uint32_t),
  };

  const std::array<VkWriteDescriptorSet, 5> writes{
      VkWriteDescriptorSet{
          .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
          .pNext = nullptr,
          .dstSet = descriptor_set,
          .dstBinding = 0,
          .dstArrayElement = 0,
          .descriptorCount = 1,
          .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
          .pImageInfo = nullptr,
          .pBufferInfo = &desc_buf0,
          .pTexelBufferView = nullptr,
      },
      VkWriteDescriptorSet{
          .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
          .pNext = nullptr,
          .dstSet = descriptor_set,
          .dstBinding = 1,
          .dstArrayElement = 0,
          .descriptorCount = 1,
          .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
          .pImageInfo = nullptr,
          .pBufferInfo = &desc_buf1,
          .pTexelBufferView = nullptr,
      },
      VkWriteDescriptorSet{
          .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
          .pNext = nullptr,
          .dstSet = descriptor_set,
          .dstBinding = 2,
          .dstArrayElement = 0,
          .descriptorCount = 1,
          .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
          .pImageInfo = nullptr,
          .pBufferInfo = &desc_buf2,
          .pTexelBufferView = nullptr,
      },
      VkWriteDescriptorSet{
          .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
          .pNext = nullptr,
          .dstSet = descriptor_set,
          .dstBinding = 3,
          .dstArrayElement = 0,
          .descriptorCount = 1,
          .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
          .pImageInfo = nullptr,
          .pBufferInfo = &desc_buf3,
          .pTexelBufferView = nullptr,
      },
      VkWriteDescriptorSet{
          .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
          .pNext = nullptr,
          .dstSet = descriptor_set,
          .dstBinding = 4,
          .dstArrayElement = 0,
          .descriptorCount = 1,
          .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
          .pImageInfo = nullptr,
          .pBufferInfo = &desc_buf4,
          .pTexelBufferView = nullptr,
      },
  };
  vkUpdateDescriptorSets(device_, static_cast<uint32_t>(writes.size()), writes.data(), 0, nullptr);

  vkCmdBindPipeline(cmd_buffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                    pipelines_.pipeline(ComputeKernelType::NeighborGrid));
  vkCmdBindDescriptorSets(cmd_buffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                          pipelines_.pipelineLayout(ComputeKernelType::NeighborGrid), 0, 1,
                          &descriptor_set, 0, nullptr);

  const PushConstantsNeighborGrid push_constants{
      .box_min = box_min,
      .box_max = box_max,
      .grid_dims = grid_dims,
      .cutoff = cutoff,
      .max_neighbors = MAX_NEIGHBORS,
      .num_atoms = num_atoms,
      .pad = 0,
  };
  vkCmdPushConstants(cmd_buffer, pipelines_.pipelineLayout(ComputeKernelType::NeighborGrid),
                     VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(push_constants), &push_constants);

  const uint32_t groups = (num_atoms + 255U) / 256U;
  vkCmdDispatch(cmd_buffer, groups, 1, 1);
}

bool VulkanComputeService::evaluateOrderParametersAndApplyHalos(
    VkBuffer atom_positions_buffer, VkBuffer atom_instance_buffer, uint32_t num_atoms,
    const std::array<float, 4> &box_min, const std::array<float, 4> &box_max, float cutoff,
    const GpuHaloSettings &settings, GpuHistogramOutput &out_histogram,
    std::vector<float> *out_descriptors) {
  if (!is_initialized_ || num_atoms == 0) {
    return false;
  }

  const std::array<int32_t, 4> grid_dims = computeGridDims(box_min, box_max, cutoff);
  const uint32_t total_cells = static_cast<uint32_t>(grid_dims[0]) *
                               static_cast<uint32_t>(grid_dims[1]) *
                               static_cast<uint32_t>(grid_dims[2]);

  ensureBuffers(num_atoms);
  ensureCellGridBuffers(total_cells);

  const VkCommandBufferAllocateInfo alloc_info{
      .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,
      .pNext = nullptr,
      .commandPool = command_pool_,
      .level = VK_COMMAND_BUFFER_LEVEL_PRIMARY,
      .commandBufferCount = 1,
  };

  VkCommandBuffer cmd_buffer = VK_NULL_HANDLE;
  if (vkAllocateCommandBuffers(device_, &alloc_info, &cmd_buffer) != VK_SUCCESS) {
    return false;
  }

  const VkCommandBufferBeginInfo begin_info{
      .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO,
      .pNext = nullptr,
      .flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT,
      .pInheritanceInfo = nullptr,
  };

  vkBeginCommandBuffer(cmd_buffer, &begin_info);

  vkCmdFillBuffer(cmd_buffer, histogram_ssbo_.buffer(), 0, (100U + 2U) * sizeof(uint32_t), 0);
  vkCmdFillBuffer(cmd_buffer, cell_head_.buffer(), 0,
                  static_cast<VkDeviceSize>(total_cells) * sizeof(uint32_t), EMPTY_CELL_SENTINEL);

  const VkMemoryBarrier barrier_clear{
      .sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER,
      .pNext = nullptr,
      .srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT,
      .dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
  };
  vkCmdPipelineBarrier(cmd_buffer, VK_PIPELINE_STAGE_TRANSFER_BIT,
                       VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &barrier_clear, 0, nullptr, 0,
                       nullptr);

  // Pass 1a: O(N) 3D Uniform Spatial Cell-List Build
  recordCellListBuildPass(cmd_buffer, atom_positions_buffer, num_atoms, box_min, box_max, grid_dims,
                          total_cells);
  emitComputeBarrier(cmd_buffer);

  // Pass 1b: O(N) 27-Cell Neighbor Grid Traversal
  recordNeighborGridPass(cmd_buffer, atom_positions_buffer, num_atoms, box_min, box_max, grid_dims,
                         cutoff, total_cells);
  emitComputeBarrier(cmd_buffer);

  // Pass 2: Steinhardt Order Parameters
  {
    VkDescriptorSetLayout layout =
        pipelines_.descriptorSetLayout(ComputeKernelType::SteinhardtOrder);
    const VkDescriptorSetAllocateInfo ds_alloc{
        .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
        .pNext = nullptr,
        .descriptorPool = descriptor_pool_,
        .descriptorSetCount = 1,
        .pSetLayouts = &layout,
    };
    VkDescriptorSet descriptor_set = VK_NULL_HANDLE;
    vkAllocateDescriptorSets(device_, &ds_alloc, &descriptor_set);

    const VkDescriptorBufferInfo desc_buf0{
        .buffer = atom_positions_buffer,
        .offset = 0,
        .range = num_atoms * sizeof(float) * 4U,
    };
    const VkDescriptorBufferInfo desc_buf1{
        .buffer = neighbor_counts_.buffer(),
        .offset = 0,
        .range = num_atoms * sizeof(uint32_t),
    };
    const VkDescriptorBufferInfo desc_buf2{
        .buffer = neighbor_indices_.buffer(),
        .offset = 0,
        .range = static_cast<VkDeviceSize>(num_atoms) * MAX_NEIGHBORS * sizeof(uint32_t),
    };
    const VkDescriptorBufferInfo desc_buf3{
        .buffer = descriptor_values_.buffer(),
        .offset = 0,
        .range = static_cast<VkDeviceSize>(num_atoms) * 4U * sizeof(float),
    };
    const VkDescriptorBufferInfo desc_buf4{
        .buffer = complex_harmonics_.buffer(),
        .offset = 0,
        .range = static_cast<VkDeviceSize>(num_atoms) * 13U * 2U * sizeof(float),
    };

    const std::array<VkWriteDescriptorSet, 5> writes{
        VkWriteDescriptorSet{
            .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
            .pNext = nullptr,
            .dstSet = descriptor_set,
            .dstBinding = 0,
            .dstArrayElement = 0,
            .descriptorCount = 1,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .pImageInfo = nullptr,
            .pBufferInfo = &desc_buf0,
            .pTexelBufferView = nullptr,
        },
        VkWriteDescriptorSet{
            .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
            .pNext = nullptr,
            .dstSet = descriptor_set,
            .dstBinding = 1,
            .dstArrayElement = 0,
            .descriptorCount = 1,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .pImageInfo = nullptr,
            .pBufferInfo = &desc_buf1,
            .pTexelBufferView = nullptr,
        },
        VkWriteDescriptorSet{
            .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
            .pNext = nullptr,
            .dstSet = descriptor_set,
            .dstBinding = 2,
            .dstArrayElement = 0,
            .descriptorCount = 1,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .pImageInfo = nullptr,
            .pBufferInfo = &desc_buf2,
            .pTexelBufferView = nullptr,
        },
        VkWriteDescriptorSet{
            .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
            .pNext = nullptr,
            .dstSet = descriptor_set,
            .dstBinding = 3,
            .dstArrayElement = 0,
            .descriptorCount = 1,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .pImageInfo = nullptr,
            .pBufferInfo = &desc_buf3,
            .pTexelBufferView = nullptr,
        },
        VkWriteDescriptorSet{
            .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
            .pNext = nullptr,
            .dstSet = descriptor_set,
            .dstBinding = 4,
            .dstArrayElement = 0,
            .descriptorCount = 1,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .pImageInfo = nullptr,
            .pBufferInfo = &desc_buf4,
            .pTexelBufferView = nullptr,
        },
    };
    vkUpdateDescriptorSets(device_, static_cast<uint32_t>(writes.size()), writes.data(), 0,
                           nullptr);

    vkCmdBindPipeline(cmd_buffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                      pipelines_.pipeline(ComputeKernelType::SteinhardtOrder));
    vkCmdBindDescriptorSets(cmd_buffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                            pipelines_.pipelineLayout(ComputeKernelType::SteinhardtOrder), 0, 1,
                            &descriptor_set, 0, nullptr);

    const PushConstantsSteinhardt push_constants{
        .box_min = box_min,
        .box_max = box_max,
        .max_neighbors = MAX_NEIGHBORS,
        .num_atoms = num_atoms,
        .cutoff = cutoff,
        .pad = 0,
    };

    vkCmdPushConstants(cmd_buffer, pipelines_.pipelineLayout(ComputeKernelType::SteinhardtOrder),
                       VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(push_constants), &push_constants);

    const uint32_t groups = (num_atoms + 255U) / 256U;
    vkCmdDispatch(cmd_buffer, groups, 1, 1);
  }

  emitComputeBarrier(cmd_buffer);

  // Pass 3: In-VRAM Halo Application
  {
    VkDescriptorSetLayout layout = pipelines_.descriptorSetLayout(ComputeKernelType::HaloApply);
    const VkDescriptorSetAllocateInfo ds_alloc{
        .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
        .pNext = nullptr,
        .descriptorPool = descriptor_pool_,
        .descriptorSetCount = 1,
        .pSetLayouts = &layout,
    };
    VkDescriptorSet descriptor_set = VK_NULL_HANDLE;
    vkAllocateDescriptorSets(device_, &ds_alloc, &descriptor_set);

    const VkDescriptorBufferInfo desc_buf0{
        .buffer = descriptor_values_.buffer(),
        .offset = 0,
        .range = num_atoms * sizeof(float),
    };
    const VkDescriptorBufferInfo desc_buf1{
        .buffer = atom_instance_buffer,
        .offset = 0,
        .range = num_atoms * (14U * sizeof(float)),
    };

    const std::array<VkWriteDescriptorSet, 2> writes{
        VkWriteDescriptorSet{
            .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
            .pNext = nullptr,
            .dstSet = descriptor_set,
            .dstBinding = 0,
            .dstArrayElement = 0,
            .descriptorCount = 1,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .pImageInfo = nullptr,
            .pBufferInfo = &desc_buf0,
            .pTexelBufferView = nullptr,
        },
        VkWriteDescriptorSet{
            .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
            .pNext = nullptr,
            .dstSet = descriptor_set,
            .dstBinding = 1,
            .dstArrayElement = 0,
            .descriptorCount = 1,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .pImageInfo = nullptr,
            .pBufferInfo = &desc_buf1,
            .pTexelBufferView = nullptr,
        },
    };
    vkUpdateDescriptorSets(device_, static_cast<uint32_t>(writes.size()), writes.data(), 0,
                           nullptr);

    vkCmdBindPipeline(cmd_buffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                      pipelines_.pipeline(ComputeKernelType::HaloApply));
    vkCmdBindDescriptorSets(cmd_buffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                            pipelines_.pipelineLayout(ComputeKernelType::HaloApply), 0, 1,
                            &descriptor_set, 0, nullptr);

    const PushConstantsHaloApply push_constants{
        .active_halo_color = settings.halo_color,
        .min_threshold = settings.min_threshold,
        .max_threshold = settings.max_threshold,
        .halo_radius_ratio = settings.halo_radius_ratio,
        .halo_intensity = settings.halo_intensity,
        .dim_non_matching = settings.dim_non_matching ? 1U : 0U,
        .hide_non_matching = settings.hide_non_matching ? 1U : 0U,
        .num_atoms = num_atoms,
        .pad = 0,
    };

    vkCmdPushConstants(cmd_buffer, pipelines_.pipelineLayout(ComputeKernelType::HaloApply),
                       VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(push_constants), &push_constants);

    const uint32_t groups = (num_atoms + 255U) / 256U;
    vkCmdDispatch(cmd_buffer, groups, 1, 1);
  }

  // Pass 4: Histogram Atomic Reduction
  {
    VkDescriptorSetLayout layout =
        pipelines_.descriptorSetLayout(ComputeKernelType::HistogramReduce);
    const VkDescriptorSetAllocateInfo ds_alloc{
        .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
        .pNext = nullptr,
        .descriptorPool = descriptor_pool_,
        .descriptorSetCount = 1,
        .pSetLayouts = &layout,
    };
    VkDescriptorSet descriptor_set = VK_NULL_HANDLE;
    vkAllocateDescriptorSets(device_, &ds_alloc, &descriptor_set);

    const VkDescriptorBufferInfo desc_buf0{
        .buffer = descriptor_values_.buffer(),
        .offset = 0,
        .range = num_atoms * sizeof(float),
    };
    const VkDescriptorBufferInfo desc_buf1{
        .buffer = histogram_ssbo_.buffer(),
        .offset = 0,
        .range = (100U + 2U) * sizeof(uint32_t),
    };

    const std::array<VkWriteDescriptorSet, 2> writes{
        VkWriteDescriptorSet{
            .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
            .pNext = nullptr,
            .dstSet = descriptor_set,
            .dstBinding = 0,
            .dstArrayElement = 0,
            .descriptorCount = 1,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .pImageInfo = nullptr,
            .pBufferInfo = &desc_buf0,
            .pTexelBufferView = nullptr,
        },
        VkWriteDescriptorSet{
            .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
            .pNext = nullptr,
            .dstSet = descriptor_set,
            .dstBinding = 1,
            .dstArrayElement = 0,
            .descriptorCount = 1,
            .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            .pImageInfo = nullptr,
            .pBufferInfo = &desc_buf1,
            .pTexelBufferView = nullptr,
        },
    };
    vkUpdateDescriptorSets(device_, static_cast<uint32_t>(writes.size()), writes.data(), 0,
                           nullptr);

    vkCmdBindPipeline(cmd_buffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                      pipelines_.pipeline(ComputeKernelType::HistogramReduce));
    vkCmdBindDescriptorSets(cmd_buffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                            pipelines_.pipelineLayout(ComputeKernelType::HistogramReduce), 0, 1,
                            &descriptor_set, 0, nullptr);

    const PushConstantsHistogram push_constants{
        .min_val = 0.0F,
        .max_val = 1.0F,
        .threshold_min = settings.min_threshold,
        .threshold_max = settings.max_threshold,
        .num_atoms = num_atoms,
    };

    vkCmdPushConstants(cmd_buffer, pipelines_.pipelineLayout(ComputeKernelType::HistogramReduce),
                       VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(push_constants), &push_constants);

    const uint32_t groups = (num_atoms + 255U) / 256U;
    vkCmdDispatch(cmd_buffer, groups, 1, 1);
  }

  const VkMemoryBarrier barrier_readback{
      .sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER,
      .pNext = nullptr,
      .srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT,
      .dstAccessMask = VK_ACCESS_TRANSFER_READ_BIT,
  };
  vkCmdPipelineBarrier(cmd_buffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                       VK_PIPELINE_STAGE_TRANSFER_BIT, 0, 1, &barrier_readback, 0, nullptr, 0,
                       nullptr);

  const VkBufferCopy copy_region{
      .srcOffset = 0,
      .dstOffset = 0,
      .size = (100U + 2U) * sizeof(uint32_t),
  };
  vkCmdCopyBuffer(cmd_buffer, histogram_ssbo_.buffer(), histogram_readback_.buffer(), 1,
                  &copy_region);

  if (out_descriptors != nullptr) {
    const BufferAllocationDesc readback_desc{
        .size = static_cast<VkDeviceSize>(num_atoms) * sizeof(float),
        .usage = VK_BUFFER_USAGE_TRANSFER_DST_BIT,
        .properties = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
    };
    descriptor_readback_.ensureCapacity(device_, physical_device_, readback_desc);
    const VkBufferCopy desc_copy{
        .srcOffset = 0,
        .dstOffset = 0,
        .size = static_cast<VkDeviceSize>(num_atoms) * sizeof(float),
    };
    vkCmdCopyBuffer(cmd_buffer, descriptor_values_.buffer(), descriptor_readback_.buffer(), 1,
                    &desc_copy);
  }

  vkEndCommandBuffer(cmd_buffer);

  const VkSubmitInfo submit_info{
      .sType = VK_STRUCTURE_TYPE_SUBMIT_INFO,
      .pNext = nullptr,
      .waitSemaphoreCount = 0,
      .pWaitSemaphores = nullptr,
      .pWaitDstStageMask = nullptr,
      .commandBufferCount = 1,
      .pCommandBuffers = &cmd_buffer,
      .signalSemaphoreCount = 0,
      .pSignalSemaphores = nullptr,
  };

  vkQueueSubmit(queue_, 1, &submit_info, VK_NULL_HANDLE);
  vkQueueWaitIdle(queue_);

  vkFreeCommandBuffers(device_, command_pool_, 1, &cmd_buffer);
  vkResetDescriptorPool(device_, descriptor_pool_, 0);

  const auto *readback_ptr = reinterpret_cast<const uint32_t *>(histogram_readback_.map());
  if (readback_ptr != nullptr) {
    std::memcpy(out_histogram.bins.data(), readback_ptr, 100 * sizeof(uint32_t));
    out_histogram.matching_count = readback_ptr[100];
    out_histogram.total_count = readback_ptr[101];
    histogram_readback_.unmap();
  }

  if (out_descriptors != nullptr) {
    out_descriptors->resize(num_atoms);
    const auto *desc_ptr = reinterpret_cast<const float *>(descriptor_readback_.map());
    if (desc_ptr != nullptr) {
      std::memcpy(out_descriptors->data(), desc_ptr,
                  static_cast<std::size_t>(num_atoms) * sizeof(float));
      descriptor_readback_.unmap();
    }
  }

  return true;
}

bool VulkanComputeService::computeDescriptorValues(
    const std::vector<std::array<float, 4>> &positions, const std::array<float, 4> &box_min,
    const std::array<float, 4> &box_max, float cutoff, app::PhaseDescriptor descriptor,
    std::vector<float> &out_values) {
  if (!is_initialized_ || positions.empty()) {
    return false;
  }
  const auto num_atoms = static_cast<uint32_t>(positions.size());
  const BufferAllocationDesc pos_desc{
      .size = static_cast<VkDeviceSize>(num_atoms) * sizeof(std::array<float, 4>),
      .usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
      .properties = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
  };
  host_positions_buf_.ensureCapacity(device_, physical_device_, pos_desc);
  host_positions_buf_.upload(std::as_bytes(std::span{positions.data(), positions.size()}));

  constexpr VkDeviceSize K_ATOM_INSTANCE_STRIDE = 56;
  const BufferAllocationDesc inst_desc{
      .size = static_cast<VkDeviceSize>(num_atoms) * K_ATOM_INSTANCE_STRIDE,
      .usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
      .properties = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
  };
  scratch_instances_buf_.ensureCapacity(device_, physical_device_, inst_desc);

  GpuHaloSettings settings{};
  settings.descriptor = descriptor;
  GpuHistogramOutput dummy_hist{};
  return evaluateOrderParametersAndApplyHalos(host_positions_buf_.buffer(),
                                              scratch_instances_buf_.buffer(), num_atoms, box_min,
                                              box_max, cutoff, settings, dummy_hist, &out_values);
}

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d::vulkan
