/**
 * @file VulkanComputeService.hpp
 * @brief Standalone service executing GPU-only order parameter evaluation and halo application.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#ifdef CORRELATION_USE_VULKAN
#include "render3d/vulkan/VulkanBuffer.hpp"
#include "render3d/vulkan/VulkanComputePipelineManager.hpp"
#include <vulkan/vulkan.h>
#endif

#include "app/viewmodel/PhaseSearchTypes.hpp"
#include <array>
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

namespace correlation::render3d::vulkan {

#ifdef CORRELATION_USE_VULKAN

struct GpuHaloSettings {
  app::PhaseDescriptor descriptor{app::PhaseDescriptor::Q6Bar};
  float min_threshold{0.46F};
  float max_threshold{0.58F};
  float halo_radius_ratio{2.70F};
  float halo_intensity{0.85F};
  std::array<float, 4> halo_color{1.0F, 0.84F, 0.0F, 1.0F};
  bool dim_non_matching{false};
  bool hide_non_matching{false};
};

struct GpuHistogramOutput {
  std::array<uint32_t, 100> bins{};
  uint32_t matching_count{0};
  uint32_t total_count{0};
};

class VulkanComputeService {
public:
  VulkanComputeService() = default;
  ~VulkanComputeService() noexcept { cleanup(); }

  VulkanComputeService(const VulkanComputeService &) = delete;
  VulkanComputeService &operator=(const VulkanComputeService &) = delete;

  [[nodiscard]] bool initialize(VkDevice device, VkPhysicalDevice physical_device,
                                VkQueue compute_queue, uint32_t queue_family_index);
  void cleanup() noexcept;

  /**
   * @brief Executes the GPU compute chain:
   * 1. Neighbor Grid construction
   * 2. Active order parameter calculation (Steinhardt, Tetrahedral, LocalEntropy, or TWF)
   * 3. Direct write to target atom instance buffer (halos, styles)
   * 4. Histogram atomic reduction
   */
  bool evaluateOrderParametersAndApplyHalos(VkBuffer atom_positions_buffer,
                                            VkBuffer atom_instance_buffer, uint32_t num_atoms,
                                            const std::array<float, 4> &box_min,
                                            const std::array<float, 4> &box_max, float cutoff,
                                            const GpuHaloSettings &settings,
                                            GpuHistogramOutput &out_histogram);

  [[nodiscard]] bool isValid() const noexcept { return is_initialized_; }

private:
  void ensureBuffers(uint32_t num_atoms);

  VkDevice device_{VK_NULL_HANDLE};
  VkPhysicalDevice physical_device_{VK_NULL_HANDLE};
  VkQueue queue_{VK_NULL_HANDLE};
  uint32_t queue_family_{0};
  VkCommandPool command_pool_{VK_NULL_HANDLE};
  VkDescriptorPool descriptor_pool_{VK_NULL_HANDLE};

  VulkanComputePipelineManager pipelines_;

  VulkanBuffer neighbor_counts_;
  VulkanBuffer neighbor_indices_;
  VulkanBuffer descriptor_values_;
  VulkanBuffer complex_harmonics_;
  VulkanBuffer histogram_ssbo_;
  VulkanBuffer histogram_readback_;

  uint32_t allocated_atom_capacity_{0};
  bool is_initialized_{false};
};

#endif // CORRELATION_USE_VULKAN

} // namespace correlation::render3d::vulkan
