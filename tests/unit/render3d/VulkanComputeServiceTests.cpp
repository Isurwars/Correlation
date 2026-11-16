/**
 * @file VulkanComputeServiceTests.cpp
 * @brief Integration tests for VulkanComputeService GPU order parameter evaluation and halos.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/vulkan/VulkanContext.hpp"
#include "render3d/vulkan/VulkanDevice.hpp"
#include "render3d/vulkan/VulkanBuffer.hpp"
#include "render3d/vulkan/compute/VulkanComputeService.hpp"
#include <gtest/gtest.h>
#include <vector>

#ifdef CORRELATION_USE_VULKAN

namespace correlation::render3d::vulkan {

namespace {

class VulkanComputeServiceTests : public ::testing::Test {
protected:
  void SetUp() override {
    if (!VulkanContext::isSupported()) {
      GTEST_SKIP() << "Vulkan driver or compatible physical device not available";
    }
  }
};

TEST_F(VulkanComputeServiceTests, InitializeAndExecuteSteinhardtOnFCC) {
  VulkanDevice device;
  ASSERT_TRUE(device.initialize());

  VulkanComputeService service;
  ASSERT_TRUE(service.initialize(device.device(), device.physicalDevice(),
                                 device.graphicsQueue(), device.queueFamilyIndex()));

  // 4 atoms forming a tetrahedral / close-packed cluster
  const std::vector<float> atom_positions = {
      0.0F, 0.0F, 0.0F, 0.0F,
      1.5F, 1.5F, 0.0F, 0.0F,
      1.5F, 0.0F, 1.5F, 0.0F,
      0.0F, 1.5F, 1.5F, 0.0F,
  };

  VulkanBuffer pos_buffer;
  BufferAllocationDesc pos_desc{
      .size = atom_positions.size() * sizeof(float),
      .usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
      .properties = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
  };
  pos_buffer.create(device.device(), device.physicalDevice(), pos_desc);
  ASSERT_TRUE(pos_buffer.isValid());
  pos_buffer.upload(std::as_bytes(std::span{atom_positions}));

  VulkanBuffer inst_buffer;
  BufferAllocationDesc inst_desc{
      .size = 4 * 14 * sizeof(float),
      .usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
      .properties = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
  };
  inst_buffer.create(device.device(), device.physicalDevice(), inst_desc);
  ASSERT_TRUE(inst_buffer.isValid());

  GpuHaloSettings settings{
      .descriptor = app::PhaseDescriptor::Q6,
      .min_threshold = 0.2F,
      .max_threshold = 0.8F,
      .halo_radius_ratio = 1.4F,
      .halo_intensity = 0.9F,
      .halo_color = {1.0F, 0.84F, 0.0F, 1.0F},
      .dim_non_matching = false,
      .hide_non_matching = false,
  };

  GpuHistogramOutput hist_out;
  const std::array<float, 4> box_min{-10.0F, -10.0F, -10.0F, 0.0F};
  const std::array<float, 4> box_max{10.0F, 10.0F, 10.0F, 0.0F};

  bool success = service.evaluateOrderParametersAndApplyHalos(
      pos_buffer.buffer(), inst_buffer.buffer(), 4, box_min, box_max, 3.2F, settings, hist_out);
  EXPECT_TRUE(success);

  EXPECT_EQ(hist_out.total_count, 4);

  pos_buffer.destroy();
  inst_buffer.destroy();
  service.cleanup();
  device.destroy();
}

} // namespace
} // namespace correlation::render3d::vulkan

#endif // CORRELATION_USE_VULKAN
