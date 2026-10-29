/**
 * @file VulkanPipelineManagerTests.cpp
 * @brief Unit tests for VulkanPipelineManager graphics pipelines, layouts, and shader modules.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/vulkan/VulkanContext.hpp"
#include "render3d/vulkan/VulkanDevice.hpp"
#include "render3d/vulkan/VulkanOffscreenTarget.hpp"
#include "render3d/vulkan/VulkanPipelineManager.hpp"
#include <gtest/gtest.h>

#ifdef CORRELATION_USE_VULKAN

namespace correlation::render3d::vulkan {
namespace {

class VulkanPipelineManagerTests : public ::testing::Test {
protected:
  void SetUp() override {
    if (!VulkanContext::isSupported()) {
      GTEST_SKIP() << "Vulkan driver or compatible physical device not available";
    }
    ASSERT_TRUE(device_.initialize());
    target_.create(device_.device(), device_.physicalDevice(), 640, 480);
    ASSERT_TRUE(target_.isValid());
  }

  void TearDown() override {
    target_.destroy();
    device_.destroy();
  }

  [[nodiscard]] VulkanDevice &device() noexcept { return device_; }
  [[nodiscard]] VulkanOffscreenTarget &target() noexcept { return target_; }

private:
  VulkanDevice device_;
  VulkanOffscreenTarget target_;
};

TEST_F(VulkanPipelineManagerTests, InitializationAndPipelineHandles) {
  VulkanPipelineManager manager;
  EXPECT_FALSE(manager.isValid());
  EXPECT_EQ(manager.spherePipeline(), VK_NULL_HANDLE);
  EXPECT_EQ(manager.cylinderPipeline(), VK_NULL_HANDLE);
  EXPECT_EQ(manager.wireframePipeline(), VK_NULL_HANDLE);
  EXPECT_EQ(manager.pipelineLayout(), VK_NULL_HANDLE);
  EXPECT_EQ(manager.descriptorSetLayout(), VK_NULL_HANDLE);

  EXPECT_TRUE(manager.initialize(device().device(), target().renderPass()));
  EXPECT_TRUE(manager.isValid());
  EXPECT_NE(manager.spherePipeline(), VK_NULL_HANDLE);
  EXPECT_NE(manager.cylinderPipeline(), VK_NULL_HANDLE);
  EXPECT_NE(manager.wireframePipeline(), VK_NULL_HANDLE);
  EXPECT_NE(manager.pipelineLayout(), VK_NULL_HANDLE);
  EXPECT_NE(manager.descriptorSetLayout(), VK_NULL_HANDLE);

  manager.destroy();
  EXPECT_FALSE(manager.isValid());
  EXPECT_EQ(manager.spherePipeline(), VK_NULL_HANDLE);
}

TEST_F(VulkanPipelineManagerTests, MoveSemantics) {
  VulkanPipelineManager manager1;
  ASSERT_TRUE(manager1.initialize(device().device(), target().renderPass()));
  VkPipeline orig_sphere = manager1.spherePipeline();
  EXPECT_NE(orig_sphere, VK_NULL_HANDLE);

  // Move constructor
  VulkanPipelineManager manager2(std::move(manager1));
  EXPECT_TRUE(manager2.isValid());
  EXPECT_EQ(manager2.spherePipeline(), orig_sphere);

  // Move assignment
  VulkanPipelineManager manager3;
  manager3 = std::move(manager2);
  EXPECT_TRUE(manager3.isValid());
  EXPECT_EQ(manager3.spherePipeline(), orig_sphere);
}

} // namespace
} // namespace correlation::render3d::vulkan

#endif // CORRELATION_USE_VULKAN
