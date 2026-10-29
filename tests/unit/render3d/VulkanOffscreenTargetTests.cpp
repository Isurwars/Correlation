/**
 * @file VulkanOffscreenTargetTests.cpp
 * @brief Unit tests for headless Vulkan offscreen render target and staging readout.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/vulkan/VulkanContext.hpp"
#include "render3d/vulkan/VulkanDevice.hpp"
#include "render3d/vulkan/VulkanOffscreenTarget.hpp"
#include <gtest/gtest.h>

namespace correlation::render3d::vulkan {
namespace {

class VulkanOffscreenTargetTests : public ::testing::Test {
protected:
  void SetUp() override {
    if (!VulkanContext::isSupported()) {
      GTEST_SKIP() << "Vulkan driver or compatible physical device not available";
    }
    ASSERT_TRUE(device_.initialize());
  }

  [[nodiscard]] VulkanDevice &device() noexcept { return device_; }

private:
  VulkanDevice device_;
};

TEST_F(VulkanOffscreenTargetTests, InitializationAndProperties) {
  VulkanOffscreenTarget target;
  EXPECT_FALSE(target.isValid());
  EXPECT_EQ(target.renderPass(), VK_NULL_HANDLE);
  EXPECT_EQ(target.framebuffer(), VK_NULL_HANDLE);

  target.create(device().device(), device().physicalDevice(), 640, 480);
  EXPECT_TRUE(target.isValid());
  EXPECT_EQ(target.width(), 640);
  EXPECT_EQ(target.height(), 480);
  EXPECT_EQ(target.extent().width, 640);
  EXPECT_EQ(target.extent().height, 480);
  EXPECT_NE(target.renderPass(), VK_NULL_HANDLE);
  EXPECT_NE(target.framebuffer(), VK_NULL_HANDLE);
  EXPECT_NE(target.colorImage(), VK_NULL_HANDLE);
  EXPECT_NE(target.depthImage(), VK_NULL_HANDLE);

  const auto span = target.readoutStagingSpan(0);
  EXPECT_EQ(span.size(), 640 * 480 * 4);

  target.destroy();
  EXPECT_FALSE(target.isValid());
}

TEST_F(VulkanOffscreenTargetTests, DynamicResize) {
  VulkanOffscreenTarget target;
  target.create(device().device(), device().physicalDevice(), 400, 300);
  EXPECT_TRUE(target.isValid());
  EXPECT_EQ(target.width(), 400);
  EXPECT_EQ(target.height(), 300);

  target.resize(800, 600);
  EXPECT_TRUE(target.isValid());
  EXPECT_EQ(target.width(), 800);
  EXPECT_EQ(target.height(), 600);

  const auto span = target.readoutStagingSpan(0);
  EXPECT_EQ(span.size(), 800 * 600 * 4);

  // Resize to identical dimensions is a no-op
  auto *const fb_before = target.framebuffer();
  target.resize(800, 600);
  EXPECT_EQ(target.framebuffer(), fb_before);
}

TEST_F(VulkanOffscreenTargetTests, RenderPassClearAndStagingReadout) {
  VulkanOffscreenTarget target;
  target.create(device().device(), device().physicalDevice(), 64, 64);
  ASSERT_TRUE(target.isValid());

  VkCommandBuffer cmd = device().beginSingleTimeCommands();
  ASSERT_NE(cmd, VK_NULL_HANDLE);

  constexpr VkClearColorValue CLEAR_COLOR = {{0.25F, 0.5F, 0.75F, 1.0F}};
  target.beginRenderPass(cmd, CLEAR_COLOR);
  VulkanOffscreenTarget::endRenderPass(cmd);
  target.copyToStaging(cmd, 0);

  device().endSingleTimeCommands(cmd);

  const auto span = target.readoutStagingSpan(0);
  ASSERT_EQ(span.size(), 64 * 64 * 4);

  // Read first pixel (RGBA)
  // 0.25 * 255 = ~64, 0.5 * 255 = ~127, 0.75 * 255 = ~191, 1.0 * 255 = 255
  EXPECT_NEAR(span[0], 64, 2);
  EXPECT_NEAR(span[1], 127, 2);
  EXPECT_NEAR(span[2], 191, 2);
  EXPECT_EQ(span[3], 255);
}

TEST_F(VulkanOffscreenTargetTests, DoubleBufferedStagingReadout) {
  VulkanOffscreenTarget target;
  target.create(device().device(), device().physicalDevice(), 32, 32);
  ASSERT_TRUE(target.isValid());

  // Buffer 0 clear with green
  {
    VkCommandBuffer cmd = device().beginSingleTimeCommands();
    constexpr VkClearColorValue CLEAR_GREEN = {{0.0F, 1.0F, 0.0F, 1.0F}};
    target.beginRenderPass(cmd, CLEAR_GREEN);
    VulkanOffscreenTarget::endRenderPass(cmd);
    target.copyToStaging(cmd, 0);
    device().endSingleTimeCommands(cmd);
  }

  // Buffer 1 clear with red
  {
    VkCommandBuffer cmd = device().beginSingleTimeCommands();
    constexpr VkClearColorValue CLEAR_RED = {{1.0F, 0.0F, 0.0F, 1.0F}};
    target.beginRenderPass(cmd, CLEAR_RED);
    VulkanOffscreenTarget::endRenderPass(cmd);
    target.copyToStaging(cmd, 1);
    device().endSingleTimeCommands(cmd);
  }

  const auto span0 = target.readoutStagingSpan(0);
  const auto span1 = target.readoutStagingSpan(1);

  EXPECT_EQ(span0.size(), 32 * 32 * 4);
  EXPECT_EQ(span1.size(), 32 * 32 * 4);

  // Span 0 should be green
  EXPECT_EQ(span0[0], 0);
  EXPECT_EQ(span0[1], 255);
  EXPECT_EQ(span0[2], 0);
  EXPECT_EQ(span0[3], 255);

  // Span 1 should be red
  EXPECT_EQ(span1[0], 255);
  EXPECT_EQ(span1[1], 0);
  EXPECT_EQ(span1[2], 0);
  EXPECT_EQ(span1[3], 255);
}

} // namespace
} // namespace correlation::render3d::vulkan
