/**
 * @file VulkanBufferManagerTests.cpp
 * @brief Unit tests for VulkanBufferManager allocations, instance buffers, and scene uploads.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/Camera3D.hpp"
#include "render3d/Scene3D.hpp"
#include "render3d/vulkan/VulkanBufferManager.hpp"
#include "render3d/vulkan/VulkanContext.hpp"
#include "render3d/vulkan/VulkanDevice.hpp"
#include "render3d/vulkan/VulkanOffscreenTarget.hpp"
#include "render3d/vulkan/VulkanPipelineManager.hpp"
#include <gtest/gtest.h>

#ifdef CORRELATION_USE_VULKAN

namespace correlation::render3d::vulkan {
namespace {

class VulkanBufferManagerTests : public ::testing::Test {
protected:
  void SetUp() override {
    if (!VulkanContext::isSupported()) {
      GTEST_SKIP() << "Vulkan driver or compatible physical device not available";
    }
    ASSERT_TRUE(device_.initialize());
    target_.create(device_.device(), device_.physicalDevice(), 640, 480);
    ASSERT_TRUE(target_.isValid());
    ASSERT_TRUE(pipeline_mgr_.initialize(device_.device(), target_.renderPass()));
  }

  void TearDown() override {
    pipeline_mgr_.destroy();
    target_.destroy();
    device_.destroy();
  }

  [[nodiscard]] VulkanDevice &device() noexcept { return device_; }
  [[nodiscard]] VulkanPipelineManager &pipelineManager() noexcept { return pipeline_mgr_; }

private:
  VulkanDevice device_;
  VulkanOffscreenTarget target_;
  VulkanPipelineManager pipeline_mgr_;
};

TEST_F(VulkanBufferManagerTests, InitializationAndProperties) {
  VulkanBufferManager buffer_mgr;
  EXPECT_FALSE(buffer_mgr.isValid());
  EXPECT_EQ(buffer_mgr.descriptorSet(), VK_NULL_HANDLE);
  EXPECT_EQ(buffer_mgr.sphereQuadBuffer(), VK_NULL_HANDLE);
  EXPECT_EQ(buffer_mgr.cylinderQuadBuffer(), VK_NULL_HANDLE);
  EXPECT_EQ(buffer_mgr.atomInstanceCount(), 0);
  EXPECT_EQ(buffer_mgr.bondInstanceCount(), 0);
  EXPECT_EQ(buffer_mgr.wireframeVertexCount(), 0);

  EXPECT_TRUE(buffer_mgr.initialize(device().device(), device().physicalDevice(),
                                    pipelineManager().descriptorSetLayout()));
  EXPECT_TRUE(buffer_mgr.isValid());
  EXPECT_NE(buffer_mgr.descriptorSet(), VK_NULL_HANDLE);
  EXPECT_NE(buffer_mgr.sphereQuadBuffer(), VK_NULL_HANDLE);
  EXPECT_NE(buffer_mgr.cylinderQuadBuffer(), VK_NULL_HANDLE);

  buffer_mgr.destroy();
  EXPECT_FALSE(buffer_mgr.isValid());
  EXPECT_EQ(buffer_mgr.descriptorSet(), VK_NULL_HANDLE);
}

TEST_F(VulkanBufferManagerTests, UploadSceneAndGeometricGrowth) {
  VulkanBufferManager buffer_mgr;
  ASSERT_TRUE(buffer_mgr.initialize(device().device(), device().physicalDevice(),
                                    pipelineManager().descriptorSetLayout()));

  Scene3D scene;
  scene.addAtom(AtomInstance{
      .position = {0.0F, 0.0F, 0.0F},
      .radius = 1.0F,
      .color = {.r = 255, .g = 0, .b = 0, .a = 255},
      .id = 1,
      .atomic_number = 29,
  });
  scene.addAtom(AtomInstance{
      .position = {2.0F, 0.0F, 0.0F},
      .radius = 1.0F,
      .color = {.r = 0, .g = 255, .b = 0, .a = 255},
      .id = 2,
      .atomic_number = 29,
  });
  scene.addBond(BondInstance{
      .start = {0.0F, 0.0F, 0.0F},
      .end = {2.0F, 0.0F, 0.0F},
      .radius = 0.2F,
      .color = {.r = 200, .g = 200, .b = 200, .a = 255},
  });
  scene.addUnitCellLine(LineSegment{
      .start = {0.0F, 0.0F, 0.0F},
      .end = {5.0F, 0.0F, 0.0F},
      .color = {.r = 255, .g = 255, .b = 255, .a = 255},
  });

  const Camera3D camera;
  const ViewportDimensions dims{.width = 800, .height = 600};

  buffer_mgr.uploadScene(scene, camera, dims);

  EXPECT_EQ(buffer_mgr.atomInstanceCount(), 2);
  EXPECT_EQ(buffer_mgr.bondInstanceCount(), 1);
  EXPECT_EQ(buffer_mgr.wireframeVertexCount(), 2);
  EXPECT_NE(buffer_mgr.atomInstanceBuffer(), VK_NULL_HANDLE);
  EXPECT_NE(buffer_mgr.bondInstanceBuffer(), VK_NULL_HANDLE);
  EXPECT_NE(buffer_mgr.wireframeBuffer(), VK_NULL_HANDLE);

  // Growth test: add 500 atoms
  for (int i = 0; i < 500; ++i) {
    scene.addAtom(AtomInstance{
        .position = {static_cast<float>(i), 1.0F, 0.0F},
        .radius = 0.5F,
        .color = {.r = 100, .g = 150, .b = 200, .a = 255},
        .id = static_cast<uint32_t>(i + 3),
        .atomic_number = 6,
    });
  }

  buffer_mgr.uploadScene(scene, camera, dims);
  EXPECT_EQ(buffer_mgr.atomInstanceCount(), 502);
  EXPECT_NE(buffer_mgr.atomInstanceBuffer(), VK_NULL_HANDLE);
}

TEST_F(VulkanBufferManagerTests, MoveSemantics) {
  VulkanBufferManager mgr1;
  ASSERT_TRUE(mgr1.initialize(device().device(), device().physicalDevice(),
                              pipelineManager().descriptorSetLayout()));
  VkDescriptorSet orig_desc = mgr1.descriptorSet();

  // Move constructor
  VulkanBufferManager mgr2(std::move(mgr1));
  EXPECT_TRUE(mgr2.isValid());
  EXPECT_EQ(mgr2.descriptorSet(), orig_desc);

  // Move assignment
  VulkanBufferManager mgr3;
  mgr3 = std::move(mgr2);
  EXPECT_TRUE(mgr3.isValid());
  EXPECT_EQ(mgr3.descriptorSet(), orig_desc);
}

} // namespace
} // namespace correlation::render3d::vulkan

#endif // CORRELATION_USE_VULKAN
