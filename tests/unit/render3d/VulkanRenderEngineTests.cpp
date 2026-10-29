/**
 * @file VulkanRenderEngineTests.cpp
 * @brief Integration tests for VulkanRenderEngine offscreen rendering and factory selection.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/Camera3D.hpp"
#include "render3d/RenderEngineFactory.hpp"
#include "render3d/Scene3D.hpp"
#include "render3d/vulkan/VulkanContext.hpp"
#include "render3d/vulkan/VulkanRenderEngine.hpp"
#include <gtest/gtest.h>

namespace correlation::render3d {
namespace {

class VulkanRenderEngineTests : public ::testing::Test {
protected:
  void SetUp() override {
    if (!VulkanContext::isSupported()) {
      GTEST_SKIP() << "Vulkan driver or compatible physical device not available";
    }
  }
};

TEST_F(VulkanRenderEngineTests, InitializationAndProperties) {
  VulkanRenderEngine engine;
  EXPECT_EQ(engine.backendName(), "Vulkan (Headless Offscreen)");
  EXPECT_TRUE(engine.isHardwareAccelerated());

  const ViewportDimensions dims{.width = 640, .height = 480};
  ASSERT_TRUE(engine.initialize(dims));
  EXPECT_EQ(engine.dimensions().width, 640);
  EXPECT_EQ(engine.dimensions().height, 480);

  const auto frame = engine.renderFrame();
  EXPECT_EQ(frame.size(), 640 * 480 * 4);
}

TEST_F(VulkanRenderEngineTests, RenderFrameEmptyScene) {
  VulkanRenderEngine engine;
  const ViewportDimensions dims{.width = 320, .height = 240};
  ASSERT_TRUE(engine.initialize(dims));

  const auto frame = engine.renderFrame();
  ASSERT_EQ(frame.size(), 320 * 240 * 4);

  // Background clear color is {0.1F, 0.1F, 0.12F, 1.0F} -> ~{25, 25, 30, 255}
  EXPECT_NEAR(frame[0], 25, 5);
  EXPECT_NEAR(frame[1], 25, 5);
  EXPECT_NEAR(frame[2], 30, 5);
  EXPECT_EQ(frame[3], 255);
}

TEST_F(VulkanRenderEngineTests, RenderFrameWithAtomsAndBonds) {
  VulkanRenderEngine engine;
  const ViewportDimensions dims{.width = 400, .height = 400};
  ASSERT_TRUE(engine.initialize(dims));

  const auto scene = std::make_shared<Scene3D>();
  scene->addAtom(AtomInstance{
      .position = {0.0F, 0.0F, 0.0F},
      .radius = 2.0F,
      .color = {.r = 255, .g = 0, .b = 0, .a = 255}, // Bright red
      .id = 1,
      .atomic_number = 29,
  });
  scene->addAtom(AtomInstance{
      .position = {3.0F, 0.0F, 0.0F},
      .radius = 1.5F,
      .color = {.r = 0, .g = 255, .b = 0, .a = 255}, // Bright green
      .id = 2,
      .atomic_number = 6,
  });
  scene->addBond(BondInstance{
      .start = {0.0F, 0.0F, 0.0F},
      .end = {3.0F, 0.0F, 0.0F},
      .radius = 0.3F,
      .color = {.r = 200, .g = 200, .b = 200, .a = 255},
  });

  Camera3D camera;
  camera.setTarget({1.5F, 0.0F, 0.0F});
  camera.orbit(0.0F, 0.0F); // Face directly
  camera.setViewport(400.0F, 400.0F);

  engine.setScene(scene);
  engine.setCamera(camera);

  const auto frame = engine.renderFrame();
  ASSERT_EQ(frame.size(), 400 * 400 * 4);

  // Verify non-background pixels are present in the frame
  bool has_foreground = false;
  for (std::size_t i = 0; i < frame.size(); i += 4) {
    const uint8_t red = frame[i];
    const uint8_t green = frame[i + 1];
    const uint8_t blue = frame[i + 2];
    // Check if significantly different from background (25, 25, 30)
    if (red > 60 || green > 60 || blue > 60) {
      has_foreground = true;
      break;
    }
  }
  EXPECT_TRUE(has_foreground);
}

TEST_F(VulkanRenderEngineTests, DynamicResize) {
  VulkanRenderEngine engine;
  ASSERT_TRUE(engine.initialize(ViewportDimensions{.width = 640, .height = 480}));

  engine.resize(ViewportDimensions{.width = 200, .height = 150});
  EXPECT_EQ(engine.dimensions().width, 200);
  EXPECT_EQ(engine.dimensions().height, 150);

  const auto frame = engine.renderFrame();
  EXPECT_EQ(frame.size(), 200 * 150 * 4);
}

TEST_F(VulkanRenderEngineTests, FactoryVulkanOnly) {
  auto engine = RenderEngineFactory::createEngine(BackendPreference::VulkanOnly,
                                                  ViewportDimensions{.width = 320, .height = 240});
  ASSERT_NE(engine, nullptr);
  EXPECT_TRUE(engine->isHardwareAccelerated());
  EXPECT_EQ(engine->backendName(), "Vulkan (Headless Offscreen)");

  const auto frame = engine->renderFrame();
  EXPECT_EQ(frame.size(), 320 * 240 * 4);
}

TEST_F(VulkanRenderEngineTests, FactoryAutoSelection) {
  auto engine = RenderEngineFactory::createEngine(BackendPreference::Auto,
                                                  ViewportDimensions{.width = 320, .height = 240});
  ASSERT_NE(engine, nullptr);
  EXPECT_TRUE(engine->isHardwareAccelerated());
  EXPECT_EQ(engine->backendName(), "Vulkan (Headless Offscreen)");
}

} // namespace
} // namespace correlation::render3d
