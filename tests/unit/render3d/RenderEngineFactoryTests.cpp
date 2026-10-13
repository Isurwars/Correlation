/**
 * @file RenderEngineFactoryTests.cpp
 * @brief Unit tests for IRenderEngine polymorphism and RenderEngineFactory fallback logic.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/Cell.hpp"
#include "render3d/RenderEngineFactory.hpp"
#include <gtest/gtest.h>

namespace correlation::render3d {

TEST(RenderEngineFactoryTests, CreateCpuEngineDirectly) {
  const auto engine = RenderEngineFactory::createEngine(BackendPreference::CpuOnly, {400, 300});
  ASSERT_NE(engine, nullptr);

  EXPECT_EQ(engine->dimensions().width, 400);
  EXPECT_EQ(engine->dimensions().height, 300);
  EXPECT_FALSE(engine->isHardwareAccelerated());
  EXPECT_EQ(engine->backendName(), "CPU (Software Rasterizer)");

  const auto frame = engine->renderFrame();
  EXPECT_EQ(frame.size(), 400 * 300 * 4);
}

TEST(RenderEngineFactoryTests, CreateAutoEngineWithFallback) {
  const auto engine = RenderEngineFactory::createEngine(BackendPreference::Auto, {320, 240});
  ASSERT_NE(engine, nullptr);

  EXPECT_EQ(engine->dimensions().width, 320);
  EXPECT_EQ(engine->dimensions().height, 240);
  EXPECT_FALSE(engine->backendName().empty());

  const auto frame = engine->renderFrame();
  EXPECT_EQ(frame.size(), 320 * 240 * 4);
}

TEST(RenderEngineFactoryTests, VulkanOnlyReturnsNullptrWhenUnavailable) {
  // If Vulkan is not supported in the current environment, VulkanOnly must safely return nullptr
  const auto engine = RenderEngineFactory::createEngine(BackendPreference::VulkanOnly, {200, 200});
  if (engine == nullptr) {
    SUCCEED();
  } else {
    EXPECT_TRUE(engine->isHardwareAccelerated());
    EXPECT_EQ(engine->backendName(), "Vulkan (Headless Offscreen)");
  }
}

TEST(RenderEngineFactoryTests, RenderSceneThroughInterface) {
  const auto engine = RenderEngineFactory::createEngine(BackendPreference::Auto, {200, 150});
  ASSERT_NE(engine, nullptr);

  core::Cell cell;
  const math::Matrix3<real_t> lattice(math::Vector3<real_t>(6.0, 0.0, 0.0),
                                      math::Vector3<real_t>(0.0, 6.0, 0.0),
                                      math::Vector3<real_t>(0.0, 0.0, 6.0));
  cell.updateLattice(lattice);
  cell.addAtom("Si", math::Vector3<real_t>(3.0, 3.0, 3.0));

  const auto scene = Scene3D::fromCell(cell, 1.2, true);
  engine->setScene(scene);

  Camera3D camera;
  camera.setViewport(200.0F, 150.0F);
  camera.fitToBox(scene->boundsMin(), scene->boundsMax());
  engine->setCamera(camera);

  const auto pixels = engine->renderFrame();
  EXPECT_EQ(pixels.size(), 200 * 150 * 4);

  // Verify non-zero pixels exist in buffer
  std::size_t painted_pixels = 0;
  for (std::size_t i = 0; i < pixels.size(); i += 4) {
    if (pixels[i] != 24 || pixels[i + 1] != 24 || pixels[i + 2] != 29) {
      ++painted_pixels;
    }
  }
  EXPECT_GT(painted_pixels, 10);
}

} // namespace correlation::render3d
