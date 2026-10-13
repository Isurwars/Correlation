/**
 * @file SoftwareRasterizerTests.cpp
 * @brief Unit tests for CPU SoftwareFramebuffer and SoftwareRasterizer depth testing.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/Cell.hpp"
#include "render3d/cpu/SoftwareFramebuffer.hpp"
#include "render3d/cpu/SoftwareRasterizer.hpp"
#include <gtest/gtest.h>

namespace correlation::render3d {

TEST(SoftwareRasterizerTests, FramebufferLifecycle) {
  SoftwareFramebuffer fb(200, 100);
  EXPECT_EQ(fb.width(), 200);
  EXPECT_EQ(fb.height(), 100);
  EXPECT_EQ(fb.colorBuffer().size(), 200 * 100 * 4);
  EXPECT_EQ(fb.depthBuffer().size(), 200 * 100);

  const RgbaColor bg{.r = 30, .g = 40, .b = 50, .a = 255};
  fb.clear(bg, 1.0F);

  EXPECT_EQ(fb.colorBuffer()[0], 30);
  EXPECT_EQ(fb.colorBuffer()[1], 40);
  EXPECT_EQ(fb.colorBuffer()[2], 50);
  EXPECT_EQ(fb.colorBuffer()[3], 255);
  EXPECT_FLOAT_EQ(fb.getDepth(0, 0), 1.0F);

  // Depth test write
  const RgbaColor pixel_color{.r = 255, .g = 0, .b = 0, .a = 255};
  const bool written = fb.testAndSetPixel(50, 50, 0.5F, pixel_color);
  EXPECT_TRUE(written);
  EXPECT_FLOAT_EQ(fb.getDepth(50, 50), 0.5F);

  // Further pixel rejected
  const RgbaColor behind_color{.r = 0, .g = 255, .b = 0, .a = 255};
  const bool rejected = fb.testAndSetPixel(50, 50, 0.8F, behind_color);
  EXPECT_FALSE(rejected);
  EXPECT_FLOAT_EQ(fb.getDepth(50, 50), 0.5F);
  EXPECT_EQ(fb.colorBuffer()[(50 * 200 + 50) * 4], 255); // Red preserved
}

TEST(SoftwareRasterizerTests, DepthOcclusionFrontOverBack) {
  SoftwareFramebuffer fb(160, 120);
  Camera3D camera;
  camera.setViewport(160.0F, 120.0F);
  camera.setTarget(math::Vector3<float>(0.0F, 0.0F, 0.0F));

  SoftwareRasterizer rasterizer;

  // Scene with two overlapping spheres: Carbon (dark) in front, Oxygen (red) behind
  Scene3D scene;
  // Camera looks from (0, 0, 20) toward (0, 0, 0)
  // Atom 1: near camera at (0, 0, 5) -> front
  scene.addAtom(AtomInstance{
      .position = math::Vector3<float>(0.0F, 0.0F, 5.0F),
      .radius = 2.0F,
      .color = getElementColor("C"), // Dark gray
      .id = 1,
      .atomic_number = 6,
  });
  // Atom 2: far from camera at (0, 0, -5) -> behind
  scene.addAtom(AtomInstance{
      .position = math::Vector3<float>(0.0F, 0.0F, -5.0F),
      .radius = 2.0F,
      .color = getElementColor("O"), // Bright Red
      .id = 2,
      .atomic_number = 8,
  });

  rasterizer.render(scene, camera, fb);

  // Center pixel (80, 60) should see front atom (C - dark gray, red component < 150)
  const std::size_t center_idx = (60 * 160 + 80) * 4;
  const uint8_t red_val = fb.colorBuffer()[center_idx];
  EXPECT_LT(red_val, 150); // Dark gray, not red
}

TEST(SoftwareRasterizerTests, RenderCellWithWireframeAndBonds) {
  core::Cell cell;
  const math::Matrix3<real_t> lattice(math::Vector3<real_t>(8.0, 0.0, 0.0),
                                      math::Vector3<real_t>(0.0, 8.0, 0.0),
                                      math::Vector3<real_t>(0.0, 0.0, 8.0));
  cell.updateLattice(lattice);
  cell.addAtom("C", math::Vector3<real_t>(4.0, 4.0, 4.0));
  cell.addAtom("H", math::Vector3<real_t>(4.0, 4.0, 5.09));

  const auto scene = Scene3D::fromCell(cell, 1.2, true);
  ASSERT_NE(scene, nullptr);

  Camera3D camera;
  camera.setViewport(200.0F, 150.0F);
  camera.fitToBox(scene->boundsMin(), scene->boundsMax());

  SoftwareFramebuffer fb(200, 150);
  SoftwareRasterizer rasterizer;
  rasterizer.render(*scene, camera, fb);

  // Count non-background pixels
  std::size_t painted_pixels = 0;
  const auto &pixels = fb.colorBuffer();
  for (std::size_t i = 0; i < pixels.size(); i += 4) {
    if (pixels[i] != 24 || pixels[i + 1] != 24 || pixels[i + 2] != 29) {
      ++painted_pixels;
    }
  }

  EXPECT_GT(painted_pixels, 50); // Scene was visibly rasterized
}

} // namespace correlation::render3d
