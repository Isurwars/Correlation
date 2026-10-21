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

namespace {

[[nodiscard]] constexpr std::size_t pixelOffset(std::size_t x_coord, std::size_t y_coord,
                                                std::size_t stride_width) noexcept {
  return (y_coord * stride_width + x_coord) * 4UL;
}

} // namespace

TEST(SoftwareRasterizerTests, FramebufferLifecycle) {
  SoftwareFramebuffer framebuffer(200, 100);
  EXPECT_EQ(framebuffer.width(), 200);
  EXPECT_EQ(framebuffer.height(), 100);
  EXPECT_EQ(framebuffer.colorBuffer().size(), 200UL * 100UL * 4UL);
  EXPECT_EQ(framebuffer.depthBuffer().size(), 200UL * 100UL);

  const RgbaColor bg_color{.r = 30, .g = 40, .b = 50, .a = 255};
  framebuffer.clear(bg_color, 1.0F);

  EXPECT_EQ(framebuffer.colorBuffer()[0], 30);
  EXPECT_EQ(framebuffer.colorBuffer()[1], 40);
  EXPECT_EQ(framebuffer.colorBuffer()[2], 50);
  EXPECT_EQ(framebuffer.colorBuffer()[3], 255);
  EXPECT_FLOAT_EQ(framebuffer.getDepth(0, 0), 1.0F);

  // Depth test write
  const RgbaColor pixel_color{.r = 255, .g = 0, .b = 0, .a = 255};
  const bool written = framebuffer.testAndSetPixel(50, 50, 0.5F, pixel_color);
  EXPECT_TRUE(written);
  EXPECT_FLOAT_EQ(framebuffer.getDepth(50, 50), 0.5F);

  // Further pixel rejected
  const RgbaColor behind_color{.r = 0, .g = 255, .b = 0, .a = 255};
  const bool rejected = framebuffer.testAndSetPixel(50, 50, 0.8F, behind_color);
  EXPECT_FALSE(rejected);
  EXPECT_FLOAT_EQ(framebuffer.getDepth(50, 50), 0.5F);
  EXPECT_EQ(framebuffer.colorBuffer()[pixelOffset(50, 50, 200)], 255); // Red preserved
}

TEST(SoftwareRasterizerTests, DepthOcclusionFrontOverBack) {
  SoftwareFramebuffer framebuffer(160, 120);
  Camera3D camera;
  camera.setViewport(160.0F, 120.0F);
  camera.setTarget(math::Vector3<float>(0.0F, 0.0F, 0.0F));

  const SoftwareRasterizer rasterizer;

  // Scene with two overlapping spheres: Carbon (dark) in front, Oxygen (red) behind
  Scene3D scene;
  // Camera looks from (0, 0, 20) toward (0, 0, 0)
  // Atom 1: near camera at (0, 0, 5) -> front
  scene.addAtom(AtomInstance{
      .position = math::Vector3<float>(0.0F, 0.0F, 5.0F),
      .radius = 2.0F,
      .covalent_radius = 2.0F,
      .color = getElementColor("C"), // Dark gray
      .id = 1,
      .atomic_number = 6,
  });
  // Atom 2: far from camera at (0, 0, -5) -> behind
  scene.addAtom(AtomInstance{
      .position = math::Vector3<float>(0.0F, 0.0F, -5.0F),
      .radius = 2.0F,
      .covalent_radius = 2.0F,
      .color = getElementColor("O"), // Bright Red
      .id = 2,
      .atomic_number = 8,
  });

  rasterizer.render(scene, camera, framebuffer);

  // Center pixel (80, 60) should see front atom (C - dark gray, red component < 150)
  const std::size_t center_idx = pixelOffset(80, 60, 160);
  const uint8_t red_val = framebuffer.colorBuffer()[center_idx];
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

  SoftwareFramebuffer framebuffer(200, 150);
  const SoftwareRasterizer rasterizer;
  rasterizer.render(*scene, camera, framebuffer);

  // Count non-background pixels
  std::size_t painted_pixels = 0;
  const auto &pixels = framebuffer.colorBuffer();
  for (std::size_t offset = 0; offset < pixels.size(); offset += 4) {
    if (pixels[offset] != 24 || pixels[offset + 1] != 24 || pixels[offset + 2] != 29) {
      ++painted_pixels;
    }
  }

  EXPECT_GT(painted_pixels, 50); // Scene was visibly rasterized
}

TEST(SoftwareRasterizerTests, BlendAdditivePixelDepthAndSaturation) {
  SoftwareFramebuffer framebuffer(100, 100);
  const RgbaColor bg_color{.r = 20, .g = 30, .b = 40, .a = 255};
  framebuffer.clear(bg_color, 1.0F);

  // Out of bounds check
  EXPECT_FALSE(framebuffer.blendAdditivePixel(100, 100, 0.5F,
                                              RgbaColor{.r = 255, .g = 255, .b = 255, .a = 255}));

  // Zero alpha does not blend
  EXPECT_FALSE(framebuffer.blendAdditivePixel(50, 50, 0.5F,
                                              RgbaColor{.r = 255, .g = 255, .b = 255, .a = 0}));

  // Successful blend at depth 0.5F
  const RgbaColor add_color{.r = 100, .g = 100, .b = 100, .a = 255};
  const bool blended = framebuffer.blendAdditivePixel(50, 50, 0.5F, add_color);
  EXPECT_TRUE(blended);

  // Depth buffer MUST NOT be modified by additive blending
  EXPECT_FLOAT_EQ(framebuffer.getDepth(50, 50), 1.0F);

  // Color values should have saturated/accumulated
  const std::size_t idx = pixelOffset(50, 50, 100);
  EXPECT_EQ(framebuffer.colorBuffer()[idx], 120);     // 20 + 100
  EXPECT_EQ(framebuffer.colorBuffer()[idx + 1], 130); // 30 + 100
  EXPECT_EQ(framebuffer.colorBuffer()[idx + 2], 140); // 40 + 100

  // Pixel behind existing occluder (if depth was already 0.3F)
  framebuffer.testAndSetPixel(20, 20, 0.3F, RgbaColor{.r = 255, .g = 0, .b = 0, .a = 255});
  EXPECT_FLOAT_EQ(framebuffer.getDepth(20, 20), 0.3F);

  // Blend at depth 0.6F (behind 0.3F) should be rejected
  const bool occluded = framebuffer.blendAdditivePixel(20, 20, 0.6F, add_color);
  EXPECT_FALSE(occluded);
}

TEST(SoftwareRasterizerTests, SphereHaloRenderingAddsAdditiveGlow) {
  SoftwareFramebuffer framebuffer(160, 120);
  Camera3D camera;
  camera.setViewport(160.0F, 120.0F);
  camera.setTarget(math::Vector3<float>(0.0F, 0.0F, 0.0F));

  const SoftwareRasterizer rasterizer;

  Scene3D scene;
  // Atom at origin with radius 2.0 and green additive halo
  scene.addAtom(AtomInstance{
      .position = math::Vector3<float>(0.0F, 0.0F, 0.0F),
      .radius = 2.0F,
      .covalent_radius = 2.0F,
      .color = RgbaColor{.r = 255, .g = 0, .b = 0, .a = 255}, // Red atom
      .id = 1,
      .atomic_number = 8,
      .halo_color = RgbaColor{.r = 0, .g = 255, .b = 0, .a = 200}, // Green halo
      .halo_radius_ratio = 1.5F,
      .dim_factor = 1.0F,
  });

  rasterizer.render(scene, camera, framebuffer);

  // 1. Center pixel (80, 60): atom core (red dominant)
  const std::size_t center_idx = pixelOffset(80, 60, 160);
  EXPECT_GT(framebuffer.colorBuffer()[center_idx], 100);

  // 2. Halo ring pixel (at x = 98, y = 60, offset 18 px > screen_r ≈ 14.5 and < halo_r ≈ 21.7)
  const std::size_t halo_idx = pixelOffset(98, 60, 160);
  // Green channel should be noticeably elevated above background (background g = 24)
  EXPECT_GT(framebuffer.colorBuffer()[halo_idx + 1], 24);

  // 3. Pixel far outside halo (at x = 120, y = 60, offset 40 px)
  const std::size_t outer_idx = pixelOffset(120, 60, 160);
  EXPECT_EQ(framebuffer.colorBuffer()[outer_idx], 24);     // Dark slate r
  EXPECT_EQ(framebuffer.colorBuffer()[outer_idx + 1], 24); // Dark slate g
  EXPECT_EQ(framebuffer.colorBuffer()[outer_idx + 2], 29); // Dark slate b
}

TEST(SoftwareRasterizerTests, SphereDimFactorModulatesCoreLuminance) {
  Camera3D camera;
  camera.setViewport(160.0F, 120.0F);
  camera.setTarget(math::Vector3<float>(0.0F, 0.0F, 0.0F));

  const SoftwareRasterizer rasterizer;

  // Scene Normal (dim_factor = 1.0)
  Scene3D scene_normal;
  scene_normal.addAtom(AtomInstance{
      .position = math::Vector3<float>(0.0F, 0.0F, 0.0F),
      .radius = 2.0F,
      .covalent_radius = 2.0F,
      .color = RgbaColor{.r = 200, .g = 200, .b = 200, .a = 255},
      .id = 1,
      .atomic_number = 6,
      .halo_color = RgbaColor{.r = 0, .g = 0, .b = 0, .a = 0},
      .halo_radius_ratio = 1.0F,
      .dim_factor = 1.0F,
  });

  SoftwareFramebuffer framebuffer_normal(160, 120);
  rasterizer.render(scene_normal, camera, framebuffer_normal);

  // Scene Dimmed (dim_factor = 0.3)
  Scene3D scene_dimmed;
  scene_dimmed.addAtom(AtomInstance{
      .position = math::Vector3<float>(0.0F, 0.0F, 0.0F),
      .radius = 2.0F,
      .covalent_radius = 2.0F,
      .color = RgbaColor{.r = 200, .g = 200, .b = 200, .a = 255},
      .id = 1,
      .atomic_number = 6,
      .halo_color = RgbaColor{.r = 0, .g = 0, .b = 0, .a = 0},
      .halo_radius_ratio = 1.0F,
      .dim_factor = 0.3F,
  });

  SoftwareFramebuffer framebuffer_dimmed(160, 120);
  rasterizer.render(scene_dimmed, camera, framebuffer_dimmed);

  const std::size_t center_idx = pixelOffset(80, 60, 160);
  const uint8_t lum_normal = framebuffer_normal.colorBuffer()[center_idx];
  const uint8_t lum_dimmed = framebuffer_dimmed.colorBuffer()[center_idx];

  EXPECT_GT(lum_normal, lum_dimmed);
  EXPECT_LT(lum_dimmed, lum_normal / 2);
}

} // namespace correlation::render3d
