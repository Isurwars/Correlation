/**
 * @file RayPickerTests.cpp
 * @brief Unit tests for ray-sphere atom picking and camera HUD metrics.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/Camera3D.hpp"
#include "render3d/RayPicker.hpp"
#include "render3d/Scene3D.hpp"
#include <gtest/gtest.h>

namespace correlation::render3d {

TEST(RayPickerTests, PicksClosestVisibleAtomInPerspectiveAndOrtho) {
  Camera3D camera;
  camera.setViewport(800.0F, 600.0F);
  camera.setTarget({0.0F, 0.0F, 0.0F});
  camera.setOrientation(0.0F, 0.0F);

  Scene3D scene;
  AtomInstance back_atom{};
  back_atom.position = {0.0F, 0.0F, -2.0F};
  back_atom.radius = 0.8F;
  back_atom.id = 0;
  scene.addAtom(back_atom);

  AtomInstance front_atom{};
  front_atom.position = {0.0F, 0.0F, 2.0F};
  front_atom.radius = 0.8F;
  front_atom.id = 1;
  scene.addAtom(front_atom);

  const auto hit_persp = pickAtom(camera, scene, 0.5F, 0.5F);
  ASSERT_TRUE(hit_persp.has_value());
  EXPECT_EQ(*hit_persp, 1U);

  camera.setProjectionMode(ProjectionMode::Orthographic);
  const auto hit_ortho = pickAtom(camera, scene, 0.5F, 0.5F);
  ASSERT_TRUE(hit_ortho.has_value());
  EXPECT_EQ(*hit_ortho, 1U);

  const auto miss = pickAtom(camera, scene, 0.02F, 0.02F);
  EXPECT_FALSE(miss.has_value());
}

TEST(RayPickerTests, IgnoresCulledAtoms) {
  Camera3D camera;
  camera.setViewport(800.0F, 600.0F);
  camera.setTarget({0.0F, 0.0F, 0.0F});
  camera.setOrientation(0.0F, 0.0F);

  Scene3D scene;
  AtomInstance hidden_front{};
  hidden_front.position = {0.0F, 0.0F, 2.0F};
  hidden_front.radius = 1.0F;
  hidden_front.dim_factor = 0.0F;
  hidden_front.halo_color = {.r = 0, .g = 0, .b = 0, .a = 0};
  scene.addAtom(hidden_front);

  AtomInstance visible_back{};
  visible_back.position = {0.0F, 0.0F, -1.0F};
  visible_back.radius = 1.0F;
  visible_back.dim_factor = 1.0F;
  scene.addAtom(visible_back);

  const auto hit = pickAtom(camera, scene, 0.5F, 0.5F);
  ASSERT_TRUE(hit.has_value());
  EXPECT_EQ(*hit, 1U);
}

TEST(RayPickerTests, ComputesValidCameraHudMetrics) {
  Camera3D camera;
  camera.setViewport(800.0F, 600.0F);
  camera.setOrientation(45.0F, 30.0F);

  const CameraHudMetrics metrics = computeCameraHudMetrics(camera, 600.0F);
  EXPECT_GE(metrics.scale_bar_px, 32.0F);
  EXPECT_LE(metrics.scale_bar_px, 140.0F);
  EXPECT_FALSE(metrics.scale_bar_label.empty());
}

} // namespace correlation::render3d
