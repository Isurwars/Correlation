/**
 * @file Camera3DTests.cpp
 * @brief Unit tests for Camera3D orbit model, panning, zooming, and bounding-box fit.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/Camera3D.hpp"
#include <gtest/gtest.h>

namespace correlation::render3d {

using correlation::real_t;

TEST(Camera3DTests, DefaultInitialization) {
  Camera3D const camera;
  EXPECT_FLOAT_EQ(camera.target().x(), 0.0F);
  EXPECT_FLOAT_EQ(camera.target().y(), 0.0F);
  EXPECT_FLOAT_EQ(camera.target().z(), 0.0F);
  EXPECT_GT(camera.distance(), 0.0F);
  EXPECT_GT(camera.fov(), 0.0F);
  EXPECT_GT(camera.aspectRatio(), 0.0F);
}

TEST(Camera3DTests, OrbitRotation) {
  Camera3D camera;
  camera.orbit(90.0F, 0.0F);
  EXPECT_NEAR(camera.azimuth(), 90.0F, 1e-4F);

  // Elevation clamp test
  camera.orbit(0.0F, 100.0F);
  EXPECT_LE(camera.elevation(), 89.0F);

  camera.orbit(0.0F, -200.0F);
  EXPECT_GE(camera.elevation(), -89.0F);
}

TEST(Camera3DTests, PanAndZoom) {
  Camera3D camera;
  const auto initial_target = camera.target();
  camera.pan(10.0F, 5.0F);
  EXPECT_NE(camera.target().x(), initial_target.x());

  const float initial_dist = camera.distance();
  camera.zoom(0.2F);
  EXPECT_LT(camera.distance(), initial_dist);

  camera.zoom(-0.5F);
  EXPECT_GT(camera.distance(), initial_dist * 0.8F);
}

TEST(Camera3DTests, FitToBox) {
  Camera3D camera;
  const math::Vector3<float> min_pt(-10.0F, -10.0F, -10.0F);
  const math::Vector3<float> max_pt(10.0F, 10.0F, 10.0F);

  camera.fitToBox(min_pt, max_pt);
  EXPECT_NEAR(camera.target().x(), 0.0F, 1e-4F);
  EXPECT_NEAR(camera.target().y(), 0.0F, 1e-4F);
  EXPECT_NEAR(camera.target().z(), 0.0F, 1e-4F);
  EXPECT_GT(camera.distance(), 20.0F);
}

TEST(Camera3DTests, MatrixGeneration) {
  Camera3D camera;
  camera.setViewport(800.0F, 600.0F);

  const auto view = camera.viewMatrix();
  const auto proj = camera.projectionMatrix();
  const auto view_proj = camera.viewProjectionMatrix();

  static_cast<void>(proj);

  const math::Vector3<float> origin(0.0F, 0.0F, 0.0F);
  const auto p_view = view.transformPoint(origin);
  EXPECT_LT(p_view.z(), 0.0F); // Point in front of camera has negative view Z

  const auto p_clip = view_proj.transformPoint(origin);
  EXPECT_GE(p_clip.z(), -1.0F);
  EXPECT_LE(p_clip.z(), 1.0F);
}

} // namespace correlation::render3d
