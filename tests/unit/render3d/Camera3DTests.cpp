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

TEST(Camera3DTests, DefaultInitialization) {
  Camera3D camera;
  EXPECT_FLOAT_EQ(camera.target().x(), 0.0f);
  EXPECT_FLOAT_EQ(camera.target().y(), 0.0f);
  EXPECT_FLOAT_EQ(camera.target().z(), 0.0f);
  EXPECT_GT(camera.distance(), 0.0f);
  EXPECT_GT(camera.fov(), 0.0f);
  EXPECT_GT(camera.aspectRatio(), 0.0f);
}

TEST(Camera3DTests, OrbitRotation) {
  Camera3D camera;
  camera.orbit(90.0f, 0.0f);
  EXPECT_NEAR(camera.azimuth(), 90.0f, 1e-4f);

  // Elevation clamp test
  camera.orbit(0.0f, 100.0f);
  EXPECT_LE(camera.elevation(), 89.0f);

  camera.orbit(0.0f, -200.0f);
  EXPECT_GE(camera.elevation(), -89.0f);
}

TEST(Camera3DTests, PanAndZoom) {
  Camera3D camera;
  const auto initial_target = camera.target();
  camera.pan(10.0f, 5.0f);
  EXPECT_NE(camera.target().x(), initial_target.x());

  const float initial_dist = camera.distance();
  camera.zoom(0.2f);
  EXPECT_LT(camera.distance(), initial_dist);

  camera.zoom(-0.5f);
  EXPECT_GT(camera.distance(), initial_dist * 0.8f);
}

TEST(Camera3DTests, FitToBox) {
  Camera3D camera;
  const math::Vector3<float> min_pt(-10.0f, -10.0f, -10.0f);
  const math::Vector3<float> max_pt(10.0f, 10.0f, 10.0f);

  camera.fitToBox(min_pt, max_pt);
  EXPECT_NEAR(camera.target().x(), 0.0f, 1e-4f);
  EXPECT_NEAR(camera.target().y(), 0.0f, 1e-4f);
  EXPECT_NEAR(camera.target().z(), 0.0f, 1e-4f);
  EXPECT_GT(camera.distance(), 20.0f);
}

TEST(Camera3DTests, MatrixGeneration) {
  Camera3D camera;
  camera.setViewport(800.0f, 600.0f);

  const auto view = camera.viewMatrix();
  const auto proj = camera.projectionMatrix();
  const auto view_proj = camera.viewProjectionMatrix();

  const math::Vector3<float> origin(0.0f, 0.0f, 0.0f);
  const auto p_view = view.transformPoint(origin);
  EXPECT_LT(p_view.z(), 0.0f); // Point in front of camera has negative view Z

  const auto p_clip = view_proj.transformPoint(origin);
  EXPECT_GE(p_clip.z(), -1.0f);
  EXPECT_LE(p_clip.z(), 1.0f);
}

} // namespace correlation::render3d
