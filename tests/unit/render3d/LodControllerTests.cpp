/**
 * @file LodControllerTests.cpp
 * @brief Unit tests for LodController level evaluation and dynamic motion policies.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/LodController.hpp"
#include <gtest/gtest.h>

namespace correlation::render3d {

TEST(LodControllerTests, SmallSystemHighLod) {
  LodController controller;
  const auto policy_static = controller.evaluatePolicy(1000, false);
  EXPECT_EQ(policy_static.level, LodLevel::High);
  EXPECT_TRUE(policy_static.render_bonds);
  EXPECT_TRUE(policy_static.render_wireframe);
  EXPECT_FLOAT_EQ(policy_static.quality_scale, 1.0F);

  const auto policy_moving = controller.evaluatePolicy(1000, true);
  EXPECT_EQ(policy_moving.level, LodLevel::High);
  EXPECT_TRUE(policy_moving.render_bonds);
}

TEST(LodControllerTests, MediumSystemLod) {
  LodController controller;
  const auto policy_static = controller.evaluatePolicy(10000, false);
  EXPECT_EQ(policy_static.level, LodLevel::Medium);
  EXPECT_TRUE(policy_static.render_bonds);
  EXPECT_FLOAT_EQ(policy_static.quality_scale, 0.8F);

  const auto policy_moving = controller.evaluatePolicy(10000, true);
  EXPECT_EQ(policy_moving.level, LodLevel::Medium);
  EXPECT_FALSE(policy_moving.render_bonds); // Bonds suppressed during motion
  EXPECT_FLOAT_EQ(policy_moving.quality_scale, 0.6F);
}

TEST(LodControllerTests, LargeSystemLowLod) {
  LodController controller;
  const auto policy = controller.evaluatePolicy(50000, false);
  EXPECT_EQ(policy.level, LodLevel::Low);
  EXPECT_FALSE(policy.render_bonds);
  EXPECT_TRUE(policy.render_wireframe);
  EXPECT_FLOAT_EQ(policy.quality_scale, 0.5F);

  const auto policy_moving = controller.evaluatePolicy(50000, true);
  EXPECT_FLOAT_EQ(policy_moving.quality_scale, 0.35F);
}

TEST(LodControllerTests, CustomThresholds) {
  LodController controller;
  controller.setHighThreshold(500);
  controller.setMediumThreshold(2000);

  const auto policy = controller.evaluatePolicy(800, false);
  EXPECT_EQ(policy.level, LodLevel::Medium);
}

} // namespace correlation::render3d
