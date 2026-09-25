// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "app/core/AppOptions.hpp"
#include "app/services/PhysicsService.hpp"
#include "core/Cell.hpp"
#include <gtest/gtest.h>

namespace {

using correlation::app::AppDefaults;
using correlation::app::PhysicsService;
using correlation::core::Cell;

TEST(PhysicsServiceTests, ComputeRecommendedTimeStepHandlesNullAndEmptyCell) {
  EXPECT_DOUBLE_EQ(PhysicsService::computeRecommendedTimeStep(nullptr), AppDefaults::TIME_STEP);

  Cell empty_cell;
  EXPECT_DOUBLE_EQ(PhysicsService::computeRecommendedTimeStep(&empty_cell), AppDefaults::TIME_STEP);
}

TEST(PhysicsServiceTests, ComputeRecommendedTimeStepWithKnownElements) {
  Cell cell;
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("O", {1.0, 1.0, 1.0});

  const real_t recommended_dt = PhysicsService::computeRecommendedTimeStep(&cell);
  EXPECT_GT(recommended_dt, 0.0);
  // Oxygen mass is ~15.9994 amu, lower than Silicon ~28.0855 amu.
  // sqrt(9 * 15.9994 / 5) ~= 5.3665 fs
  EXPECT_NEAR(recommended_dt, 5.3665, 0.01);
}

TEST(PhysicsServiceTests, ScaleBondCutoffsValidAndEdgeCases) {
  const correlation::analysis::BondCutoffMatrix empty_mat;
  EXPECT_TRUE(PhysicsService::scaleBondCutoffs(empty_mat, 1.5).empty());

  const correlation::analysis::BondCutoffMatrix mat = {
      {{.min_sq = 1.0, .max_sq = 4.0}, {.min_sq = 2.0, .max_sq = 5.0}},
      {{.min_sq = 2.0, .max_sq = 5.0}, {.min_sq = 1.5, .max_sq = 4.5}}};

  // Zero or negative scale factor returns unchanged matrix
  auto non_scaled = PhysicsService::scaleBondCutoffs(mat, 0.0);
  EXPECT_DOUBLE_EQ(non_scaled[0][0].min_sq, 1.0);
  EXPECT_DOUBLE_EQ(non_scaled[0][0].max_sq, 4.0);

  auto scaled = PhysicsService::scaleBondCutoffs(mat, 2.0);
  EXPECT_DOUBLE_EQ(scaled[0][0].min_sq, 4.0);
  EXPECT_DOUBLE_EQ(scaled[0][0].max_sq, 16.0);
  EXPECT_DOUBLE_EQ(scaled[0][1].min_sq, 8.0);
  EXPECT_DOUBLE_EQ(scaled[0][1].max_sq, 20.0);
}

TEST(PhysicsServiceTests, BuildUniformBondCutoffs) {
  EXPECT_TRUE(PhysicsService::buildUniformBondCutoffs(0, 1.0, 3.0).empty());

  const auto uniform_mat = PhysicsService::buildUniformBondCutoffs(2, 1.5, 3.0);
  ASSERT_EQ(uniform_mat.size(), 2U);
  ASSERT_EQ(uniform_mat[0].size(), 2U);
  ASSERT_EQ(uniform_mat[1].size(), 2U);

  EXPECT_NEAR(uniform_mat[0][0].min_sq, 2.25, 1e-5);
  EXPECT_NEAR(uniform_mat[0][0].max_sq, 9.0, 1e-5);
  EXPECT_NEAR(uniform_mat[0][1].min_sq, 2.25, 1e-5);
  EXPECT_NEAR(uniform_mat[0][1].max_sq, 9.0, 1e-5);
  EXPECT_NEAR(uniform_mat[1][0].min_sq, 2.25, 1e-5);
  EXPECT_NEAR(uniform_mat[1][0].max_sq, 9.0, 1e-5);
  EXPECT_NEAR(uniform_mat[1][1].min_sq, 2.25, 1e-5);
  EXPECT_NEAR(uniform_mat[1][1].max_sq, 9.0, 1e-5);
}

} // namespace
