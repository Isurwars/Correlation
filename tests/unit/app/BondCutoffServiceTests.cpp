// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only

#include "app/services/BondCutoffService.hpp"
#include "app/services/TrajectoryLoader.hpp"

#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace {

std::string getTestDataDir() {
  std::vector<std::string> const candidates = {
      "../../tests/data/",
      "../tests/data/",
      "tests/data/",
      "data/",
  };
  for (const auto &dir : candidates) {
    if (std::filesystem::exists(dir + "xyz/clean.xyz")) {
      return dir;
    }
  }
  return "../../tests/data/";
}

TEST(BondCutoffServiceTests, NullTrajectoryReturnsEmptyCutoffs) {
  EXPECT_TRUE(correlation::app::BondCutoffService::getRecommendedBondCutoffs(nullptr).empty());
  EXPECT_DOUBLE_EQ(correlation::app::BondCutoffService::getBondCutoff(nullptr, 0, 0), 0.0);
  EXPECT_DOUBLE_EQ(correlation::app::BondCutoffService::getMinBondCutoff(nullptr, 0, 0), 0.0);
}

TEST(BondCutoffServiceTests, ComputesRecommendedCutoffsForLoadedTrajectory) {
  correlation::app::TrajectoryLoader loader;
  std::string const filepath = getTestDataDir() + "xyz/clean.xyz";

  auto res = loader.loadFile(filepath);
  ASSERT_TRUE(res.has_value());

  auto cutoffs =
      correlation::app::BondCutoffService::getRecommendedBondCutoffs(loader.trajectoryMut());
  EXPECT_FALSE(cutoffs.empty());
  EXPECT_GT(correlation::app::BondCutoffService::getBondCutoff(loader.trajectory(), 0, 0), 0.0);
}

TEST(BondCutoffServiceTests, AppliesScaledCutoffs) {
  correlation::app::TrajectoryLoader loader;
  std::string const filepath = getTestDataDir() + "xyz/clean.xyz";

  auto res = loader.loadFile(filepath);
  ASSERT_TRUE(res.has_value());

  auto base_cutoffs =
      correlation::app::BondCutoffService::getRecommendedBondCutoffs(loader.trajectoryMut());
  auto scaled_cutoffs =
      correlation::app::BondCutoffService::applyScaledBondCutoffs(loader.trajectoryMut(), 1.5);

  EXPECT_FALSE(scaled_cutoffs.empty());
  EXPECT_GT(scaled_cutoffs[0][0].max_sq, base_cutoffs[0][0].max_sq);
}

TEST(BondCutoffServiceTests, SetsUniformCutoffs) {
  correlation::app::TrajectoryLoader loader;
  std::string const filepath = getTestDataDir() + "xyz/clean.xyz";

  auto res = loader.loadFile(filepath);
  ASSERT_TRUE(res.has_value());

  auto uniform_cutoffs = correlation::app::BondCutoffService::setUniformBondCutoff(
      loader.cell(), loader.trajectoryMut(), 1.0, 3.5);

  ASSERT_FALSE(uniform_cutoffs.empty());
  EXPECT_DOUBLE_EQ(uniform_cutoffs[0][0].min_sq, 1.0 * 1.0);
  EXPECT_DOUBLE_EQ(uniform_cutoffs[0][0].max_sq, 3.5 * 3.5);
}

} // namespace
