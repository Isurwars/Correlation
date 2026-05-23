// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/MSDCalculator.hpp"
#include "core/Trajectory.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {

using namespace correlation::calculators;
using namespace correlation::core;

TEST(MSDCalculatorTests, ComputesCorrectMSDAndDeff) {
  // Construct a trajectory of 5 frames with 1 atom moving linearly
  // x(t) = t, y(t) = 0, z(t) = 0
  std::vector<Cell> frames;
  for (int t = 0; t < 5; ++t) {
    Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
    cell.addAtom("Si", {static_cast<double>(t), 0.0, 0.0});
    frames.push_back(std::move(cell));
  }

  double time_step = 1.0; // 1 fs
  Trajectory traj(std::move(frames), time_step);

  // Act
  auto results = MSDCalculator::calculate(traj, 2, 0, 5);

  // Assert
  ASSERT_TRUE(results.find("MSD") != results.end());
  ASSERT_TRUE(results.find("D_eff") != results.end());

  const auto &msd_hist = results.at("MSD");
  const auto &deff_hist = results.at("D_eff");

  // Check lag 0
  EXPECT_DOUBLE_EQ(msd_hist.partials.at("Total")[0], 0.0);
  EXPECT_DOUBLE_EQ(deff_hist.partials.at("Total")[0], 0.0);

  // Check lag 1
  // MSD(1) = 1^2 = 1.0
  // D_eff(1) = 1.0 / (6 * 1 * 1.0) = 1.0 / 6.0
  EXPECT_NEAR(msd_hist.partials.at("Total")[1], 1.0, 1e-6);
  EXPECT_NEAR(deff_hist.partials.at("Total")[1], 1.0 / 6.0, 1e-6);

  // Check lag 2
  // MSD(2) = 2^2 = 4.0
  // D_eff(2) = 4.0 / (6 * 2 * 1.0) = 4.0 / 12.0 = 1.0 / 3.0
  EXPECT_NEAR(msd_hist.partials.at("Total")[2], 4.0, 1e-6);
  EXPECT_NEAR(deff_hist.partials.at("Total")[2], 1.0 / 3.0, 1e-6);
}

} // namespace correlation::testing
