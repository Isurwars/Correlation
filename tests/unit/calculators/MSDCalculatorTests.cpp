// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DynamicsAnalyzer.hpp"
#include "calculators/MSDCalculator.hpp"
#include "core/Trajectory.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {

using namespace correlation::calculators;
using namespace correlation::core;
using namespace correlation::analysis;

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
  auto results = MSDCalculator::calculate(traj, MaxFrames{2}, StartFrame{0}, EndFrame{5});

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

// --- Extreme / Edge-Case Tests ---

TEST(MSDCalculatorTests, StationaryAtomsMSDIsZero) {
  // All atoms remain stationary — MSD should be 0 at all lags
  // Use addFrame() to bypass the constructor's automatic deduplication.
  Trajectory traj;
  traj.setTimeStep(1.0);
  for (int t = 0; t < 5; ++t) {
    Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
    cell.addAtom("Ar", {5.0, 5.0, 5.0}); // Same position every frame
    traj.addFrame(cell);
  }

  auto results = MSDCalculator::calculate(traj, MaxFrames{3}, StartFrame{0}, EndFrame{5});

  ASSERT_TRUE(results.find("MSD") != results.end());
  const auto &msd = results.at("MSD").partials.at("Total");

  for (size_t i = 0; i < msd.size(); ++i) {
    EXPECT_NEAR(msd[i], 0.0, 1e-10) << "MSD should be 0 at lag " << i;
  }
}

TEST(MSDCalculatorTests, MultiAtomAveraging) {
  // Two atoms: one moves +x, one moves +y. MSD should average contributions.
  std::vector<Cell> frames;
  for (int t = 0; t < 4; ++t) {
    Cell cell({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0});
    cell.addAtom("Si", {static_cast<double>(t), 5.0, 5.0}); // Moves +x
    cell.addAtom("Si", {5.0, static_cast<double>(t), 5.0}); // Moves +y
    frames.push_back(std::move(cell));
  }

  Trajectory traj(std::move(frames), 1.0);
  auto results = MSDCalculator::calculate(traj, MaxFrames{2}, StartFrame{0}, EndFrame{4});

  ASSERT_TRUE(results.find("MSD") != results.end());
  const auto &msd = results.at("MSD").partials.at("Total");

  // At lag 0: MSD = 0
  EXPECT_NEAR(msd[0], 0.0, 1e-6);
  // At lag 1: each atom displaced by 1.0 -> squared = 1.0, average = 1.0
  EXPECT_NEAR(msd[1], 1.0, 1e-6);
  // At lag 2: each displaced by 2.0 -> squared = 4.0, average = 4.0
  EXPECT_NEAR(msd[2], 4.0, 1e-6);
}

TEST(MSDCalculatorTests, FrameRangeSubset) {
  // Create 10 frames, but only compute MSD over frames 3-7
  std::vector<Cell> frames;
  for (int t = 0; t < 10; ++t) {
    Cell cell({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0});
    cell.addAtom("Si", {static_cast<double>(t), 0.0, 0.0});
    frames.push_back(std::move(cell));
  }

  Trajectory traj(std::move(frames), 1.0);
  // Restrict to frames [3, 7) with max correlation lag 2
  auto results = MSDCalculator::calculate(traj, MaxFrames{2}, StartFrame{3}, EndFrame{7});

  ASSERT_TRUE(results.find("MSD") != results.end());
  const auto &msd = results.at("MSD").partials.at("Total");

  // Lag 0 should always be 0
  EXPECT_NEAR(msd[0], 0.0, 1e-6);
  // Lag 1: displacement is always 1.0 (uniform motion) -> MSD = 1.0
  EXPECT_NEAR(msd[1], 1.0, 1e-6);
}

TEST(MSDCalculatorTests, ComputeDiffusionCoefficientMSD) {
  // Create a linear MSD trajectory: MSD(t) = 6.0 * t
  // Slope is 6.0, so D should be 1.0
  std::vector<double> time = {0.0, 1.0, 2.0, 3.0, 4.0};
  std::vector<double> msd = {0.0, 6.0, 12.0, 18.0, 24.0};

  double d = correlation::analysis::DynamicsAnalyzer::computeDiffusionCoefficientMSD(time, msd);
  EXPECT_NEAR(d, 1.0, 1e-6);
}

TEST(MSDCalculatorTests, DynamicsAnalyzerMSDNonPhysicalInputs) {
  // Mismatched size / empty inputs
  std::vector<double> time_empty = {};
  std::vector<double> msd_empty = {};
  EXPECT_DOUBLE_EQ(correlation::analysis::DynamicsAnalyzer::computeDiffusionCoefficientMSD(time_empty, msd_empty), 0.0);

  std::vector<double> time_small = {0.0};
  std::vector<double> msd_small = {0.0};
  EXPECT_DOUBLE_EQ(correlation::analysis::DynamicsAnalyzer::computeDiffusionCoefficientMSD(time_small, msd_small), 0.0);

  // Negative slope (non-physical diffusion)
  std::vector<double> time = {0.0, 1.0, 2.0, 3.0, 4.0};
  std::vector<double> msd_neg_slope = {24.0, 18.0, 12.0, 6.0, 0.0};
  EXPECT_DOUBLE_EQ(correlation::analysis::DynamicsAnalyzer::computeDiffusionCoefficientMSD(time, msd_neg_slope), 0.0);
}

} // namespace correlation::testing
