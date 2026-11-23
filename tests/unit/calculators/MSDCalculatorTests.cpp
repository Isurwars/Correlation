// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DynamicsAnalyzer.hpp"
#include "calculators/temporal/MSDCalculator.hpp"
#include "core/Trajectory.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {

using namespace correlation::calculators;
using namespace correlation::core;
using namespace correlation::analysis;
using correlation::real_t;

namespace {
void verifyMSDValues(const Histogram &msd_hist, const Histogram &deff_hist) {
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

void verifyQuadraticMsdCurve(const std::vector<real_t> &msd, double lag1_val, double lag2_val) {
  ASSERT_EQ(msd.size(), 3U);
  EXPECT_NEAR(msd[0], 0.0, 1e-6);
  EXPECT_NEAR(msd[1], lag1_val, 1e-5);
  EXPECT_NEAR(msd[2], lag2_val, 1e-5);
}

void verifyPartialMsdMetadata(const Histogram &msd_hist, const Histogram &deff_hist) {
  EXPECT_EQ(msd_hist.file_suffix, "_MSD");
  EXPECT_EQ(deff_hist.file_suffix, "_D_eff");
  ASSERT_TRUE(msd_hist.partials.contains("Si"));
  ASSERT_TRUE(msd_hist.partials.contains("O"));
  ASSERT_TRUE(deff_hist.partials.contains("Si"));
  ASSERT_TRUE(deff_hist.partials.contains("O"));
}

void verifyPartialMsdLag1(const Histogram &msd_hist, const Histogram &deff_hist) {
  EXPECT_NEAR(msd_hist.partials.at("Si")[1], 0.25, 1e-5);
  EXPECT_NEAR(msd_hist.partials.at("O")[1], 1.0, 1e-5);
  EXPECT_NEAR(msd_hist.partials.at("Total")[1], 0.625, 1e-5);
  EXPECT_NEAR(deff_hist.partials.at("Si")[1], 0.25 / 6.0, 1e-5);
  EXPECT_NEAR(deff_hist.partials.at("O")[1], 1.0 / 6.0, 1e-5);
}
} // namespace

TEST(MSDCalculatorTests, ComputesCorrectMSDAndDeff) {
  // Construct a trajectory of 5 frames with 1 atom moving linearly
  // x(t) = t, y(t) = 0, z(t) = 0
  std::vector<Cell> frames;
  for (int i = 0; i < 5; ++i) {
    Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
    cell.addAtom("Si", {static_cast<real_t>(i), 0.0, 0.0});
    frames.push_back(std::move(cell));
  }

  const real_t time_step = 1.0; // 1 fs
  const Trajectory traj(std::move(frames), time_step);

  // Act
  const auto results = MSDCalculator::calculate(traj, MaxFrames{2}, StartFrame{0}, EndFrame{5});

  // Assert
  ASSERT_TRUE(results.contains("MSD"));
  ASSERT_TRUE(results.contains("D_eff"));

  const auto &msd_hist = results.at("MSD");
  const auto &deff_hist = results.at("D_eff");

  verifyMSDValues(msd_hist, deff_hist);
}

// --- Extreme / Edge-Case Tests ---

TEST(MSDCalculatorTests, StationaryAtomsMSDIsZero) {
  // All atoms remain stationary — MSD should be 0 at all lags
  // Use addFrame() to bypass the constructor's automatic deduplication.
  Trajectory traj;
  traj.setTimeStep(1.0);
  for (int i = 0; i < 5; ++i) {
    Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
    cell.addAtom("Ar", {5.0, 5.0, 5.0}); // Same position every frame
    traj.addFrame(cell);
  }

  auto results = MSDCalculator::calculate(traj, MaxFrames{3}, StartFrame{0}, EndFrame{5});

  ASSERT_TRUE(results.contains("MSD"));
  const auto &msd = results.at("MSD").partials.at("Total");

  for (size_t i = 0; i < msd.size(); ++i) {
    EXPECT_NEAR(msd[i], 0.0, 1e-10) << "MSD should be 0 at lag " << i;
  }
}

TEST(MSDCalculatorTests, MultiAtomAveraging) {
  // Two equal-mass atoms: one moves +x, one moves +y.
  // With default CoM drift removal, v_com = (+0.5, +0.5, 0), so relative displacements are
  // (+0.5, -0.5, 0) and (-0.5, +0.5, 0) -> squared displacement = 0.5 * lag^2.
  std::vector<Cell> frames;
  for (int i = 0; i < 4; ++i) {
    Cell cell({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0});
    cell.addAtom("Si", {static_cast<real_t>(i), 5.0, 5.0}); // Moves +x
    cell.addAtom("Si", {5.0, static_cast<real_t>(i), 5.0}); // Moves +y
    frames.push_back(std::move(cell));
  }

  const Trajectory traj(std::move(frames), 1.0);
  const auto results = MSDCalculator::calculate(traj, MaxFrames{2}, StartFrame{0}, EndFrame{4});

  ASSERT_TRUE(results.contains("MSD"));
  verifyQuadraticMsdCurve(results.at("MSD").partials.at("Total"), 0.5, 2.0);

  // Without CoM drift removal, lab-frame squared displacement = 1.0 * lag^2.
  const auto lab_msd =
      DynamicsAnalyzer::calculateMSD(traj, MaxFrames{2}, StartFrame{0}, EndFrame{4}, false);
  verifyQuadraticMsdCurve(lab_msd, 1.0, 4.0);
}

TEST(MSDCalculatorTests, CoMDriftRemovalEliminatesRigidTranslation) {
  // Two atoms translating rigidly together at v = (1, 0, 0): CoM-corrected MSD must be 0.
  std::vector<Cell> frames;
  for (int i = 0; i < 6; ++i) {
    Cell cell({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0});
    cell.addAtom("Si", {static_cast<real_t>(i), 2.0, 2.0});
    cell.addAtom("O", {static_cast<real_t>(i) + static_cast<real_t>(2.0), 2.0, 2.0});
    frames.push_back(std::move(cell));
  }

  const Trajectory traj(std::move(frames), 1.0);
  const auto results = MSDCalculator::calculate(traj, MaxFrames{3}, StartFrame{0}, EndFrame{6});
  ASSERT_TRUE(results.contains("MSD"));
  for (const real_t val : results.at("MSD").partials.at("Total")) {
    EXPECT_NEAR(val, 0.0, 1e-6);
  }
}

TEST(MSDCalculatorTests, PerElementPartialMSDAndDeffFileSuffix) {
  // Counter-moving pairs of Si (±0.5 Å/fs) and O (±1.0 Å/fs) have zero net CoM drift.
  std::vector<Cell> frames;
  for (int i = 0; i < 6; ++i) {
    const auto step = static_cast<real_t>(i);
    Cell cell({40.0, 0.0, 0.0}, {0.0, 40.0, 0.0}, {0.0, 0.0, 40.0});
    cell.addAtom("Si", {static_cast<real_t>(10.0) + (static_cast<real_t>(0.5) * step), 5.0, 5.0});
    cell.addAtom("Si", {static_cast<real_t>(10.0) - (static_cast<real_t>(0.5) * step), 5.0, 5.0});
    cell.addAtom("O", {static_cast<real_t>(20.0) + step, 15.0, 5.0});
    cell.addAtom("O", {static_cast<real_t>(20.0) - step, 15.0, 5.0});
    frames.push_back(std::move(cell));
  }

  const Trajectory traj(std::move(frames), 1.0);
  const auto results = MSDCalculator::calculate(traj, MaxFrames{2}, StartFrame{0}, EndFrame{6});
  ASSERT_TRUE(results.contains("MSD"));
  ASSERT_TRUE(results.contains("D_eff"));

  const auto &msd_hist = results.at("MSD");
  const auto &deff_hist = results.at("D_eff");
  verifyPartialMsdMetadata(msd_hist, deff_hist);
  verifyPartialMsdLag1(msd_hist, deff_hist);
}

TEST(MSDCalculatorTests, CenteredOLSPreservesPrecisionOnLargeTimeOffsets) {
  // Large baseline time offset where uncentered float OLS suffers catastrophic cancellation.
  std::vector<real_t> time(200);
  std::vector<real_t> msd(200);
  for (size_t i = 0; i < 200; ++i) {
    const double t_val = 100000.0 + static_cast<double>(i);
    time[i] = static_cast<real_t>(t_val);
    msd[i] = static_cast<real_t>(6.0 * 0.25 * static_cast<double>(i));
  }
  const real_t d_coeff = DynamicsAnalyzer::computeDiffusionCoefficientMSD(time, msd);
  EXPECT_NEAR(d_coeff, 0.25, 1e-3);
}

TEST(MSDCalculatorTests, FrameRangeSubset) {
  // Create 10 frames, but only compute MSD over frames 3-7
  std::vector<Cell> frames;
  for (int i = 0; i < 10; ++i) {
    Cell cell({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0});
    cell.addAtom("Si", {static_cast<real_t>(i), 0.0, 0.0});
    frames.push_back(std::move(cell));
  }

  const Trajectory traj(std::move(frames), 1.0);
  // Restrict to frames [3, 7) with max correlation lag 2
  auto results = MSDCalculator::calculate(traj, MaxFrames{2}, StartFrame{3}, EndFrame{7});

  ASSERT_TRUE(results.contains("MSD"));
  const auto &msd = results.at("MSD").partials.at("Total");

  // Lag 0 should always be 0
  EXPECT_NEAR(msd[0], 0.0, 1e-6);
  // Lag 1: displacement is always 1.0 (uniform motion) -> MSD = 1.0
  EXPECT_NEAR(msd[1], 1.0, 1e-6);
}

TEST(MSDCalculatorTests, ComputeDiffusionCoefficientMSD) {
  // Create a linear MSD trajectory: MSD(t) = 6.0 * t
  // Slope is 6.0, so D should be 1.0
  std::vector<real_t> const time = {0.0, 1.0, 2.0, 3.0, 4.0};
  std::vector<real_t> const msd = {0.0, 6.0, 12.0, 18.0, 24.0};

  real_t const diffusion_coefficient =
      correlation::analysis::DynamicsAnalyzer::computeDiffusionCoefficientMSD(time, msd);
  EXPECT_NEAR(diffusion_coefficient, 1.0, 1e-6);
}

TEST(MSDCalculatorTests, DynamicsAnalyzerMSDNonPhysicalInputs) {
  // Mismatched size / empty inputs
  const std::vector<real_t> time_empty = {};
  const std::vector<real_t> msd_empty = {};
  EXPECT_DOUBLE_EQ(correlation::analysis::DynamicsAnalyzer::computeDiffusionCoefficientMSD(
                       time_empty, msd_empty),
                   0.0);

  const std::vector<real_t> time_small = {0.0};
  const std::vector<real_t> msd_small = {0.0};
  EXPECT_DOUBLE_EQ(correlation::analysis::DynamicsAnalyzer::computeDiffusionCoefficientMSD(
                       time_small, msd_small),
                   0.0);

  // Negative slope (non-physical diffusion)
  const std::vector<real_t> time = {0.0, 1.0, 2.0, 3.0, 4.0};
  const std::vector<real_t> msd_neg_slope = {24.0, 18.0, 12.0, 6.0, 0.0};
  EXPECT_DOUBLE_EQ(
      correlation::analysis::DynamicsAnalyzer::computeDiffusionCoefficientMSD(time, msd_neg_slope),
      0.0);
}

} // namespace correlation::testing
