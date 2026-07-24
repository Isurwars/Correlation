// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "analysis/DynamicsAnalyzer.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::analysis {
namespace {
// Test fixture for VACF and VDOS tests.
class VACFTests : public ::testing::Test {
protected:
  // No special setup needed for VACF usually, or different from Test06
};
} // namespace

TEST_F(VACFTests, CalculateVACF_and_VDOS) {
  correlation::core::Cell cell({10, 10, 10, 90, 90, 90});
  cell.addAtom("Ar", {0, 0, 0});
  correlation::core::Trajectory trajectory;
  trajectory.addFrame(cell);
  trajectory.addFrame(cell); // Static
  trajectory.calculateVelocities();

  DistributionFunctions dists(cell, 0.0, {{0.0}});

  dists.calculateVACF(trajectory, correlation::analysis::MaxFrames{1});
  EXPECT_NO_THROW(dists.getHistogram("VACF"));
  const auto &vacf = dists.getHistogram("VACF").partials.at("Total");

  correlation::core::Trajectory tMoving;
  correlation::core::Cell c_1({10, 10, 10, 90, 90, 90});
  c_1.addAtom("Ar", {0.0, 0.0, 0.0});
  c_1.addAtom("Ar", {0.0, 0.0, 0.0});

  correlation::core::Cell c_2({10, 10, 10, 90, 90, 90});
  c_2.addAtom("Ar", {1.0, 0.0, 0.0});
  c_2.addAtom("Ar", {-1.0, 0.0, 0.0});

  correlation::core::Cell c_3({10, 10, 10, 90, 90, 90});
  c_3.addAtom("Ar", {2.0, 0.0, 0.0});
  c_3.addAtom("Ar", {-2.0, 0.0, 0.0});

  tMoving.addFrame(c_1);
  tMoving.addFrame(c_2);
  tMoving.addFrame(c_3);
  tMoving.setTimeStep(1.0);
  tMoving.calculateVelocities();

  dists.calculateVACF(tMoving, correlation::analysis::MaxFrames{1});
  const auto &vacf2 = dists.getHistogram("Normalized VACF").partials.at("Total");
  EXPECT_NEAR(vacf2[0], 1.0, 1e-6); // t=0
  EXPECT_NEAR(vacf2[1], 1.0, 1e-6); // t=1, const velocity

  // VDOS
  dists.calculateVDOS();
  EXPECT_NO_THROW(dists.getHistogram("VDOS"));
  const auto &vdos_hist = dists.getHistogram("VDOS");
  EXPECT_TRUE(vdos_hist.partials.count("Frequency_cm_1"));
}

TEST_F(VACFTests, CalculateVACF_WithFrameRange) {
  correlation::core::Trajectory tRange;
  tRange.setTimeStep(1.0);

  for (int i = 0; i < 10; ++i) {
    correlation::core::Cell cell({10, 10, 10, 90, 90, 90});
    cell.addAtom("Ar", {static_cast<real_t>(i), 0.0, 0.0});
    cell.addAtom("Ar", {-static_cast<real_t>(i), 0.0, 0.0});
    tRange.addFrame(cell);
  }
  tRange.calculateVelocities();

  correlation::core::Cell base_cell({10, 10, 10, 90, 90, 90});
  base_cell.addAtom("Ar", {0, 0, 0});
  base_cell.addAtom("Ar", {0, 0, 0});
  DistributionFunctions dists(base_cell, 0.0, {{0.0}});

  // Restrict VACF to frames 2 through 7 (5 frames total)
  dists.calculateVACF(tRange, correlation::analysis::MaxFrames{3}, correlation::analysis::StartFrame{2},
                      correlation::analysis::EndFrame{7});

  EXPECT_NO_THROW(dists.getHistogram("VACF"));
  const auto &vacf = dists.getHistogram("VACF").partials.at("Total");

  // Since we requested max correlation frames = 3, size should be 4 (includes
  // lag 0)
  EXPECT_EQ(vacf.size(), 4);

  // Velocity is constant (+1.0 and -1.0), so dot product is (1*1) + (-1*-1)
  // = 2.0 per frame pair. Averaged over 2 atoms, vacf[lag] should be 1.0 for
  // all lags within the constrained region.
  EXPECT_NEAR(vacf[0], 1.0, 1e-6);
  EXPECT_NEAR(vacf[1], 1.0, 1e-6);
  EXPECT_NEAR(vacf[2], 1.0, 1e-6);
  EXPECT_NEAR(vacf[3], 1.0, 1e-6);
}

TEST_F(VACFTests, CalculateVACF_GasLike) {
  correlation::core::Trajectory tGas;
  tGas.setTimeStep(1.0);

  // We want to simulate a gas-like behavior where VACF decays exponentially
  // without negative regions. We'll provide a sequence of positions such that
  // the central-difference velocity decreases strictly. r(0) = 0.0 r(1) = 1.0
  // r(2) = 1.75
  // r(3) = 2.25
  // r(4) = 2.5
  // r(5) = 2.625
  // r(6) = 2.6875
  // r(7) = 2.71875
  // r(8) = 2.734375

  std::vector<real_t> positions = {0.0, 1.0, 1.75, 2.25, 2.5, 2.625, 2.6875, 2.71875, 2.734375};

  for (real_t pos : positions) {
    correlation::core::Cell cell({10, 10, 10, 90, 90, 90});
    cell.addAtom("Ar", {pos, 0.0, 0.0});
    cell.addAtom("Ar", {-pos, 0.0, 0.0}); // To balance COM
    tGas.addFrame(cell);
  }
  tGas.calculateVelocities();

  correlation::core::Cell base_cell({10, 10, 10, 90, 90, 90});
  base_cell.addAtom("Ar", {0, 0, 0});
  base_cell.addAtom("Ar", {0, 0, 0});
  DistributionFunctions dists(base_cell, 0.0, {{0.0}});

  dists.calculateVACF(tGas, correlation::analysis::MaxFrames{4}); // Calculate up to 4 lags
  EXPECT_NO_THROW(dists.getHistogram("Normalized VACF"));

  const auto &vacf = dists.getHistogram("Normalized VACF").partials.at("Total");

  // Since velocities strictly decrease, the normalized VACF should strictly
  // decrease but remain positive (gas-like monotonic decay).
  EXPECT_GT(vacf.size(), 4);
  EXPECT_NEAR(vacf[0], 1.0, 1e-6);
  EXPECT_LT(vacf[1], vacf[0]);
  EXPECT_LT(vacf[2], vacf[1]);
  EXPECT_LT(vacf[3], vacf[2]);
  EXPECT_LT(vacf[4], vacf[3]);
  EXPECT_GT(vacf[4], 0.0);
}

TEST_F(VACFTests, ComputeDiffusionCoefficientVACF_and_RelaxationTime) {
  std::vector<real_t> time = {0.0, 1.0, 2.0};
  std::vector<real_t> vacf = {3.0, 3.0, 3.0};
  std::vector<real_t> norm_vacf = {1.0, 1.0, 1.0};

  real_t vacf_diffusion = DynamicsAnalyzer::computeDiffusionCoefficientVACF(time, vacf);
  real_t relaxation_time = DynamicsAnalyzer::computeRelaxationTime(time, norm_vacf);

  EXPECT_NEAR(vacf_diffusion, 2.0, 1e-6);
  EXPECT_NEAR(relaxation_time, 2.0, 1e-6);
}

TEST_F(VACFTests, DistributionFunctionsDynamicProperties) {
  correlation::core::Cell cell({10, 10, 10, 90, 90, 90});
  cell.addAtom("Ar", {0.0, 0.0, 0.0});
  DistributionFunctions dists(cell, 0.0, {{0.0}});

  // Verify initial state is 0.0
  EXPECT_DOUBLE_EQ(dists.getDiffusionCoefficientMSD(), 0.0);
  EXPECT_DOUBLE_EQ(dists.getDiffusionCoefficientVACF(), 0.0);
  EXPECT_DOUBLE_EQ(dists.getRelaxationTime(), 0.0);
  EXPECT_DOUBLE_EQ(dists.getDeborahNumber(), 0.0);

  // Set values
  dists.setDiffusionCoefficientMSD(1.23);
  dists.setDiffusionCoefficientVACF(4.56);
  dists.setRelaxationTime(7.89);
  dists.setDeborahNumber(0.12);

  // Verify updated values
  EXPECT_THAT(dists.getDiffusionCoefficientMSD(), correlation::testing::IsRealEq(1.23));
  EXPECT_THAT(dists.getDiffusionCoefficientVACF(), correlation::testing::IsRealEq(4.56));
  EXPECT_THAT(dists.getRelaxationTime(), correlation::testing::IsRealEq(7.89));
  EXPECT_THAT(dists.getDeborahNumber(), correlation::testing::IsRealEq(0.12));
}

TEST_F(VACFTests, DynamicsAnalyzerNonPhysicalInputs) {
  // Test mismatched size / empty inputs
  std::vector<real_t> time_empty = {};
  std::vector<real_t> vacf_empty = {};
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeDiffusionCoefficientVACF(time_empty, vacf_empty), 0.0);
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeRelaxationTime(time_empty, vacf_empty), 0.0);

  std::vector<real_t> time_small = {0.0};
  std::vector<real_t> vacf_small = {1.0};
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeDiffusionCoefficientVACF(time_small, vacf_small), 0.0);
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeRelaxationTime(time_small, vacf_small), 0.0);

  std::vector<real_t> time_mismatch = {0.0, 1.0};
  std::vector<real_t> vacf_mismatch = {1.0, 1.0, 1.0};
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeDiffusionCoefficientVACF(time_mismatch, vacf_mismatch), 0.0);
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeRelaxationTime(time_mismatch, vacf_mismatch), 0.0);

  // Test non-increasing time values (dt <= 0)
  std::vector<real_t> time_non_inc = {0.0, 0.0, 1.0};
  std::vector<real_t> vacf_valid = {1.0, 1.0, 1.0};
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeDiffusionCoefficientVACF(time_non_inc, vacf_valid), 0.0);
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeRelaxationTime(time_non_inc, vacf_valid), 0.0);

  // Test negative result handling (Green-Kubo integral < 0)
  std::vector<real_t> time_valid = {0.0, 1.0, 2.0};
  std::vector<real_t> vacf_neg = {-5.0, -10.0, -5.0};
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeDiffusionCoefficientVACF(time_valid, vacf_neg), 0.0);
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeRelaxationTime(time_valid, vacf_neg), 0.0);
}

TEST_F(VACFTests, DistributionFunctionsNonPhysicalOptions) {
  correlation::core::Cell cell({10, 10, 10, 90, 90, 90});
  cell.addAtom("Ar", {0.0, 0.0, 0.0});
  DistributionFunctions dists(cell, 0.0, {{0.0}});

  // calculateRDF guards
  EXPECT_THROW(dists.calculateRDF({.r_max = -5.0, .r_bin_width = 0.05}), std::invalid_argument);
  EXPECT_THROW(dists.calculateRDF({.r_max = 20.0, .r_bin_width = -0.05}), std::invalid_argument);
  EXPECT_THROW(dists.calculateRDF({.r_max = 0.0, .r_bin_width = 0.05}), std::invalid_argument);
  EXPECT_THROW(dists.calculateRDF({.r_max = 20.0, .r_bin_width = 0.0}), std::invalid_argument);

  // calculatePAD guards
  EXPECT_THROW(dists.calculatePAD(-1.0), std::invalid_argument);
  EXPECT_THROW(dists.calculatePAD(0.0), std::invalid_argument);

  // calculateDAD guards
  EXPECT_THROW(dists.calculateDAD(-1.0), std::invalid_argument);
  EXPECT_THROW(dists.calculateDAD(0.0), std::invalid_argument);

  // calculateXRD guards
  EXPECT_THROW(dists.calculateXRD({.lambda = -1.0, .theta_min = 5.0, .theta_max = 90.0, .bin_width = 1.0}),
               std::invalid_argument);
  EXPECT_THROW(dists.calculateXRD({.lambda = 1.54, .theta_min = -5.0, .theta_max = 90.0, .bin_width = 1.0}),
               std::invalid_argument);
  EXPECT_THROW(dists.calculateXRD({.lambda = 1.54, .theta_min = 5.0, .theta_max = -90.0, .bin_width = 1.0}),
               std::invalid_argument);
  EXPECT_THROW(dists.calculateXRD({.lambda = 1.54, .theta_min = 90.0, .theta_max = 5.0, .bin_width = 1.0}),
               std::invalid_argument);
  EXPECT_THROW(dists.calculateXRD({.lambda = 1.54, .theta_min = 5.0, .theta_max = 90.0, .bin_width = -1.0}),
               std::invalid_argument);
}

} // namespace correlation::analysis
