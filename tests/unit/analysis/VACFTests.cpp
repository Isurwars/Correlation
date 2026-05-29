// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "analysis/DynamicsAnalyzer.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::analysis {

// Test fixture for VACF and VDOS tests.
class VACFTests : public ::testing::Test {
protected:
  // No special setup needed for VACF usually, or different from Test06
};

TEST_F(VACFTests, CalculateVACF_and_VDOS) {
  correlation::core::Cell c({10, 10, 10, 90, 90, 90});
  c.addAtom("Ar", {0, 0, 0});
  correlation::core::Trajectory t;
  t.addFrame(c);
  t.addFrame(c); // Static
  t.calculateVelocities();

  DistributionFunctions df(c, 0.0, {{0.0}});

  df.calculateVACF(t, 1);
  EXPECT_NO_THROW(df.getHistogram("VACF"));
  const auto &vacf = df.getHistogram("VACF").partials.at("Total");
  // Should be 1.0 normalized (if constant 0 velocity? Wait 0 velocity -> ??)
  // If static, position constant -> velocity 0. Correlation of 0 with 0 is 0.
  // Let's give it velocity.

  correlation::core::Trajectory tMoving;
  correlation::core::Cell c1 = c;

  correlation::core::Cell c2({10, 10, 10, 90, 90, 90});
  // correlation::core::Atom 1 moves +1.0 in x
  c2.addAtom("Ar", {1.0, 0.0, 0.0});
  // correlation::core::Atom 2 moves -1.0 in x (balancing COM)
  c2.addAtom("Ar", {-1.0, 0.0, 0.0});

  correlation::core::Cell c3({10, 10, 10, 90, 90, 90});
  // correlation::core::Atom 1 moves to +2.0
  c3.addAtom("Ar", {2.0, 0.0, 0.0});
  // correlation::core::Atom 2 moves to -2.0
  c3.addAtom("Ar", {-2.0, 0.0, 0.0});

  // Need to update c1 (frame 0) to have 2 atoms at 0
  c1 = correlation::core::Cell({10, 10, 10, 90, 90, 90});
  c1.addAtom("Ar", {0.0, 0.0, 0.0});
  c1.addAtom("Ar", {0.0, 0.0, 0.0});

  tMoving.addFrame(c1);
  tMoving.addFrame(c2);
  tMoving.addFrame(c3);
  tMoving.setTimeStep(1.0);
  tMoving.calculateVelocities();

  df.calculateVACF(tMoving, 1);
  const auto &vacf2 = df.getHistogram("Normalized VACF").partials.at("Total");
  EXPECT_NEAR(vacf2[0], 1.0, 1e-6); // t=0
  EXPECT_NEAR(vacf2[1], 1.0, 1e-6); // t=1, const velocity

  // VDOS
  df.calculateVDOS();
  EXPECT_NO_THROW(df.getHistogram("VDOS"));
  const auto &vdos_hist = df.getHistogram("VDOS");
  EXPECT_TRUE(vdos_hist.partials.count("Frequency_cm_1"));
}

TEST_F(VACFTests, CalculateVACF_WithFrameRange) {
  correlation::core::Trajectory tRange;
  tRange.setTimeStep(1.0);

  // Create 10 frames
  for (int i = 0; i < 10; ++i) {
    correlation::core::Cell c({10, 10, 10, 90, 90, 90});
    c.addAtom("Ar", {static_cast<double>(i), 0.0, 0.0});
    c.addAtom("Ar", {-static_cast<double>(i), 0.0, 0.0});
    tRange.addFrame(c);
  }
  tRange.calculateVelocities();

  correlation::core::Cell base_cell({10, 10, 10, 90, 90, 90});
  base_cell.addAtom("Ar", {0, 0, 0});
  base_cell.addAtom("Ar", {0, 0, 0});
  DistributionFunctions df(base_cell, 0.0, {{0.0}});

  // Restrict VACF to frames 2 through 7 (5 frames total)
  df.calculateVACF(tRange, 3, 2, 7);

  EXPECT_NO_THROW(df.getHistogram("VACF"));
  const auto &vacf = df.getHistogram("VACF").partials.at("Total");

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

  std::vector<double> positions = {0.0, 1.0, 1.75, 2.25, 2.5, 2.625, 2.6875, 2.71875, 2.734375};

  for (double x : positions) {
    correlation::core::Cell c({10, 10, 10, 90, 90, 90});
    c.addAtom("Ar", {x, 0.0, 0.0});
    c.addAtom("Ar", {-x, 0.0, 0.0}); // To balance COM
    tGas.addFrame(c);
  }
  tGas.calculateVelocities();

  correlation::core::Cell base_cell({10, 10, 10, 90, 90, 90});
  base_cell.addAtom("Ar", {0, 0, 0});
  base_cell.addAtom("Ar", {0, 0, 0});
  DistributionFunctions df(base_cell, 0.0, {{0.0}});

  df.calculateVACF(tGas, 4); // Calculate up to 4 lags
  EXPECT_NO_THROW(df.getHistogram("Normalized VACF"));

  const auto &vacf = df.getHistogram("Normalized VACF").partials.at("Total");

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
  std::vector<double> time = {0.0, 1.0, 2.0};
  std::vector<double> vacf = {3.0, 3.0, 3.0};
  std::vector<double> norm_vacf = {1.0, 1.0, 1.0};

  double d = DynamicsAnalyzer::computeDiffusionCoefficientVACF(time, vacf);
  double tau = DynamicsAnalyzer::computeRelaxationTime(time, norm_vacf);

  EXPECT_NEAR(d, 2.0, 1e-6);
  EXPECT_NEAR(tau, 2.0, 1e-6);
}

TEST_F(VACFTests, DistributionFunctionsDynamicProperties) {
  correlation::core::Cell c({10, 10, 10, 90, 90, 90});
  c.addAtom("Ar", {0.0, 0.0, 0.0});
  DistributionFunctions df(c, 0.0, {{0.0}});

  // Verify initial state is 0.0
  EXPECT_DOUBLE_EQ(df.getDiffusionCoefficientMSD(), 0.0);
  EXPECT_DOUBLE_EQ(df.getDiffusionCoefficientVACF(), 0.0);
  EXPECT_DOUBLE_EQ(df.getRelaxationTime(), 0.0);
  EXPECT_DOUBLE_EQ(df.getDeborahNumber(), 0.0);

  // Set values
  df.setDiffusionCoefficientMSD(1.23);
  df.setDiffusionCoefficientVACF(4.56);
  df.setRelaxationTime(7.89);
  df.setDeborahNumber(0.12);

  // Verify updated values
  EXPECT_DOUBLE_EQ(df.getDiffusionCoefficientMSD(), 1.23);
  EXPECT_DOUBLE_EQ(df.getDiffusionCoefficientVACF(), 4.56);
  EXPECT_DOUBLE_EQ(df.getRelaxationTime(), 7.89);
  EXPECT_DOUBLE_EQ(df.getDeborahNumber(), 0.12);
}

TEST_F(VACFTests, DynamicsAnalyzerNonPhysicalInputs) {
  // Test mismatched size / empty inputs
  std::vector<double> time_empty = {};
  std::vector<double> vacf_empty = {};
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeDiffusionCoefficientVACF(time_empty, vacf_empty), 0.0);
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeRelaxationTime(time_empty, vacf_empty), 0.0);

  std::vector<double> time_small = {0.0};
  std::vector<double> vacf_small = {1.0};
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeDiffusionCoefficientVACF(time_small, vacf_small), 0.0);
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeRelaxationTime(time_small, vacf_small), 0.0);

  std::vector<double> time_mismatch = {0.0, 1.0};
  std::vector<double> vacf_mismatch = {1.0, 1.0, 1.0};
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeDiffusionCoefficientVACF(time_mismatch, vacf_mismatch), 0.0);
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeRelaxationTime(time_mismatch, vacf_mismatch), 0.0);

  // Test non-increasing time values (dt <= 0)
  std::vector<double> time_non_inc = {0.0, 0.0, 1.0};
  std::vector<double> vacf_valid = {1.0, 1.0, 1.0};
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeDiffusionCoefficientVACF(time_non_inc, vacf_valid), 0.0);
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeRelaxationTime(time_non_inc, vacf_valid), 0.0);

  // Test negative result handling (Green-Kubo integral < 0)
  std::vector<double> time_valid = {0.0, 1.0, 2.0};
  std::vector<double> vacf_neg = {-5.0, -10.0, -5.0};
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeDiffusionCoefficientVACF(time_valid, vacf_neg), 0.0);
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeRelaxationTime(time_valid, vacf_neg), 0.0);
}

TEST_F(VACFTests, DistributionFunctionsNonPhysicalOptions) {
  correlation::core::Cell c({10, 10, 10, 90, 90, 90});
  c.addAtom("Ar", {0.0, 0.0, 0.0});
  DistributionFunctions df(c, 0.0, {{0.0}});

  // calculateRDF guards
  EXPECT_THROW(df.calculateRDF(-5.0, 0.05), std::invalid_argument);
  EXPECT_THROW(df.calculateRDF(20.0, -0.05), std::invalid_argument);
  EXPECT_THROW(df.calculateRDF(0.0, 0.05), std::invalid_argument);
  EXPECT_THROW(df.calculateRDF(20.0, 0.0), std::invalid_argument);

  // calculatePAD guards
  EXPECT_THROW(df.calculatePAD(-1.0), std::invalid_argument);
  EXPECT_THROW(df.calculatePAD(0.0), std::invalid_argument);

  // calculateDAD guards
  EXPECT_THROW(df.calculateDAD(-1.0), std::invalid_argument);
  EXPECT_THROW(df.calculateDAD(0.0), std::invalid_argument);

  // calculateXRD guards
  EXPECT_THROW(df.calculateXRD(-1.0, 5.0, 90.0, 1.0), std::invalid_argument);
  EXPECT_THROW(df.calculateXRD(1.54, -5.0, 90.0, 1.0), std::invalid_argument);
  EXPECT_THROW(df.calculateXRD(1.54, 5.0, -90.0, 1.0), std::invalid_argument);
  EXPECT_THROW(df.calculateXRD(1.54, 90.0, 5.0, 1.0), std::invalid_argument);
  EXPECT_THROW(df.calculateXRD(1.54, 5.0, 90.0, -1.0), std::invalid_argument);
}

} // namespace correlation::analysis
