// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>
#include <vector>

#include "../include/Cell.hpp"
#include "../include/DistributionFunctions.hpp"
#include "../include/Trajectory.hpp"

// Test fixture for XRD tests.
class Test08_XRD : public ::testing::Test {
protected:
  void SetUp() override {
    // A simple cubic cell containing two atoms
    cell_ = Cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
    cell_.addAtom("Ar", {5.0, 5.0, 5.0});
    cell_.addAtom("Ar", {6.5, 5.0, 5.0}); // Distance 1.5
  }

  void updateTrajectory() {
    trajectory_ = Trajectory();
    trajectory_.addFrame(cell_);
    trajectory_.precomputeBondCutoffs();
  }

  Cell cell_{};
  Trajectory trajectory_;
};

TEST_F(Test08_XRD, CalculateXRD) {
  updateTrajectory();
  DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // Needs RDF first
  df.calculateRDF(5.0, 0.05);

  // Calling calculateXRD
  df.calculateXRD(1.5406, 5.0, 90.0, 0.5);

  EXPECT_NO_THROW(df.getHistogram("XRD"));
  const auto &hist = df.getHistogram("XRD");
  EXPECT_FALSE(hist.bins.empty());

  // Check bins range
  EXPECT_GE(hist.bins.front(), 5.0);
  EXPECT_LE(hist.bins.back(), 90.0);
}

TEST_F(Test08_XRD, CalculateXRD_ThrowsIfNoRDF) {
  updateTrajectory();
  DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // Calling calculateXRD before calculateRDF should throw
  EXPECT_THROW(df.calculateXRD(1.5406, 5.0, 90.0, 0.5), std::logic_error);
}

TEST_F(Test08_XRD, CalculateXRD_InvalidBinWidth) {
  updateTrajectory();
  DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  df.calculateRDF(5.0, 0.05);

  // Calling calculateXRD with invalid bin width should throw
  EXPECT_THROW(df.calculateXRD(1.5406, 5.0, 90.0, 0.0), std::invalid_argument);
  EXPECT_THROW(df.calculateXRD(1.5406, 5.0, 90.0, -0.5), std::invalid_argument);
}

TEST_F(Test08_XRD, CalculateXRD_IntensityIsZeroAtThetaZero) {
  updateTrajectory();
  DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // Needs RDF first
  df.calculateRDF(5.0, 0.05);

  // Calling calculateXRD with range starting at 0.0
  df.calculateXRD(1.5406, 0.0, 10.0, 0.5);

  EXPECT_NO_THROW(df.getHistogram("XRD"));
  const auto &hist = df.getHistogram("XRD");
  const auto &intensities = hist.partials.at("Total");

  // Since Q(theta=0) = 0, Q < 1e-6 will branch and intensity should be exactly
  // 0.0
  EXPECT_EQ(hist.bins.front(), 0.0);
  EXPECT_DOUBLE_EQ(intensities.front(), 0.0);
}

TEST_F(Test08_XRD, CalculateXRDCubicCell) {
  // Simple cubic cell a = 3.0, 3x3x3 supercell
  Cell cubic_cell({9.0, 9.0, 9.0, 90.0, 90.0, 90.0});
  for (int x = 0; x < 2; ++x) {
    for (int y = 0; y < 2; ++y) {
      for (int z = 0; z < 2; ++z) {
        cubic_cell.addAtom("Ar", {x * 3.0, y * 3.0, z * 3.0});
      }
    }
  }

  Trajectory traj;
  traj.addFrame(cubic_cell);
  traj.precomputeBondCutoffs();

  DistributionFunctions df(cubic_cell, 4.5, traj.getBondCutoffsSQ());
  df.calculateRDF(4.5, 0.05);

  // Calculate XRD for Cu K-alpha
  df.calculateXRD(1.5406, 10.0, 90.0, 0.1);

  EXPECT_NO_THROW(df.getHistogram("XRD"));
  const auto &hist = df.getHistogram("XRD");
  EXPECT_FALSE(hist.bins.empty());

  const auto &total_intensity = hist.partials.at("Total");
  double max_intensity = 0.0;
  for (double val : total_intensity) {
    if (val > max_intensity) {
      max_intensity = val;
    }
  }

  // We expect non-zero max intensity due to Bragg scattering
  EXPECT_GT(max_intensity, 0.0);
}
