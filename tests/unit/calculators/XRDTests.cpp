// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "calculators/XRDCalculator.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::analysis {
namespace {
// Test fixture for XRD tests.
class XRDTests : public ::testing::Test {
protected:
  void SetUp() override {
    // A simple cubic cell containing two atoms
    cell_ = correlation::core::Cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
    cell_.addAtom("Ar", {5.0, 5.0, 5.0});
    cell_.addAtom("Ar", {6.5, 5.0, 5.0}); // Distance 1.5
  }

  void updateTrajectory() {
    trajectory_ = correlation::core::Trajectory();
    trajectory_.addFrame(cell_);
    trajectory_.precomputeBondCutoffs();
  }

public:
  correlation::core::Cell cell_;
  correlation::core::Trajectory trajectory_;
};
} // namespace

TEST_F(XRDTests, CalculateXRD) {
  updateTrajectory();
  DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // Needs RDF first
  dists.calculateRDF(5.0, 0.05);

  // Calling calculateXRD
  dists.calculateXRD(1.5406, 5.0, 90.0, 0.5);

  EXPECT_NO_THROW(dists.getHistogram("XRD"));
  const auto &hist = dists.getHistogram("XRD");
  EXPECT_FALSE(hist.bins.empty());

  // Check bins range
  EXPECT_GE(hist.bins.front(), 5.0);
  EXPECT_LE(hist.bins.back(), 90.0);
}

TEST_F(XRDTests, CalculateXRD_ThrowsIfNoRDF) {
  updateTrajectory();
  DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // Calling calculateXRD before calculateRDF should throw
  EXPECT_THROW(dists.calculateXRD(1.5406, 5.0, 90.0, 0.5), std::logic_error);
}

TEST_F(XRDTests, CalculateXRD_InvalidBinWidth) {
  updateTrajectory();
  DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  dists.calculateRDF(5.0, 0.05);

  // Calling calculateXRD with invalid bin width should throw
  EXPECT_THROW(dists.calculateXRD(1.5406, 5.0, 90.0, 0.0), std::invalid_argument);
  EXPECT_THROW(dists.calculateXRD(1.5406, 5.0, 90.0, -0.5), std::invalid_argument);
}

TEST_F(XRDTests, CalculateXRD_IntensityIsZeroAtThetaZero) {
  updateTrajectory();
  DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // Needs RDF first
  dists.calculateRDF(5.0, 0.05);

  // Calling calculateXRD with range starting at 0.0
  dists.calculateXRD(1.5406, 0.0, 10.0, 0.5);

  EXPECT_NO_THROW(dists.getHistogram("XRD"));
  const auto &hist = dists.getHistogram("XRD");
  const auto &intensities = hist.partials.at("Total");

  // Since Q(theta=0) = 0, Q < 1e-6 will branch and intensity should be exactly
  // 0.0
  EXPECT_EQ(hist.bins.front(), 0.0);
  EXPECT_DOUBLE_EQ(intensities.front(), 0.0);
}

TEST_F(XRDTests, CalculateXRDCubicCell) {
  // Simple cubic cell a = 3.0, 3x3x3 supercell
  correlation::core::Cell cubic_cell({9.0, 9.0, 9.0, 90.0, 90.0, 90.0});
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        cubic_cell.addAtom("Ar", {i * 3.0, j * 3.0, k * 3.0});
      }
    }
  }

  correlation::core::Trajectory traj;
  traj.addFrame(cubic_cell);
  traj.precomputeBondCutoffs();

  DistributionFunctions dists(cubic_cell, 4.5, traj.getBondCutoffsSQ());
  dists.calculateRDF(4.5, 0.05);

  // Calculate XRD for Cu K-alpha
  dists.calculateXRD(1.5406, 10.0, 90.0, 0.1);

  EXPECT_NO_THROW(dists.getHistogram("XRD"));
  const auto &hist = dists.getHistogram("XRD");
  EXPECT_FALSE(hist.bins.empty());

  const auto &total_intensity = hist.partials.at("Total");
  real_t max_intensity = static_cast<real_t>(0.0);
  real_t max_theta = static_cast<real_t>(0.0);
  for (size_t i = 0; i < total_intensity.size(); ++i) {
    if (hist.bins[i] > static_cast<real_t>(15.0) && total_intensity[i] > max_intensity) {
      max_intensity = total_intensity[i];
      max_theta = hist.bins[i];
    }
  }

  // We expect non-zero max intensity due to Bragg scattering
  EXPECT_GT(max_intensity, 0.0);

  // For a Simple Cubic cell with a = 3.0 A and Cu K-alpha (1.5406 A)
  // the first peak (100) theoretically appears at 2*theta ~ 29.75 degrees.
  // However, due to the small r_max (4.5 A) resulting in broad peaks, and
  // the rapidly decreasing atomic form factor for Ar, the apparent maximum
  // shifts to lower angles (~28.5 degrees).
  EXPECT_NEAR(max_theta, 28.5, 1.0);
}

TEST_F(XRDTests, CalculateXRD_InvalidInputsThrow) {
  updateTrajectory();
  DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  dists.calculateRDF(5.0, 0.05);

  EXPECT_THROW(dists.calculateXRD(0.0, 5.0, 90.0, 0.5), std::invalid_argument);
  EXPECT_THROW(dists.calculateXRD(-1.0, 5.0, 90.0, 0.5), std::invalid_argument);

  Histogram empty_gr;
  EXPECT_THROW(correlation::calculators::XRDCalculator::calculate(
                   empty_gr, cell_, {}, correlation::calculators::Wavelength{1.5406},
                   correlation::calculators::MinTheta{5.0}, correlation::calculators::MaxTheta{90.0},
                   correlation::calculators::BinWidth{0.5}),
               std::invalid_argument);
}

} // namespace correlation::analysis
