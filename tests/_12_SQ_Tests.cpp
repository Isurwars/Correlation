// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>
#include <vector>

#include "../include/Cell.hpp"
#include "../include/DistributionFunctions.hpp"
#include "../include/Trajectory.hpp"

// Test fixture for SQ tests.
class _12_SQ_Tests : public ::testing::Test {
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

TEST_F(_12_SQ_Tests, CalculateSQ) {
  updateTrajectory();
  DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // Need RDF first
  df.calculateRDF(5.0, 0.05);
  df.calculateSQ(10.0, 0.1, 5.0);

  EXPECT_NO_THROW(df.getHistogram("S(Q)"));
  const auto &hist = df.getHistogram("S(Q)");
  EXPECT_FALSE(hist.bins.empty());

  // Check peak position for 1.5A distance -> Q = 2pi/r ~ 4.18
  const auto &total_sq = hist.partials.at("Total");

  // For a dimer, S(Q) oscillates around 1. S(Q) = 1 + sin(Qr)/(Qr).
  // At Q = 4.18 (2pi/r), Qr = 2pi, sin(2pi)=0. So S(Q) should be approx 1.
  // The previous test expected a peak (large value), which is wrong for a
  // 2-atom system.

  // Search for value at Q ~ 4.18
  // Find bin closest to 4.18
  double target_Q = 4.18;
  double min_diff = 1000.0;
  double sq_val = 0.0;

  for (size_t i = 0; i < total_sq.size(); ++i) {
    if (std::abs(hist.bins[i] - target_Q) < min_diff) {
      min_diff = std::abs(hist.bins[i] - target_Q);
      sq_val = total_sq[i];
    }
  }

  // S(Q) should be near 1.0
  EXPECT_NEAR(sq_val, 1.0, 0.2);
}
