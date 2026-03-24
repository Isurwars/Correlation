// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>
#include <cmath>

#include "../include/Cell.hpp"
#include "../include/DistributionFunctions.hpp"
#include "../include/calculators/StructureFactorCalculator.hpp"

class _23_StructureFactorCalculator_Tests : public ::testing::Test {
protected:
  void SetUp() override {}
};

// For a simple cubic crystal, S(Q) should show peaks at the allowed Bragg
// reflections. The first peak for a simple cubic lattice is at Q = 2*pi/a.
TEST_F(_23_StructureFactorCalculator_Tests, CalculatesSimpleCubicBraggPeak) {
  // Build a 2x2x2 simple cubic supercell with lattice constant a=3.0 Å
  const double a = 3.0;
  Cell cell(std::array<double, 6>{2*a, 2*a, 2*a, 90.0, 90.0, 90.0});

  // 8 atoms at corners of a 2x2x2 supercell
  for (int ix = 0; ix < 2; ++ix)
    for (int iy = 0; iy < 2; ++iy)
      for (int iz = 0; iz < 2; ++iz)
        cell.addAtom("Si", {ix * a, iy * a, iz * a});

  DistributionFunctions df(cell);
  StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 5.0;
  settings.q_bin_width = 0.05;

  calc.calculateFrame(df, settings);

  const auto &hist = df.getHistogram("S(Q)_pw");
  EXPECT_FALSE(hist.bins.empty());
  EXPECT_TRUE(hist.partials.count("Total"));

  // The first Bragg peak for a SC lattice at a=3 is at Q = 2*pi/3 ≈ 2.094 Å^{-1}
  // Find the max S(Q) value in the range [1.8, 2.4]
  double peak_sq = 0.0;
  for (size_t i = 0; i < hist.bins.size(); ++i) {
    if (hist.bins[i] >= 1.8 && hist.bins[i] <= 2.4) {
      peak_sq = std::max(peak_sq, hist.partials.at("Total")[i]);
    }
  }
  // At a Bragg peak, S(Q) should be much greater than 1
  EXPECT_GT(peak_sq, 2.0);
}

// For a single atom, S(Q) = 1 for all Q.
TEST_F(_23_StructureFactorCalculator_Tests, SingleAtomGivesOne) {
  Cell cell(std::array<double, 6>{10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  DistributionFunctions df(cell);
  StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 5.0;
  settings.q_bin_width = 0.1;

  calc.calculateFrame(df, settings);
  const auto &hist = df.getHistogram("S(Q)_pw");

  ASSERT_FALSE(hist.bins.empty());
  for (size_t i = 0; i < hist.bins.size(); ++i) {
    // S(Q) = 1/N * |sum exp(i q.r)|^2 = 1/1 * |exp(i*0)|^2 = 1.0
    if (hist.partials.at("Total")[i] > 0.0) {
      EXPECT_NEAR(hist.partials.at("Total")[i], 1.0, 1e-6);
    }
  }
}
