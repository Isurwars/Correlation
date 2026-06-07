// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "calculators/SDFCalculator.hpp"
#include <gtest/gtest.h>

using namespace correlation::calculators;
using namespace correlation::analysis;

TEST(SDFCalculatorTests, CalculateSDF) {
  correlation::core::Cell cell;
  cell.setLatticeParameters({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});

  DistributionFunctions df(cell);
  AnalysisSettings settings;
  settings.r_bin_width = 1.0;

  SDFCalculator calc;
  calc.calculateFrame(df, settings);

  const auto &hist = df.getHistogram("SDF");
  EXPECT_EQ(hist.bins.size(), 1000); // 10x10x10 grid
  EXPECT_TRUE(hist.partials.count("Ar") > 0);
  EXPECT_TRUE(hist.partials.count("Total") > 0);
}

// --- Strengthened & Edge-Case Tests ---

TEST(SDFCalculatorTests, MultiAtomSDFHasNonzeroDensity) {
  // Two Ar atoms at known positions — SDF should have non-zero density
  correlation::core::Cell cell;
  cell.setLatticeParameters({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {2.0, 5.0, 5.0});
  cell.addAtom("Ar", {8.0, 5.0, 5.0});

  DistributionFunctions df(cell);
  AnalysisSettings settings;
  settings.r_bin_width = 1.0;

  SDFCalculator calc;
  calc.calculateFrame(df, settings);

  const auto &hist = df.getHistogram("SDF");
  EXPECT_EQ(hist.bins.size(), 1000);

  // With 2 atoms, the total SDF should have non-zero values somewhere
  const auto &total = hist.partials.at("Total");
  double sum = 0.0;
  for (double val : total) {
    sum += val;
  }
  EXPECT_GT(sum, 0.0) << "SDF with 2 atoms should have non-zero density";

  // Per-element partial should also have values
  const auto &ar_partial = hist.partials.at("Ar");
  double ar_sum = 0.0;
  for (double val : ar_partial) {
    ar_sum += val;
  }
  EXPECT_GT(ar_sum, 0.0) << "Ar partial SDF should have non-zero density";
}
