// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/ClusterCalculator.hpp"
#include "core/Cell.hpp"
#include <gtest/gtest.h>

using namespace correlation::calculators;
using namespace correlation::analysis;

TEST(ClusterCalculatorTests, BasicClustering) {
  correlation::core::Cell cell;
  cell.setLatticeParameters({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});

  // Cluster 1: 3 atoms
  cell.addAtom("Ar", {2.0, 2.0, 2.0});
  cell.addAtom("Ar", {2.0, 3.0, 2.0}); // Distance 1.0 from atom 0
  cell.addAtom("Ar", {3.0, 2.0, 2.0}); // Distance 1.0 from atom 0

  // Cluster 2: 2 atoms
  cell.addAtom("Ar", {7.0, 7.0, 7.0});
  cell.addAtom("Ar", {7.0, 8.0, 7.0}); // Distance 1.0 from atom 3

  // Cluster 3: 1 atom
  cell.addAtom("Ar", {15.0, 15.0, 15.0});

  DistributionFunctions df(cell);
  AnalysisSettings settings;

  double cutoff = 1.5;
  StructureAnalyzer analyzer(cell, cutoff, {{cutoff * cutoff}}, false);
  df.setStructureAnalyzer(&analyzer);

  ClusterCalculator calc;
  calc.calculateFrame(df, settings);

  const auto &hist = df.getHistogram("Cluster Size");
  // Max size is 3, so bins should be [1, 2, 3]
  ASSERT_EQ(hist.bins.size(), 3);
  EXPECT_DOUBLE_EQ(hist.bins[0], 1.0);
  EXPECT_DOUBLE_EQ(hist.bins[1], 2.0);
  EXPECT_DOUBLE_EQ(hist.bins[2], 3.0);

  EXPECT_TRUE(hist.partials.count("Total") > 0);

  const auto &total = hist.partials.at("Total");
  ASSERT_EQ(total.size(), 3);
  EXPECT_DOUBLE_EQ(total[0], 1.0); // One cluster of size 1
  EXPECT_DOUBLE_EQ(total[1], 1.0); // One cluster of size 2
  EXPECT_DOUBLE_EQ(total[2], 1.0); // One cluster of size 3
}

TEST(ClusterCalculatorTests, SingleGiantCluster) {
  correlation::core::Cell cell;
  cell.setLatticeParameters({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});

  // Create a chain of 10 atoms
  for (int i = 0; i < 10; ++i) {
    cell.addAtom("C", {static_cast<double>(i) * 1.0 + 0.5, 5.0, 5.0});
  }

  DistributionFunctions df(cell);
  AnalysisSettings settings;

  double cutoff = 1.5;
  StructureAnalyzer analyzer(cell, cutoff, {{cutoff * cutoff}}, false);
  df.setStructureAnalyzer(&analyzer);

  ClusterCalculator calc;
  calc.calculateFrame(df, settings);

  const auto &hist = df.getHistogram("Cluster Size");
  EXPECT_EQ(hist.bins.size(), 10);
  const auto &total = hist.partials.at("Total");

  for (size_t i = 0; i < 9; ++i) {
    EXPECT_DOUBLE_EQ(total[i], 0.0);
  }
  EXPECT_DOUBLE_EQ(total[9], 1.0); // One cluster of size 10
}

TEST(ClusterCalculatorTests, EmptyCell) {
  correlation::core::Cell cell;
  cell.setLatticeParameters({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  EXPECT_THROW(DistributionFunctions df(cell), std::invalid_argument);
}
