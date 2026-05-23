// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/CNCalculator.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {

using namespace correlation::calculators;
using namespace correlation::analysis;
using namespace correlation::core;

TEST(CNCalculatorTests, CalculatesCorrectCoordinationNumbers) {
  // Construct a cubic cell
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});

  // Add 2 atoms: C (0,0,0) and H (1.0, 0.0, 0.0)
  cell.addAtom("C", {0.0, 0.0, 0.0});
  cell.addAtom("H", {1.0, 0.0, 0.0});

  // Define cutoffs
  double neighbor_cutoff = 1.5;
  std::vector<std::vector<double>> bond_cutoffs_sq = {
      {2.25, 2.25}, // C-C, C-H
      {2.25, 2.25}  // H-C, H-H
  };

  // Build StructureAnalyzer to compute neighbor list
  StructureAnalyzer analyzer(cell, neighbor_cutoff, bond_cutoffs_sq, true);

  // Act
  auto hist = CNCalculator::calculate(cell, &analyzer);

  // Assert
  // Bin index 1 corresponds to coordination number 1
  ASSERT_TRUE(hist.partials.find("C-H") != hist.partials.end());
  EXPECT_DOUBLE_EQ(hist.partials.at("C-H")[1], 1.0);

  ASSERT_TRUE(hist.partials.find("H-C") != hist.partials.end());
  EXPECT_DOUBLE_EQ(hist.partials.at("H-C")[1], 1.0);

  ASSERT_TRUE(hist.partials.find("C-Any") != hist.partials.end());
  EXPECT_DOUBLE_EQ(hist.partials.at("C-Any")[1], 1.0);

  ASSERT_TRUE(hist.partials.find("Any-Any") != hist.partials.end());
  EXPECT_DOUBLE_EQ(hist.partials.at("Any-Any")[1], 2.0); // 1 C + 1 H = 2 atoms with CN=1
}

} // namespace correlation::testing
