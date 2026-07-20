// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/StructureAnalyzer.hpp"
#include "calculators/CNCalculator.hpp"
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
  real_t const neighbor_cutoff = 1.5;
  std::vector<std::vector<real_t>> const bond_cutoffs_sq = {
      {2.25, 2.25}, // C-C, C-H
      {2.25, 2.25}  // H-C, H-H
  };

  // Build StructureAnalyzer to compute neighbor list
  StructureAnalyzer const analyzer(cell, neighbor_cutoff, bond_cutoffs_sq, true);

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

// --- Extreme / Edge-Case Tests ---

TEST(CNCalculatorTests, IsolatedAtomHasZeroCoordination) {
  // Single atom with no neighbors within cutoff
  Cell cell({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0});
  cell.addAtom("Ar", {10.0, 10.0, 10.0});

  real_t const neighbor_cutoff = 3.0;
  std::vector<std::vector<real_t>> const bond_cutoffs_sq = {{9.0}};

  StructureAnalyzer const analyzer(cell, neighbor_cutoff, bond_cutoffs_sq, true);
  auto hist = CNCalculator::calculate(cell, &analyzer);

  // With no neighbors, the CN calculator doesn't create an "Ar-Ar" key at all
  // (it only creates keys for atoms that HAVE neighbors of a given type).
  // But "Any-Any" is always created and should be all zeros for an isolated atom.
  ASSERT_TRUE(hist.partials.find("Any-Any") != hist.partials.end());
  const auto &any_any = hist.partials.at("Any-Any");
  real_t total = 0.0;
  for (real_t const v : any_any) {
    total += v;
  }
  EXPECT_DOUBLE_EQ(total, 0.0) << "Isolated atom should contribute no CN counts";
}

TEST(CNCalculatorTests, HighCoordinationFCC) {
  // FCC unit cell: each atom has 12 nearest neighbors
  Cell cell({1.0, 1.0, 1.0, 90.0, 90.0, 90.0});
  cell.addAtom("Cu", {0.0, 0.0, 0.0});
  cell.addAtom("Cu", {0.5, 0.5, 0.0});
  cell.addAtom("Cu", {0.5, 0.0, 0.5});
  cell.addAtom("Cu", {0.0, 0.5, 0.5});

  // FCC nearest-neighbor distance = a/sqrt(2) ≈ 0.707
  // Use cutoff slightly above that
  real_t const neighbor_cutoff = 0.75;
  std::vector<std::vector<real_t>> const bond_cutoffs_sq = {{0.75 * 0.75}};

  // Need to consider periodic images (ignore_periodic_self_interactions = false)
  StructureAnalyzer const analyzer(cell, neighbor_cutoff, bond_cutoffs_sq, false);
  auto hist = CNCalculator::calculate(cell, &analyzer);

  // In FCC, each atom should have CN=12
  ASSERT_TRUE(hist.partials.find("Cu-Cu") != hist.partials.end());
  const auto &cu_cn = hist.partials.at("Cu-Cu");
  ASSERT_GT(cu_cn.size(), 12);
  // All 4 atoms should have CN=12
  EXPECT_DOUBLE_EQ(cu_cn[12], 4.0);
}

} // namespace correlation::testing
