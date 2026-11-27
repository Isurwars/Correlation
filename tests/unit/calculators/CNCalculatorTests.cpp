// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/StructureAnalyzer.hpp"
#include "calculators/spatial/CNCalculator.hpp"
#include "core/Cell.hpp"

#include "../../CrystalTestHelper.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {

using namespace correlation::calculators;
using namespace correlation::analysis;
using namespace correlation::core;
using correlation::real_t;

TEST(CNCalculatorTests, CalculatesCorrectCoordinationNumbers) {
  // Construct a cubic cell
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});

  // Add 2 atoms: C (0,0,0) and H (1.0, 0.0, 0.0)
  cell.addAtom("C", {0.0, 0.0, 0.0});
  cell.addAtom("H", {1.0, 0.0, 0.0});

  // Define cutoffs
  real_t const neighbor_cutoff = 1.5;
  BondCutoffMatrix const bond_cutoffs = {
      {
          {.min_sq = static_cast<real_t>(0.36), .max_sq = 2.25},
          {.min_sq = static_cast<real_t>(0.36), .max_sq = 2.25},
      }, // C-C, C-H
      {
          {.min_sq = static_cast<real_t>(0.36), .max_sq = 2.25},
          {.min_sq = static_cast<real_t>(0.36), .max_sq = 2.25},
      }, // H-C, H-H
  };

  // Build StructureAnalyzer to compute neighbor list
  StructureAnalyzer const analyzer(cell, neighbor_cutoff, bond_cutoffs, true);

  // Act
  auto hist = CNCalculator::calculate(cell, &analyzer);

  // Assert
  // Bin index 1 corresponds to coordination number 1
  ASSERT_TRUE(hist.partials.contains("C-H"));
  EXPECT_DOUBLE_EQ(hist.partials.at("C-H")[1], 1.0);

  ASSERT_TRUE(hist.partials.contains("H-C"));
  EXPECT_DOUBLE_EQ(hist.partials.at("H-C")[1], 1.0);

  ASSERT_TRUE(hist.partials.contains("C-Any"));
  EXPECT_DOUBLE_EQ(hist.partials.at("C-Any")[1], 1.0);

  ASSERT_TRUE(hist.partials.contains("Any-Any"));
  EXPECT_DOUBLE_EQ(hist.partials.at("Any-Any")[1], 2.0); // 1 C + 1 H = 2 atoms with CN=1
}

// --- Extreme / Edge-Case Tests ---

TEST(CNCalculatorTests, IsolatedAtomHasZeroCoordination) {
  // Single atom with no neighbors within cutoff
  Cell cell({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0});
  cell.addAtom("Ar", {10.0, 10.0, 10.0});

  real_t const neighbor_cutoff = 3.0;
  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = static_cast<real_t>(0.36), .max_sq = 9.0}}};

  StructureAnalyzer const analyzer(cell, neighbor_cutoff, bond_cutoffs, true);
  const auto hist = CNCalculator::calculate(cell, &analyzer);

  // An isolated atom has 0 neighbors, so it is counted in bin 0 of Ar-Ar, Ar-Any, Any-Any, Total.
  ASSERT_TRUE(hist.partials.contains("Ar-Ar"));
  EXPECT_DOUBLE_EQ(hist.partials.at("Ar-Ar")[0], 1.0);
  ASSERT_TRUE(hist.partials.contains("Ar-Any"));
  EXPECT_DOUBLE_EQ(hist.partials.at("Ar-Any")[0], 1.0);
  ASSERT_TRUE(hist.partials.contains("Any-Any"));
  EXPECT_DOUBLE_EQ(hist.partials.at("Any-Any")[0], 1.0);
  ASSERT_TRUE(hist.partials.contains("Total"));
  EXPECT_DOUBLE_EQ(hist.partials.at("Total")[0], 1.0);
}

TEST(CNCalculatorTests, MultiSpeciesAnyAndTotalAggregatePerAtom) {
  Cell cell({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0});
  // Central Si atom bonded to 2 Si neighbors and 2 O neighbors -> total CN = 4
  cell.addAtom("Si", {10.0, 10.0, 10.0});
  cell.addAtom("Si", {11.0, 10.0, 10.0});
  cell.addAtom("Si", {9.0, 10.0, 10.0});
  cell.addAtom("O", {10.0, 11.0, 10.0});
  cell.addAtom("O", {10.0, 9.0, 10.0});

  real_t const neighbor_cutoff = static_cast<real_t>(1.2);
  BondCutoffMatrix const bond_cutoffs = {
      {
          {.min_sq = static_cast<real_t>(0.25), .max_sq = static_cast<real_t>(1.44)},
          {.min_sq = static_cast<real_t>(0.25), .max_sq = static_cast<real_t>(1.44)},
      },
      {
          {.min_sq = static_cast<real_t>(0.25), .max_sq = static_cast<real_t>(1.44)},
          {.min_sq = static_cast<real_t>(0.25), .max_sq = static_cast<real_t>(1.44)},
      },
  };

  StructureAnalyzer const analyzer(cell, neighbor_cutoff, bond_cutoffs, true);
  const auto hist = CNCalculator::calculate(cell, &analyzer);

  // Central Si has 2 Si + 2 O = 4 total neighbors (not two atoms with CN=2).
  ASSERT_TRUE(hist.partials.contains("Si-Any"));
  ASSERT_GT(hist.partials.at("Si-Any").size(), 4U);
  EXPECT_DOUBLE_EQ(hist.partials.at("Si-Any")[4], 1.0);
  EXPECT_DOUBLE_EQ(hist.partials.at("Si-Any")[2], 0.0);
  EXPECT_DOUBLE_EQ(hist.partials.at("Si-Any")[1], 2.0); // The 2 outer Si atoms each have CN=1

  // Total counts must sum strictly to atomCount() (5 atoms) without double-counting.
  real_t total_atoms_counted = 0.0;
  for (real_t const val : hist.partials.at("Total")) {
    total_atoms_counted += val;
  }
  EXPECT_DOUBLE_EQ(total_atoms_counted, 5.0);
}

TEST(CNCalculatorTests, HighCoordinationFCC) {
  // FCC unit cell: each atom has 12 nearest neighbors
  const auto cell = correlation::testing::crystals::createFCCCell(1.0, "Cu", 1, 1, 1);

  // FCC nearest-neighbor distance = a/sqrt(2) ≈ 0.707
  // Use cutoff slightly above that
  real_t const neighbor_cutoff = 0.75;
  BondCutoffMatrix const bond_cutoffs = {
      {{.min_sq = static_cast<real_t>(0.36), .max_sq = 0.75 * 0.75}},
  };

  // Need to consider periodic images (ignore_periodic_self_interactions = false)
  StructureAnalyzer const analyzer(cell, neighbor_cutoff, bond_cutoffs, false);
  const auto hist = CNCalculator::calculate(cell, &analyzer);

  // In FCC, each atom should have CN=12
  ASSERT_TRUE(hist.partials.contains("Cu-Cu"));
  const auto &cu_cn = hist.partials.at("Cu-Cu");
  ASSERT_GT(cu_cn.size(), 12);
  // All 4 atoms should have CN=12
  EXPECT_DOUBLE_EQ(cu_cn[12], 4.0);
}

} // namespace correlation::testing
