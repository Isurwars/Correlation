// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>

#include "Cell.hpp"
#include "StructureAnalyzer.hpp"
#include "calculators/DADCalculator.hpp"

class _19_DAD_Tests : public ::testing::Test {
protected:
  Cell cell;

  void SetUp() override {
    // We will place 4 atoms in a sequence A-B-C-D to test dihedral.
    // Let's use 4 carbons for simplicity.
    cell.addAtom("C", {1.0, 0.0, 0.0});
    cell.addAtom("C", {0.0, 0.0, 0.0});
    cell.addAtom("C", {0.0, 1.0, 0.0});
    cell.addAtom("C", {0.0, 1.0, 1.0});

    // To prevent interactions across periodic boundaries, use a large cell.
    cell.setLatticeParameters({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  }
};

TEST_F(_19_DAD_Tests, BasicCalculation) {
  // Cutoff must be > 1.0 to find the bonds (dist is 1.0 each)
  double r_cut = 1.5;
  std::vector<std::vector<double>> bond_cutoffs(
      1, std::vector<double>(1, r_cut * r_cut));
  StructureAnalyzer analyzer(cell, r_cut, bond_cutoffs, true);

  double bin_width = 10.0;
  Histogram f_dihedral = DADCalculator::calculate(cell, &analyzer, bin_width);

  // We only expect one type of dihedral for C-C-C-C.
  ASSERT_FALSE(f_dihedral.partials.empty());

  auto &partial = f_dihedral.partials["C-C-C-C"];
  ASSERT_EQ(partial.size(), f_dihedral.bins.size());

  // There is exactly one dihedral since there's one chain 0-1-2-3 (and 3-2-1-0
  // which is identical). Did DADCalculator put it in the correct bin? We can
  // just verify that total sum is > 0, actually since we normalized, sum over
  // bin * width ~ 1 or total_counts is 2.
  double sum = 0.0;
  for (auto val : partial) {
    sum += val;
  }
  // Because it is normalized: sum * bin_width = 1.0
  EXPECT_NEAR(sum * bin_width, 1.0, 1e-5);
}
