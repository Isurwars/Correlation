// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/StructureAnalyzer.hpp"
#include "calculators/DADCalculator.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>

namespace correlation::analysis {

class _18_DAD_Tests : public ::testing::Test {
protected:
  correlation::core::Cell cell;

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

TEST_F(_18_DAD_Tests, BasicCalculation) {
  // Cutoff must be > 1.0 to find the bonds (dist is 1.0 each)
  double r_cut = 1.5;
  std::vector<std::vector<double>> bond_cutoffs(
      1, std::vector<double>(1, r_cut * r_cut));
  StructureAnalyzer analyzer(cell, r_cut, bond_cutoffs, true);

  double bin_width = 10.0;
  Histogram f_dihedral = correlation::calculators::DADCalculator::calculate(
      cell, &analyzer, bin_width);

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

TEST_F(_18_DAD_Tests, IcosahedronAnglesDAD) {
  correlation::core::Cell cell_iso;
  cell_iso.setLatticeParameters({30.0, 30.0, 30.0, 90.0, 90.0, 90.0});
  cell_iso.addAtom("Si", {10.0, 10.0, 10.0}); // Center
  double phi = (1.0 + std::sqrt(5.0)) / 2.0;
  std::vector<std::vector<double>> verts = {
      {0, 1, phi}, {0, 1, -phi}, {0, -1, phi}, {0, -1, -phi},
      {1, phi, 0}, {1, -phi, 0}, {-1, phi, 0}, {-1, -phi, 0},
      {phi, 0, 1}, {phi, 0, -1}, {-phi, 0, 1}, {-phi, 0, -1}};

  for (const auto &v : verts) {
    cell_iso.addAtom("Si", {10.0 + v[0], 10.0 + v[1], 10.0 + v[2]});
  }

  double r_cut = 2.5;
  std::vector<std::vector<double>> bond_cutoffs(
      1, std::vector<double>(1, r_cut * r_cut));
  StructureAnalyzer analyzer(cell_iso, r_cut, bond_cutoffs, true);

  double bin_width = 1.0;
  Histogram f_dihedral = correlation::calculators::DADCalculator::calculate(
      cell_iso, &analyzer, bin_width);

  ASSERT_FALSE(f_dihedral.partials.empty());
  auto &partial = f_dihedral.partials["Si-Si-Si-Si"];

  // DAD expects multiple angles due to Center-Vertex and Vertex-Vertex chains
  std::vector<double> expected_angles = {0.0,   31.7,  36.0,   63.4,  72.0,
                                         100.0, 108.0, 138.19, 144.0, 180.0};

  for (double target : expected_angles) {
    bool found = false;
    for (size_t i = 0; i < partial.size(); ++i) {
      if (partial[i] > 1e-4) {
        double angle = std::abs(f_dihedral.bins[i]);
        if (std::abs(angle - target) < 2.0) { // Tolerance considering binning
          found = true;
          break;
        }
      }
    }
    EXPECT_TRUE(found) << "Should find DAD peak near " << target << " degrees";
  }
}
} // namespace correlation::analysis
