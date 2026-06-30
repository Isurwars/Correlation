// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/StructureAnalyzer.hpp"
#include "calculators/DADCalculator.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>
#include <numbers>

namespace correlation::analysis {
namespace {
class DADTests : public ::testing::Test {
public:
  correlation::core::Cell cell_;

protected:
  void SetUp() override {
    // We will place 4 atoms in a sequence A-B-C-D to test dihedral.
    // Let's use 4 carbons for simplicity.
    cell_.addAtom("C", {1.0, 0.0, 0.0});
    cell_.addAtom("C", {0.0, 0.0, 0.0});
    cell_.addAtom("C", {0.0, 1.0, 0.0});
    cell_.addAtom("C", {0.0, 1.0, 1.0});

    // To prevent interactions across periodic boundaries, use a large cell.
    cell_.setLatticeParameters({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  }
};
} // namespace

TEST_F(DADTests, BasicCalculation) {
  // Cutoff must be > 1.0 to find the bonds (dist is 1.0 each)
  double const r_cut = 1.5;
  std::vector<std::vector<double>> const bond_cutoffs(1, std::vector<double>(1, r_cut * r_cut));
  StructureAnalyzer const analyzer(cell_, r_cut, bond_cutoffs, true);

  double const bin_width = 10.0;
  Histogram f_dihedral = correlation::calculators::DADCalculator::calculate(cell_, &analyzer, bin_width);

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

TEST_F(DADTests, IcosahedronAnglesDAD) {
  correlation::core::Cell cell_iso;
  cell_iso.setLatticeParameters({30.0, 30.0, 30.0, 90.0, 90.0, 90.0});
  cell_iso.addAtom("Si", {10.0, 10.0, 10.0}); // Center
  double phi = std::numbers::phi;
  std::vector<std::vector<double>> const vertices = {{0, 1, phi}, {0, 1, -phi}, {0, -1, phi}, {0, -1, -phi},
                                                     {1, phi, 0}, {1, -phi, 0}, {-1, phi, 0}, {-1, -phi, 0},
                                                     {phi, 0, 1}, {phi, 0, -1}, {-phi, 0, 1}, {-phi, 0, -1}};

  for (const auto &vertex : vertices) {
    cell_iso.addAtom("Si", {10.0 + vertex[0], 10.0 + vertex[1], 10.0 + vertex[2]});
  }

  double const r_cut = 2.5;
  std::vector<std::vector<double>> const bond_cutoffs(1, std::vector<double>(1, r_cut * r_cut));
  StructureAnalyzer const analyzer(cell_iso, r_cut, bond_cutoffs, true);

  double const bin_width = 1.0;
  Histogram f_dihedral = correlation::calculators::DADCalculator::calculate(cell_iso, &analyzer, bin_width);

  ASSERT_FALSE(f_dihedral.partials.empty());
  auto &partial = f_dihedral.partials["Si-Si-Si-Si"];

  // DAD expects multiple angles due to Center-Vertex and Vertex-Vertex chains
  std::vector<double> const expected_angles = {0.0, 31.7, 36.0, 63.4, 72.0, 100.0, 108.0, 138.19, 144.0, 180.0};

  for (double const target : expected_angles) {
    bool found = false;
    for (size_t i = 0; i < partial.size(); ++i) {
      if (partial[i] > 1e-4) {
        double const angle = std::abs(f_dihedral.bins[i]);
        if (std::abs(angle - target) < 2.0) { // Tolerance considering binning
          found = true;
          break;
        }
      }
    }
    EXPECT_TRUE(found) << "Should find DAD peak near " << target << " degrees";
  }
}

TEST_F(DADTests, NullNeighborsThrows) {
  EXPECT_THROW({ correlation::calculators::DADCalculator::calculate(cell_, nullptr, 10.0); }, std::logic_error);
}

} // namespace correlation::analysis
