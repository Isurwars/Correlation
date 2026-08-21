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
  real_t const r_cut = 1.5;
  BondCutoffMatrix const bond_cutoffs(1, std::vector<BondCutoffRange>(1, BondCutoffRange{0.36, r_cut * r_cut}));
  StructureAnalyzer const analyzer(cell_, r_cut, bond_cutoffs, true);

  real_t const bin_width = 10.0;
  auto results = correlation::calculators::DADCalculator::calculate(cell_, &analyzer, bin_width);
  ASSERT_TRUE(results.contains("DAD"));
  ASSERT_TRUE(results.contains("DAD_raw"));

  Histogram f_dihedral = std::move(results["DAD"]);
  Histogram f_dihedral_raw = std::move(results["DAD_raw"]);

  // We only expect one type of dihedral for C-C-C-C.
  ASSERT_FALSE(f_dihedral.partials.empty());

  auto &partial = f_dihedral.partials["C-C-C-C"];
  ASSERT_EQ(partial.size(), f_dihedral.bins.size());

  // There is exactly one dihedral since there's one chain 0-1-2-3 (and 3-2-1-0
  // which is identical). Did DADCalculator put it in the correct bin? We can
  // just verify that total sum is > 0, actually since we normalized, sum over
  // bin * width ~ 1 or total_counts is 2.
  real_t sum = 0.0;
  for (auto val : partial) {
    sum += val;
  }
  // Because it is normalized: sum * bin_width = 1.0
  EXPECT_NEAR(sum * bin_width, 1.0, 1e-5);

  EXPECT_EQ(f_dihedral_raw.y_unit, "counts");
  ASSERT_TRUE(f_dihedral_raw.partials.contains("C-C-C-C"));
  real_t raw_sum = 0.0;
  for (auto val : f_dihedral_raw.partials["C-C-C-C"]) {
    raw_sum += val;
  }
  EXPECT_GT(raw_sum, 0.0);
}

TEST_F(DADTests, IcosahedronAnglesDAD) {
  correlation::core::Cell cell_iso;
  cell_iso.setLatticeParameters({30.0, 30.0, 30.0, 90.0, 90.0, 90.0});
  cell_iso.addAtom("Si", {10.0, 10.0, 10.0}); // Center
  real_t const phi = std::numbers::phi;
  std::vector<std::vector<real_t>> const vertices = {{0, 1, phi}, {0, 1, -phi}, {0, -1, phi}, {0, -1, -phi},
                                                     {1, phi, 0}, {1, -phi, 0}, {-1, phi, 0}, {-1, -phi, 0},
                                                     {phi, 0, 1}, {phi, 0, -1}, {-phi, 0, 1}, {-phi, 0, -1}};

  for (const auto &vertex : vertices) {
    cell_iso.addAtom("Si", correlation::math::Vector3<real_t>(static_cast<real_t>(10.0) + vertex[0],
                                                              static_cast<real_t>(10.0) + vertex[1],
                                                              static_cast<real_t>(10.0) + vertex[2]));
  }

  real_t const r_cut = 2.5;
  BondCutoffMatrix const bond_cutoffs(1, std::vector<BondCutoffRange>(1, BondCutoffRange{0.36, r_cut * r_cut}));
  StructureAnalyzer const analyzer(cell_iso, r_cut, bond_cutoffs, true);

  real_t const bin_width = 1.0;
  auto results = correlation::calculators::DADCalculator::calculate(cell_iso, &analyzer, bin_width);
  ASSERT_TRUE(results.contains("DAD"));
  Histogram f_dihedral = std::move(results["DAD"]);

  ASSERT_FALSE(f_dihedral.partials.empty());
  auto &partial = f_dihedral.partials["Si-Si-Si-Si"];

  // DAD expects multiple angles due to Center-Vertex and Vertex-Vertex chains
  std::vector<real_t> const expected_angles = {
      static_cast<real_t>(0.0),   static_cast<real_t>(31.7),  static_cast<real_t>(36.0),  static_cast<real_t>(63.4),
      static_cast<real_t>(72.0),  static_cast<real_t>(100.0), static_cast<real_t>(108.0), static_cast<real_t>(138.19),
      static_cast<real_t>(144.0), static_cast<real_t>(180.0)};

  for (real_t const target : expected_angles) {
    bool found = false;
    for (size_t i = 0; i < partial.size(); ++i) {
      if (partial[i] > static_cast<real_t>(1e-4)) {
        real_t const angle = std::abs(f_dihedral.bins[i]);
        if (std::abs(angle - target) < static_cast<real_t>(2.0)) { // Tolerance considering binning
          found = true;
          break;
        }
      }
    }
    EXPECT_TRUE(found) << "Target dihedral angle: " << target << " was not found!";
  }
}

TEST_F(DADTests, NullNeighborsThrows) {
  EXPECT_THROW({ correlation::calculators::DADCalculator::calculate(cell_, nullptr, 10.0); }, std::logic_error);
}

} // namespace correlation::analysis
