/**
 * @file ChiralityCalculatorTests.cpp
 * @brief Unit tests for ChiralityCalculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/ChiralityCalculator.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::analysis {
namespace {

class ChiralityCalculatorTests : public ::testing::Test {
protected:
  // Helper to extract the peak bin index for a histogram
  static real_t getPeakValue(const Histogram &hist) {
    const auto &vals = hist.partials.at("Total");
    real_t max_val = -1.0;
    size_t max_bin = 0;
    for (size_t i = 0; i < vals.size(); ++i) {
      if (vals[i] > max_val) {
        max_val = vals[i];
        max_bin = i;
      }
    }
    return hist.bins[max_bin];
  }
};

TEST_F(ChiralityCalculatorTests, AchiralCoplanarMotif) {
  // Central atom at (0,0,0) and 3 coplanar neighbors in XY plane (z = 0)
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0}); // Central atom (index 0)
  cell.addAtom("Ar", {6.0, 5.0, 5.0}); // Neighbor 1: r1 = (1.0, 0, 0), d = 1.0
  cell.addAtom("Ar", {5.0, 6.1, 5.0}); // Neighbor 2: r2 = (0, 1.1, 0), d = 1.1
  cell.addAtom("Ar", {3.8, 3.8, 5.0}); // Neighbor 3: r3 = (-1.2, -1.2, 0), d = 1.697

  StructureAnalyzer const analyzer(cell, 2.5, {{2.5 * 2.5}}, false);
  real_t const chi = correlation::calculators::ChiralityCalculator::computeSingleAtomChirality(0, cell, &analyzer);

  EXPECT_NEAR(chi, 0.0, 1e-7);
}

TEST_F(ChiralityCalculatorTests, ChiralRightHandedMotif) {
  // Central atom at (5.0, 5.0, 5.0) and 3 neighbors forming a right-handed basis
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0}); // Central atom (index 0)
  cell.addAtom("Ar", {6.0, 5.0, 5.0}); // Neighbor 1: r1 = (1.0, 0, 0), d = 1.0
  cell.addAtom("Ar", {5.0, 6.1, 5.0}); // Neighbor 2: r2 = (0, 1.1, 0), d = 1.1
  cell.addAtom("Ar", {5.0, 5.0, 6.2}); // Neighbor 3: r3 = (0, 0, 1.2), d = 1.2

  StructureAnalyzer const analyzer(cell, 2.5, {{2.5 * 2.5}}, false);
  real_t const chi = correlation::calculators::ChiralityCalculator::computeSingleAtomChirality(0, cell, &analyzer);

  EXPECT_NEAR(chi, 1.0, 1e-7);
}

TEST_F(ChiralityCalculatorTests, ChiralLeftHandedMotif) {
  // Central atom at (5.0, 5.0, 5.0) and 3 neighbors forming a left-handed basis
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0}); // Central atom (index 0)
  cell.addAtom("Ar", {6.0, 5.0, 5.0}); // Neighbor 1: r1 = (1.0, 0, 0), d = 1.0
  cell.addAtom("Ar", {5.0, 6.1, 5.0}); // Neighbor 2: r2 = (0, 1.1, 0), d = 1.1
  cell.addAtom("Ar", {5.0, 5.0, 3.8}); // Neighbor 3: r3 = (0, 0, -1.2), d = 1.2

  StructureAnalyzer const analyzer(cell, 2.5, {{2.5 * 2.5}}, false);
  real_t const chi = correlation::calculators::ChiralityCalculator::computeSingleAtomChirality(0, cell, &analyzer);

  EXPECT_NEAR(chi, -1.0, 1e-7);
}

TEST_F(ChiralityCalculatorTests, HistogramDistribution) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});
  cell.addAtom("Ar", {6.0, 5.0, 5.0});
  cell.addAtom("Ar", {5.0, 6.1, 5.0});
  cell.addAtom("Ar", {5.0, 5.0, 6.2});

  StructureAnalyzer const analyzer(cell, 2.5, {{2.5 * 2.5}}, false);
  auto hist = correlation::calculators::ChiralityCalculator::calculate(cell, &analyzer);

  EXPECT_EQ(hist.title, "Chiral Order Parameter Distribution");
  EXPECT_EQ(hist.x_label, "Chiral Order Parameter (\\chi)");
  EXPECT_EQ(hist.file_suffix, "_cop");
  EXPECT_FALSE(hist.partials.empty());

  // There are 4 atoms total.
  // Only atom 0 has 3 neighbors (chirality = 1.0).
  // The other 3 atoms have only 1 neighbor (chirality = 0.0).
  // Thus, the probability distribution should have non-zero contributions at 0.0 and 1.0.
  real_t sum = 0.0;
  for (const auto &val : hist.partials.at("Total")) {
    sum += val;
  }
  EXPECT_NEAR(sum, 1.0, 1e-6); // Probability sums to 1
}

} // namespace
} // namespace correlation::analysis
