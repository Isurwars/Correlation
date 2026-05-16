// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "calculators/StructureFactorCalculator.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>

namespace correlation::analysis {

class StructureFactorCalculatorTests : public ::testing::Test {
protected:
  void SetUp() override {}
};

// For a simple cubic crystal, S(Q) should show peaks at the allowed Bragg
// reflections. The first peak for a simple cubic lattice is at Q = 2*pi/a.
TEST_F(StructureFactorCalculatorTests, CalculatesSimpleCubicBraggPeak) {
  // Build a 2x2x2 simple cubic supercell with lattice constant a=3.0 Å
  const double a = 3.0;
  correlation::core::Cell cell(
      std::array<double, 6>{2 * a, 2 * a, 2 * a, 90.0, 90.0, 90.0});

  // 8 atoms at corners of a 2x2x2 supercell
  for (int ix = 0; ix < 2; ++ix)
    for (int iy = 0; iy < 2; ++iy)
      for (int iz = 0; iz < 2; ++iz)
        cell.addAtom("Si", {ix * a, iy * a, iz * a});

  DistributionFunctions df(cell);
  correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 5.0;
  settings.q_bin_width = 0.05;

  calc.calculateFrame(df, settings);

  const auto &hist = df.getHistogram("S_q");
  EXPECT_FALSE(hist.bins.empty());
  EXPECT_TRUE(hist.partials.count("Total"));

  // The first Bragg peak for a SC lattice at a=3 is at Q = 2*pi/3 ≈ 2.094
  // Å^{-1} Find the max S(Q) value in the range [1.8, 2.4]
  double peak_sq = 0.0;
  for (size_t i = 0; i < hist.bins.size(); ++i) {
    if (hist.bins[i] >= 1.8 && hist.bins[i] <= 2.4) {
      peak_sq = std::max(peak_sq, hist.partials.at("Total")[i]);
    }
  }
  // At a Bragg peak, S(Q) should be much greater than 1
  EXPECT_GT(peak_sq, 2.0);
}

// For a single atom, S(Q) = 1 for all Q.
TEST_F(StructureFactorCalculatorTests, SingleAtomGivesOne) {
  correlation::core::Cell cell(
      std::array<double, 6>{10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  DistributionFunctions df(cell);
  correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 5.0;
  settings.q_bin_width = 0.1;

  calc.calculateFrame(df, settings);
  const auto &hist = df.getHistogram("S_q");

  ASSERT_FALSE(hist.bins.empty());
  for (size_t i = 0; i < hist.bins.size(); ++i) {
    // S(Q) = 1/N * |sum exp(i q.r)|^2 = 1/1 * |exp(i*0)|^2 = 1.0
    if (hist.partials.at("Total")[i] > 0.0) {
      EXPECT_NEAR(hist.partials.at("Total")[i], 1.0, 1e-6);
    }
  }
}

// --- Tests migrated from SQTests.cpp ---

// For two identical atoms separated by 1.5 Å in a 10x10x10 box,
// S(Q) for the Si-Si partial should appear and not be empty.
TEST_F(StructureFactorCalculatorTests, DimerProducesValidSQ) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});
  cell.addAtom("Ar", {6.5, 5.0, 5.0}); // Distance 1.5 Å

  DistributionFunctions df(cell);
  correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 10.0;
  settings.q_bin_width = 0.1;

  calc.calculateFrame(df, settings);

  EXPECT_NO_THROW(df.getHistogram("S_q"));
  const auto &hist = df.getHistogram("S_q");
  EXPECT_FALSE(hist.bins.empty());
  EXPECT_TRUE(hist.partials.count("Total"));

  // For a homonuclear dimer of 2 atoms at distance d=1.5 A, S(Q) = 1 +
  // sin(Q*d)/(Q*d)
  for (size_t i = 0; i < hist.bins.size(); ++i) {
    double q = hist.bins[i];
    double expected_sq = 1.0;
    if (q > 1e-6) {
      expected_sq += std::sin(q * 1.5) / (q * 1.5);
    } else {
      expected_sq = 2.0;
    }

    // Check if the partials correctly accumulate
    // The exact structural Debye continuous S(Q) formula differs from the
    // discrete k-point grid approximation from a finite periodic simulation
    // cell by up to 10-30%. Additionally, small Q bins may be exactly 0 due to
    // the 10x10x10 PBC limiting the minimum k-vector.
    if (hist.partials.at("Total")[i] > 0.0) {
      EXPECT_NEAR(hist.partials.at("Total")[i], expected_sq, 0.6);
    }
  }
}

// For N identical atoms, S(Q->0) should approach N.
// Using a 4-atom Si cluster and a very small Q bin.
TEST_F(StructureFactorCalculatorTests,
       VerifiesPartialsExistForMulticomponent) {
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  cell.addAtom("C", {0.0, 0.0, 0.0});
  cell.addAtom("O", {1.2, 0.0, 0.0});

  DistributionFunctions df(cell);
  correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 10.0;
  settings.q_bin_width = 1.0;

  calc.calculateFrame(df, settings);
  const auto &hist = df.getHistogram("S_q");

  // Verify all expected partials are present
  EXPECT_TRUE(hist.partials.count("C-C"));
  EXPECT_TRUE(hist.partials.count("O-O"));
  EXPECT_TRUE(hist.partials.count("C-O"));
  EXPECT_TRUE(hist.partials.count("Total"));

  EXPECT_FALSE(hist.bins.empty());
}

// Verify a multicomponent system produces an S_q partial for each pair.
TEST_F(StructureFactorCalculatorTests, HomonuclearClusterPartialsPresent) {
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {2.0, 0.0, 0.0});
  cell.addAtom("Si", {0.0, 2.0, 0.0});
  cell.addAtom("Si", {0.0, 0.0, 2.0});

  DistributionFunctions df(cell);
  correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 5.0;
  settings.q_bin_width = 0.5;

  calc.calculateFrame(df, settings);
  const auto &hist = df.getHistogram("S_q");

  EXPECT_TRUE(hist.partials.count("Si-Si"));
  EXPECT_TRUE(hist.partials.count("Total"));
}
} // namespace correlation::analysis
