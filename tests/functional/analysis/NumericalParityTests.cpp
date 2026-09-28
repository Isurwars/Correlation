// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

/**
 * @file NumericalParityTests.cpp
 * @brief Numerical accuracy and float/double parity tests for canonical crystallographic systems.
 */

#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/order/SteinhardtCalculator.hpp"
#include "calculators/scattering/StructureFactorCalculator.hpp"
#include "calculators/spatial/CNCalculator.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/Precision.hpp"

#include "../../CrystalTestHelper.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <numbers>
#include <utility>
#include <vector>

namespace correlation::testing {

namespace {

using correlation::analysis::DistributionFunctions;
using correlation::analysis::StructureAnalyzer;
using correlation::calculators::CNCalculator;
using correlation::calculators::SteinhardtCalculator;
using correlation::calculators::StructureFactorCalculator;
using correlation::core::Cell;
using correlation::core::Trajectory;

/// Asserts numerical parity against analytical reference with tolerance tiering.
void assertRelativeParity(double actual, double expected, double max_rel_err_single,
                          double max_rel_err_double) {
  double const allowed_rel_err = is_single_precision ? max_rel_err_single : max_rel_err_double;
  double const abs_err = std::abs(actual - expected);
  double const denom = std::max(std::abs(expected), 1e-12);
  double const rel_err = abs_err / denom;

  EXPECT_LE(rel_err, allowed_rel_err)
      << "Actual: " << actual << ", Expected: " << expected << ", RelError: " << rel_err
      << " (Limit: " << allowed_rel_err << ")";
}

/// Asserts parity on binned histogram peak coordinates accounting for bin width discretization.
void assertBinnedParity(double actual, double expected, double bin_width, double max_rel_err_single,
                        double max_rel_err_double) {
  double const allowed_rel_err = is_single_precision ? max_rel_err_single : max_rel_err_double;
  double const abs_err = std::abs(actual - expected);
  double const allowed_abs_err =
      std::max(allowed_rel_err * std::abs(expected), (bin_width / 2.0) + 1e-5);

  EXPECT_LE(abs_err, allowed_abs_err)
      << "Actual: " << actual << ", Expected: " << expected << ", AbsError: " << abs_err
      << " (Allowed: " << allowed_abs_err << ", BinWidth: " << bin_width << ")";
}

/// Helper to locate peak (r, max_val) within a designated bin range.
[[nodiscard]] std::pair<real_t, real_t> findPeakInRange(const std::vector<real_t> &bins,
                                                        const std::vector<real_t> &values,
                                                        real_t min_x, real_t max_x) {
  real_t peak_x = 0;
  real_t max_val = 0;
  for (std::size_t i = 0; i < bins.size() && i < values.size(); ++i) {
    if (bins[i] >= min_x && bins[i] <= max_x && values[i] > max_val) {
      max_val = values[i];
      peak_x = bins[i];
    }
  }
  return {peak_x, max_val};
}

/// Asserts that a coordination distribution contains a single isolated peak at expected_cn.
void assertSingleCoordinationPeak(const std::vector<real_t> &partial, std::size_t expected_cn,
                                  double expected_count) {
  ASSERT_GT(partial.size(), expected_cn);
  for (std::size_t i = 0; i < partial.size(); ++i) {
    double const expected = (i == expected_cn) ? expected_count : 0.0;
    EXPECT_DOUBLE_EQ(static_cast<double>(partial[i]), expected)
        << "Mismatch at coordination bin " << i;
  }
}

/// Finds the bin center of the last non-zero active bin in a 1D histogram.
[[nodiscard]] real_t findLastActiveBinCenter(const std::vector<real_t> &data, real_t min_val,
                                             real_t bin_step) {
  real_t center = 0;
  for (std::size_t bin = 0; bin < data.size(); ++bin) {
    if (data[bin] > 0) {
      center = min_val + (static_cast<real_t>(bin) + static_cast<real_t>(0.5)) * bin_step;
    }
  }
  return center;
}

class NumericalParityTests : public ::testing::Test {
protected:
  static constexpr double SINGLE_PRECISION_REL_TOL = 1e-4;
  static constexpr double DOUBLE_PRECISION_REL_TOL = 1e-12;
};

} // namespace

TEST_F(NumericalParityTests, DiamondSiliconCoordinationParity) {
  // Canonical Diamond cubic Silicon: a = 5.43 Angstroms, 2x2x2 supercell (64 atoms)
  constexpr auto LAT_A = static_cast<real_t>(5.43);
  Cell cell = crystals::createDiamondCell(LAT_A, "Si", 2, 2, 2);
  ASSERT_EQ(cell.atomCount(), 64);

  Trajectory traj;
  traj.addFrame(cell);
  traj.precomputeBondCutoffs();

  // First neighbor distance in Diamond Si is (sqrt(3)/4) * a ≈ 2.351253 Angstroms
  constexpr auto NEIGHBOR_CUTOFF = static_cast<real_t>(2.60);
  StructureAnalyzer const analyzer(cell, NEIGHBOR_CUTOFF, traj.getBondCutoffsSQ());

  auto hist = CNCalculator::calculate(cell, &analyzer);

  // In ideal diamond cubic Si, every single atom has strictly 4 nearest neighbors
  ASSERT_TRUE(hist.partials.contains("Si-Si"));
  const auto &partial = hist.partials.at("Si-Si");

  assertSingleCoordinationPeak(partial, 4U, 64.0);
}

TEST_F(NumericalParityTests, DiamondSiliconPADTetrahedralAngleParity) {
  // Diamond Silicon tetrahedral angle is acos(-1/3) ≈ 109.47122 degrees
  constexpr auto LAT_A = static_cast<real_t>(5.43);
  Cell cell = crystals::createDiamondCell(LAT_A, "Si", 2, 2, 2);

  Trajectory traj;
  traj.addFrame(cell);
  traj.precomputeBondCutoffs();

  constexpr auto NEIGHBOR_CUTOFF = static_cast<real_t>(2.60);
  DistributionFunctions dists(cell, NEIGHBOR_CUTOFF, traj.getBondCutoffsSQ());

  // Calculate PAD with high resolution (0.05 degree binning)
  dists.calculatePAD(0.05);

  const auto &hist = dists.getHistogram("PAD");
  ASSERT_TRUE(hist.partials.contains("Si-Si-Si"));
  const auto &partial = hist.partials.at("Si-Si-Si");

  auto [peak_angle, peak_intensity] =
      findPeakInRange(hist.bins, partial, static_cast<real_t>(105.0), static_cast<real_t>(115.0));

  ASSERT_GT(peak_intensity, 0.0);

  // Theoretical tetrahedral angle in degrees
  constexpr double ANALYTICAL_TETRAHEDRAL_DEG = 109.471220634;
  assertBinnedParity(static_cast<double>(peak_angle), ANALYTICAL_TETRAHEDRAL_DEG, 0.05,
                     SINGLE_PRECISION_REL_TOL, DOUBLE_PRECISION_REL_TOL);
}

TEST_F(NumericalParityTests, BccIronRDFPeakPositionsParity) {
  // BCC Iron: a = 2.8665 Angstroms, 3x3x3 supercell (54 atoms)
  constexpr auto LAT_A = static_cast<real_t>(2.8665);
  Cell cell = crystals::createBCCCell(LAT_A, "Fe", 3, 3, 3);
  ASSERT_EQ(cell.atomCount(), 54);

  Trajectory traj;
  traj.addFrame(cell);
  traj.precomputeBondCutoffs();

  constexpr auto MAX_R = static_cast<real_t>(4.0);
  DistributionFunctions dists(cell, MAX_R, traj.getBondCutoffsSQ());

  // Fine RDF bin width: 0.005 Angstroms
  dists.calculateRDF({.r_max = MAX_R, .r_bin_width = static_cast<real_t>(0.005)});

  const auto &hist = dists.getHistogram("g_r");
  ASSERT_TRUE(hist.partials.contains("Fe-Fe"));
  const auto &partial = hist.partials.at("Fe-Fe");

  // 1st neighbor shell: (sqrt(3)/2) * a ≈ 2.482463 Angstroms
  constexpr double ANALYTICAL_D1 = (std::numbers::sqrt3 / 2.0) * 2.8665;
  auto [peak_r1, max1] =
      findPeakInRange(hist.bins, partial, static_cast<real_t>(2.40), static_cast<real_t>(2.55));
  ASSERT_GT(max1, 0.0);
  assertBinnedParity(static_cast<double>(peak_r1), ANALYTICAL_D1, 0.005, SINGLE_PRECISION_REL_TOL,
                     DOUBLE_PRECISION_REL_TOL);

  // 2nd neighbor shell: a = 2.8665 Angstroms
  constexpr double ANALYTICAL_D2 = 2.8665;
  auto [peak_r2, max2] =
      findPeakInRange(hist.bins, partial, static_cast<real_t>(2.80), static_cast<real_t>(2.95));
  ASSERT_GT(max2, 0.0);
  assertBinnedParity(static_cast<double>(peak_r2), ANALYTICAL_D2, 0.005, SINGLE_PRECISION_REL_TOL,
                     DOUBLE_PRECISION_REL_TOL);
}

TEST_F(NumericalParityTests, BccIronSteinhardtOrderParametersParity) {
  // BCC Iron with 1st + 2nd shell cutoff (3.1 Angstroms covers 8 + 6 = 14 neighbors)
  constexpr auto LAT_A = static_cast<real_t>(2.8665);
  Cell cell = crystals::createBCCCell(LAT_A, "Fe", 2, 2, 2);

  Trajectory traj;
  traj.addFrame(cell);
  traj.precomputeBondCutoffs();

  constexpr auto NEIGHBOR_CUTOFF = static_cast<real_t>(3.1);
  StructureAnalyzer const analyzer(cell, NEIGHBOR_CUTOFF, traj.getBondCutoffsSQ());

  auto hists = SteinhardtCalculator::calculate(cell, &analyzer);

  ASSERT_TRUE(hists.contains("Q4"));
  ASSERT_TRUE(hists.contains("Q6"));
  ASSERT_TRUE(hists.contains("W6_hat"));

  // Check peak bin for Q4, Q6, W6_hat
  const auto &q4_data = hists.at("Q4").partials.at("Total");
  const auto &q6_data = hists.at("Q6").partials.at("Total");
  const auto &w6_data = hists.at("W6_hat").partials.at("Total");

  constexpr std::size_t Q_BINS = 100;
  constexpr auto D_Q = static_cast<real_t>(1.0 / static_cast<double>(Q_BINS));
  constexpr std::size_t W_BINS = 100;
  constexpr auto W_MIN = static_cast<real_t>(-0.2);
  constexpr auto W_MAX = static_cast<real_t>(0.2);
  constexpr auto D_W = (W_MAX - W_MIN) / static_cast<real_t>(W_BINS);

  auto const q4_val = findLastActiveBinCenter(q4_data, 0, D_Q);
  auto const q6_val = findLastActiveBinCenter(q6_data, 0, D_Q);
  auto const w6_val = findLastActiveBinCenter(w6_data, W_MIN, D_W);

  // Analytical values for 14-neighbor BCC shell:
  // Q4 ≈ 0.036369, Q6 ≈ 0.510688, W6_hat ≈ 0.013161
  constexpr double EXPECTED_Q4 = 0.036369;
  constexpr double EXPECTED_Q6 = 0.510688;
  constexpr double EXPECTED_W6 = 0.013161;

  // Histogram binning discretization resolution is ~0.01
  EXPECT_NEAR(static_cast<double>(q4_val), EXPECTED_Q4, 0.015);
  EXPECT_NEAR(static_cast<double>(q6_val), EXPECTED_Q6, 0.015);
  EXPECT_NEAR(static_cast<double>(w6_val), EXPECTED_W6, 0.015);
}

TEST_F(NumericalParityTests, BccIronStructureFactorBraggPeakParity) {
  // BCC Iron S(Q) first Bragg peak (110)
  // Q(110) = (2 * pi / a) * sqrt(2) ≈ 3.099496 Angstrom^-1
  constexpr auto LAT_A = static_cast<real_t>(2.8665);
  Cell cell = crystals::createBCCCell(LAT_A, "Fe", 3, 3, 3);

  DistributionFunctions dists(cell);
  StructureFactorCalculator calc;
  analysis::AnalysisSettings settings;
  settings.q_max = static_cast<real_t>(5.0);
  settings.q_bin_width = static_cast<real_t>(0.02);

  calc.calculateFrame(dists, settings);

  const auto &hist = dists.getHistogram("S_q");
  ASSERT_FALSE(hist.bins.empty());
  ASSERT_TRUE(hist.partials.contains("Total"));
  const auto &total_sq = hist.partials.at("Total");

  auto [peak_q, max_sq] =
      findPeakInRange(hist.bins, total_sq, static_cast<real_t>(2.90), static_cast<real_t>(3.30));

  ASSERT_GT(max_sq, 0.0);

  constexpr double ANALYTICAL_Q110 = (2.0 * std::numbers::pi / 2.8665) * std::numbers::sqrt2;
  assertBinnedParity(static_cast<double>(peak_q), ANALYTICAL_Q110, 0.02, SINGLE_PRECISION_REL_TOL,
                     DOUBLE_PRECISION_REL_TOL);
}

} // namespace correlation::testing
