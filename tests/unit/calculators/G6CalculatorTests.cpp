/**
 * @file G6CalculatorTests.cpp
 * @brief Unit tests for G6Calculator and CorrelationLengthFit.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "../../CrystalTestHelper.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/spatial/G6Calculator.hpp"
#include "core/Cell.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

namespace correlation::testing {

using namespace correlation::calculators;
using namespace correlation::core;
using namespace correlation::analysis;

TEST(G6CalculatorTests, FitsExponentialDecayAccurately) {
  constexpr size_t num_bins = 10;
  std::vector<real_t> r_bins(num_bins);
  std::vector<real_t> g6_values(num_bins);

  constexpr real_t target_xi = 3.5;
  for (size_t pt_idx = 0; pt_idx < num_bins; ++pt_idx) {
    r_bins[pt_idx] = static_cast<real_t>(pt_idx + 1);
    g6_values[pt_idx] = static_cast<real_t>(std::exp(-r_bins[pt_idx] / target_xi));
  }

  const auto fit = G6Calculator::fitCorrelationLength(r_bins, g6_values);

  EXPECT_FALSE(fit.is_plateau);
  EXPECT_NEAR(fit.xi_6, target_xi, 0.05);
  EXPECT_GT(fit.r_squared, 0.99);
}

TEST(G6CalculatorTests, DetectsCrystallinePlateauForConstantG6) {
  constexpr size_t num_bins = 8;
  std::vector<real_t> r_bins(num_bins);
  std::vector<real_t> g6_values(num_bins, 0.85);

  for (size_t pt_idx = 0; pt_idx < num_bins; ++pt_idx) {
    r_bins[pt_idx] = static_cast<real_t>(pt_idx + 1) * static_cast<real_t>(1.5);
  }

  const auto fit = G6Calculator::fitCorrelationLength(r_bins, g6_values);

  EXPECT_TRUE(fit.is_plateau);
  EXPECT_TRUE(std::isinf(fit.xi_6));
}

TEST(G6CalculatorTests, HandlesInsufficientPointsGracefully) {
  std::vector<real_t> r_bins = {1.0};
  std::vector<real_t> g6_values = {0.5};

  const auto fit = G6Calculator::fitCorrelationLength(r_bins, g6_values);

  EXPECT_DOUBLE_EQ(fit.xi_6, 0.0);
  EXPECT_DOUBLE_EQ(fit.r_squared, 0.0);
  EXPECT_FALSE(fit.is_plateau);
}

TEST(G6CalculatorTests, ThrowsOnInvalidInputs) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Cu", {0.0, 0.0, 0.0});
  cell.addAtom("Cu", {2.0, 0.0, 0.0});

  // Null neighbors
  EXPECT_THROW((void)G6Calculator::calculate(cell, nullptr, 5.0, 0.1), std::logic_error);

  const StructureAnalyzer analyzer(cell, 3.0, {{{.min_sq = 0.0, .max_sq = 9.0}}}, true);

  // Invalid r_max
  EXPECT_THROW((void)G6Calculator::calculate(cell, &analyzer, 0.0, 0.1), std::invalid_argument);
  EXPECT_THROW((void)G6Calculator::calculate(cell, &analyzer, -1.0, 0.1), std::invalid_argument);

  // Invalid r_bin_width
  EXPECT_THROW((void)G6Calculator::calculate(cell, &analyzer, 5.0, 0.0), std::invalid_argument);
  EXPECT_THROW((void)G6Calculator::calculate(cell, &analyzer, 5.0, -0.1), std::invalid_argument);
}

TEST(G6CalculatorTests, EvaluatesFCCLatticeCorrelation) {
  // Build a 2x2x2 FCC Copper supercell (32 atoms)
  const Cell cell = crystals::createFCCCell(3.615, "Cu", 2, 2, 2);
  ASSERT_EQ(cell.atomCount(), 32);

  // Neighbor cutoff capturing first coordination shell (~2.556 A)
  const StructureAnalyzer analyzer(cell, 3.0, {{{.min_sq = 0.0, .max_sq = 9.0}}}, true);

  const auto hist = G6Calculator::calculate(cell, &analyzer, 5.0, 0.1);

  EXPECT_EQ(hist.title, "g₆(r) — Bond-Orientational Correlation Function");
  EXPECT_EQ(hist.x_label, "r");
  EXPECT_EQ(hist.y_label, "g₆(r)");
  EXPECT_EQ(hist.file_suffix, "g6_r");

  ASSERT_TRUE(hist.partials.contains("Total"));
  const auto &g6_curve = hist.partials.at("Total");
  ASSERT_EQ(g6_curve.size(), hist.bins.size());

  // Find the first neighbor shell bin (~2.55 A)
  const auto bin_idx = static_cast<size_t>(2.55 / 0.1);
  ASSERT_LT(bin_idx, g6_curve.size());

  // In an ideal FCC crystal, bond coherence between nearest neighbors is close to 1.0
  EXPECT_GT(g6_curve[bin_idx], 0.7);

  // FCC crystal should reflect crystalline plateau
  EXPECT_NE(hist.description.find("plateau"), std::string::npos);
}

TEST(G6CalculatorTests, CalculateFramePopulatesDistributionFunctions) {
  const Cell cell = crystals::createFCCCell(3.615, "Cu", 2, 2, 2);
  DistributionFunctions dists(cell, 3.0, {{{.min_sq = 0.0, .max_sq = 9.0}}});

  AnalysisSettings settings;
  settings.r_max = 5.0;
  settings.r_bin_width = 0.1;

  const G6Calculator calc;
  calc.calculateFrame(dists, settings);

  EXPECT_TRUE(dists.getAllHistograms().contains("G6_r"));
  const auto &hist = dists.getHistogram("G6_r");
  EXPECT_FALSE(hist.bins.empty());
  EXPECT_TRUE(hist.partials.contains("Total"));
}

} // namespace correlation::testing
