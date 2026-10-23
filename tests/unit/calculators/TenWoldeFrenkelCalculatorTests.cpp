// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../../CrystalTestHelper.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/order/SteinhardtCalculator.hpp"
#include "calculators/order/TenWoldeFrenkelCalculator.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::calculators {
namespace {

class TenWoldeFrenkelCalculatorTests : public ::testing::Test {};

TEST_F(TenWoldeFrenkelCalculatorTests, NullNeighborsThrows) {
  const auto cell = correlation::testing::crystals::createFCCCell(1.0, "Cu", 1, 1, 1);
  EXPECT_THROW(
      { static_cast<void>(TenWoldeFrenkelCalculator::calculatePerAtom(cell, nullptr)); },
      std::logic_error);
}

TEST_F(TenWoldeFrenkelCalculatorTests, EmptyCell) {
  const core::Cell empty_cell;
  const analysis::StructureAnalyzer analyzer(empty_cell, 1.0, {}, false);
  const auto per_atom = TenWoldeFrenkelCalculator::calculatePerAtom(empty_cell, &analyzer);
  EXPECT_TRUE(per_atom.empty());

  const auto hist = TenWoldeFrenkelCalculator::calculate(empty_cell, &analyzer);
  EXPECT_TRUE(hist.bins.empty());
}

TEST_F(TenWoldeFrenkelCalculatorTests, PerfectFCC) {
  // FCC Cu cell, 4 atoms per unit cell. 1x1x1 supercell has coordination 12.
  const auto cell = correlation::testing::crystals::createFCCCell(1.0, "Cu", 1, 1, 1);
  const analysis::StructureAnalyzer analyzer(cell, 0.8, {{{.min_sq = 0.36, .max_sq = 0.8 * 0.8}}},
                                             false);
  ASSERT_EQ(analyzer.neighborGraph().getNeighbors(0).size(), 12UL);

  const auto n_conn = TenWoldeFrenkelCalculator::calculatePerAtom(cell, &analyzer);
  ASSERT_EQ(n_conn.size(), cell.atomCount());

  // In perfect FCC, all 12 bonds are mutually coherent with chi_ij >= 0.7
  for (const auto val : n_conn) {
    EXPECT_FLOAT_EQ(static_cast<float>(val), 12.0F);
  }

  // Verify histogram generation
  const auto hist = TenWoldeFrenkelCalculator::calculate(cell, &analyzer);
  ASSERT_FALSE(hist.bins.empty());
  ASSERT_TRUE(hist.partials.contains("Total"));
  ASSERT_TRUE(hist.partials.contains("Cu"));

  // Bin 12 should hold 100% of the distribution
  EXPECT_FLOAT_EQ(static_cast<float>(hist.partials.at("Total")[12]), 1.0F);
  EXPECT_FLOAT_EQ(static_cast<float>(hist.partials.at("Cu")[12]), 1.0F);

  // Description should state solid-like fraction 1.0000
  EXPECT_NE(hist.description.find("solid-like fraction: 1.0000"), std::string::npos);
}

TEST_F(TenWoldeFrenkelCalculatorTests, PrecomputedSteinhardtEquivalence) {
  const auto cell = correlation::testing::crystals::createFCCCell(1.0, "Cu", 2, 2, 2);
  const analysis::StructureAnalyzer analyzer(cell, 0.8, {{{.min_sq = 0.36, .max_sq = 0.8 * 0.8}}},
                                             false);

  // 1. Without pre-computed Steinhardt (calculates internally)
  const auto n_conn_auto = TenWoldeFrenkelCalculator::calculatePerAtom(cell, &analyzer);

  // 2. With pre-computed Steinhardt keeping q6m
  const SteinhardtCalculator::SteinhardtOptions opts{.keep_q6m = true};
  const auto steinhardt_params = SteinhardtCalculator::calculatePerAtom(cell, &analyzer, opts);
  ASSERT_TRUE(steinhardt_params.q6m_normalized.has_value());

  const TenWoldeFrenkelParams default_params;
  const auto n_conn_precomputed = TenWoldeFrenkelCalculator::calculatePerAtom(
      cell, &analyzer, default_params, &steinhardt_params);

  ASSERT_EQ(n_conn_auto.size(), n_conn_precomputed.size());
  for (std::size_t i = 0; i < n_conn_auto.size(); ++i) {
    EXPECT_FLOAT_EQ(static_cast<float>(n_conn_auto[i]), static_cast<float>(n_conn_precomputed[i]));
  }
}

TEST_F(TenWoldeFrenkelCalculatorTests, ThresholdFiltering) {
  const auto cell = correlation::testing::crystals::createFCCCell(1.0, "Cu", 1, 1, 1);
  const analysis::StructureAnalyzer analyzer(cell, 0.8, {{{.min_sq = 0.36, .max_sq = 0.8 * 0.8}}},
                                             false);

  // Very strict chi_cut = 0.9999 should reject bonds unless perfectly aligned
  const TenWoldeFrenkelParams strict_params{.chi_cut = 0.9999, .nconn_min = 12};
  const auto n_conn_strict =
      TenWoldeFrenkelCalculator::calculatePerAtom(cell, &analyzer, strict_params);
  ASSERT_EQ(n_conn_strict.size(), cell.atomCount());
  for (const auto val : n_conn_strict) {
    EXPECT_LE(val, 12.0);
  }

  // Lenient chi_cut = 0.1 should keep all 12 bonds
  const TenWoldeFrenkelParams lenient_params{.chi_cut = 0.1, .nconn_min = 7};
  const auto n_conn_lenient =
      TenWoldeFrenkelCalculator::calculatePerAtom(cell, &analyzer, lenient_params);
  ASSERT_EQ(n_conn_lenient.size(), cell.atomCount());
  for (const auto val : n_conn_lenient) {
    EXPECT_FLOAT_EQ(static_cast<float>(val), 12.0F);
  }
}

TEST_F(TenWoldeFrenkelCalculatorTests, CalculateFrameIntegration) {
  const auto cell = correlation::testing::crystals::createFCCCell(1.0, "Cu", 1, 1, 1);
  analysis::DistributionFunctions dists(cell, 0.8, {{{.min_sq = 0.36, .max_sq = 0.8 * 0.8}}});
  analysis::AnalysisSettings settings;
  settings.twf_chi_cut = 0.7;
  settings.twf_nconn_min = 7;

  const TenWoldeFrenkelCalculator calc;
  calc.calculateFrame(dists, settings);

  EXPECT_TRUE(dists.getAllHistograms().contains("TWF_Nconn"));
  const auto &hist = dists.getAllHistograms().at("TWF_Nconn");
  EXPECT_FALSE(hist.bins.empty());
  EXPECT_TRUE(hist.partials.contains("Total"));
}

} // namespace
} // namespace correlation::calculators
