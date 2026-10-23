// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../../CrystalTestHelper.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "calculators/order/TetrahedralCalculator.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::calculators {
namespace {

class TetrahedralCalculatorTests : public ::testing::Test {};

TEST_F(TetrahedralCalculatorTests, EmptyCell) {
  const core::Cell empty_cell;
  const auto per_atom = TetrahedralCalculator::calculatePerAtom(empty_cell);
  EXPECT_TRUE(per_atom.empty());

  const auto hist = TetrahedralCalculator::calculate(empty_cell);
  EXPECT_FALSE(hist.bins.empty());
  EXPECT_TRUE(hist.partials.empty());
}

TEST_F(TetrahedralCalculatorTests, UnderCoordinatedCell) {
  // Cell with 3 atoms: fewer than 4 neighbors exist, expect 0.0 for all atoms
  core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {1.0, 0.0, 0.0});
  cell.addAtom("Si", {0.0, 1.0, 0.0});

  const auto q_tet = TetrahedralCalculator::calculatePerAtom(cell);
  ASSERT_EQ(q_tet.size(), 3UL);
  for (const auto val : q_tet) {
    EXPECT_FLOAT_EQ(static_cast<float>(val), 0.0F);
  }
}

TEST_F(TetrahedralCalculatorTests, DiamondSiliconPerAtomValues) {
  // Diamond Silicon (lat_a = 5.43 Å) has coordination 4 with ideal tetrahedral angles
  // cos(psi) = -1/3 => q_tet == 1.0 for all bulk atoms
  const auto cell = correlation::testing::crystals::createDiamondCell(5.43, "Si", 2, 2, 2);
  ASSERT_GT(cell.atomCount(), 0UL);

  // First coordination shell for Diamond Si (cutoff 3.0 Å encompasses 4 nearest neighbors)
  const analysis::StructureAnalyzer analyzer(cell, 3.0, {{{.min_sq = 1.0, .max_sq = 9.0}}}, false);
  const auto q_tet = TetrahedralCalculator::calculatePerAtom(cell, &analyzer);
  ASSERT_EQ(q_tet.size(), cell.atomCount());

  for (const auto val : q_tet) {
    EXPECT_NEAR(static_cast<double>(val), 1.0, 1e-3);
  }
}

TEST_F(TetrahedralCalculatorTests, DiamondSiliconHistogramDistribution) {
  const auto cell = correlation::testing::crystals::createDiamondCell(5.43, "Si", 2, 2, 2);
  const analysis::StructureAnalyzer analyzer(cell, 3.0, {{{.min_sq = 1.0, .max_sq = 9.0}}}, false);

  const auto hist = TetrahedralCalculator::calculate(cell, &analyzer);
  ASSERT_FALSE(hist.bins.empty());
  ASSERT_TRUE(hist.partials.contains("Total"));
  ASSERT_TRUE(hist.partials.contains("Si"));

  // The last bin (around 1.0) should hold ~100% of the distribution
  const auto &total_dist = hist.partials.at("Total");
  EXPECT_NEAR(static_cast<double>(total_dist.back()), 1.0, 1e-2);
}

TEST_F(TetrahedralCalculatorTests, DiamondSiliconWithoutAnalyzer) {
  // Test fallback to Cell minimum image search when analyzer == nullptr
  const auto cell = correlation::testing::crystals::createDiamondCell(5.43, "Si", 2, 2, 2);
  const auto q_tet = TetrahedralCalculator::calculatePerAtom(cell, nullptr);

  ASSERT_EQ(q_tet.size(), cell.atomCount());
  for (const auto val : q_tet) {
    EXPECT_NEAR(static_cast<double>(val), 1.0, 1e-3);
  }
}

TEST_F(TetrahedralCalculatorTests, FactoryRegistration) {
  const auto *calc_by_name = CalculatorFactory::instance().getCalculator("Tetrahedral");
  ASSERT_NE(calc_by_name, nullptr);
  EXPECT_EQ(calc_by_name->getName(), "Tetrahedral");
  EXPECT_EQ(calc_by_name->getShortName(), "q_tet");

  const auto *calc_by_short = CalculatorFactory::instance().getCalculator("q_tet");
  ASSERT_NE(calc_by_short, nullptr);
  EXPECT_EQ(calc_by_short, calc_by_name);
}

TEST_F(TetrahedralCalculatorTests, CalculateFrameDistribution) {
  const auto cell = correlation::testing::crystals::createDiamondCell(5.43, "Si", 1, 1, 1);
  analysis::DistributionFunctions dists(cell);
  const TetrahedralCalculator calc;
  const analysis::AnalysisSettings settings;

  calc.calculateFrame(dists, settings);

  EXPECT_TRUE(dists.getAllHistograms().contains("Q_tet"));
  const auto &hist = dists.getAllHistograms().at("Q_tet");
  EXPECT_EQ(hist.title, "Chau-Hardwick Tetrahedral Order Parameter");
  EXPECT_FALSE(hist.bins.empty());
}

} // namespace
} // namespace correlation::calculators
