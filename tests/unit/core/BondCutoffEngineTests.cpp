// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/BondCutoffEngine.hpp"
#include "core/Cell.hpp"
#include "physics/PhysicalData.hpp"

#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::core;

TEST(BondCutoffEngineTests, ComputeFromElementsCalculatesCorrectRanges) {
  const std::vector<Element> elements{
      Element("H"),
      Element("O"),
  };

  const auto cutoffs = BondCutoffEngine::compute(elements);
  ASSERT_EQ(cutoffs.size(), 2);
  ASSERT_EQ(cutoffs[0].size(), 2);
  ASSERT_EQ(cutoffs[1].size(), 2);

  const real_t r_h = correlation::physics::getCovalentRadius("H");
  const real_t r_o = correlation::physics::getCovalentRadius("O");
  const real_t expected_h_o = (r_h + r_o) * static_cast<real_t>(1.3);

  EXPECT_NEAR(std::sqrt(cutoffs[0][1].max_sq), expected_h_o, 1e-4);
  EXPECT_NEAR(std::sqrt(cutoffs[1][0].max_sq), expected_h_o, 1e-4);
}

TEST(BondCutoffEngineTests, ComputeFromCellUsesUniqueElements) {
  Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("O", {1.0, 1.0, 1.0});

  const auto cutoffs = BondCutoffEngine::compute(cell);
  EXPECT_EQ(cutoffs.size(), 2);
}

TEST(BondCutoffEngineTests, ValidateThrowsOnNegativeOrInvertedBounds) {
  const correlation::analysis::BondCutoffMatrix negative_bound{
      {correlation::analysis::BondCutoffRange{.min_sq = -1.0, .max_sq = 2.0}},
  };
  EXPECT_THROW(BondCutoffEngine::validate(negative_bound), std::invalid_argument);

  const correlation::analysis::BondCutoffMatrix inverted_bounds{
      {correlation::analysis::BondCutoffRange{.min_sq = 5.0, .max_sq = 2.0}},
  };
  EXPECT_THROW(BondCutoffEngine::validate(inverted_bounds), std::invalid_argument);

  const correlation::analysis::BondCutoffMatrix valid_bounds{
      {correlation::analysis::BondCutoffRange{.min_sq = 1.0, .max_sq = 3.0}},
  };
  EXPECT_NO_THROW(BondCutoffEngine::validate(valid_bounds));
}

TEST(BondCutoffEngineTests, AccessorsHandleInAndOutOfBounds) {
  const correlation::analysis::BondCutoffMatrix matrix{
      {correlation::analysis::BondCutoffRange{.min_sq = 4.0, .max_sq = 9.0}},
  };

  EXPECT_NEAR(BondCutoffEngine::getBondCutoffSQ(matrix, 0, 0), 9.0, 1e-5);
  EXPECT_NEAR(BondCutoffEngine::getMinBondCutoffSQ(matrix, 0, 0), 4.0, 1e-5);
  EXPECT_NEAR(BondCutoffEngine::getBondCutoff(matrix, 0, 0), 3.0, 1e-5);
  EXPECT_NEAR(BondCutoffEngine::getMinBondCutoff(matrix, 0, 0), 2.0, 1e-5);

  // Out of bounds returns default 0 range
  EXPECT_NEAR(BondCutoffEngine::getBondCutoffSQ(matrix, 1, 0), 0.0, 1e-5);
  EXPECT_NEAR(BondCutoffEngine::getMinBondCutoffSQ(matrix, 0, 1), 0.0, 1e-5);
}

} // namespace correlation::testing
