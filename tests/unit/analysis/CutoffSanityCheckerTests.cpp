/**
 * @file CutoffSanityCheckerTests.cpp
 * @brief Unit tests for CutoffSanityChecker validation engine and repair actions.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/CutoffSanityChecker.hpp"
#include "core/Atom.hpp"
#include "core/Cell.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <vector>

using namespace correlation::analysis;
using namespace correlation::core;

namespace {

class CutoffSanityCheckerTests : public ::testing::Test {
protected:
  Cell cell_{{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}}; // L_min = 10.0, L_min/2 = 5.0
  std::vector<Element> elements_{
      Element{.symbol = "Si", .id = ElementID{0}},
      Element{.symbol = "O", .id = ElementID{1}},
  };
};

TEST_F(CutoffSanityCheckerTests, ReturnsWarningOnEmptyElements) {
  const BondCutoffMatrix empty_matrix;
  const auto report = CutoffSanityChecker::checkCutoffs(empty_matrix, {}, &cell_);

  EXPECT_EQ(report.severity, SanitySeverity::Warning);
  EXPECT_FALSE(report.has_fix);
  EXPECT_EQ(report.title, "No Elements");
}

TEST_F(CutoffSanityCheckerTests, ReturnsErrorOnDimensionMismatch) {
  const BondCutoffMatrix matrix(1); // Mismatch with elements_.size() == 2
  const auto report = CutoffSanityChecker::checkCutoffs(matrix, elements_, &cell_);

  EXPECT_EQ(report.severity, SanitySeverity::Error);
  EXPECT_TRUE(report.has_fix);
  EXPECT_EQ(report.fix_action, SanityFixAction::ResetToCovalent);
}

TEST_F(CutoffSanityCheckerTests, DetectsOptimalCutoffs) {
  BondCutoffMatrix matrix(2);
  // Si-Si (0, 0): [1.5, 3.0]
  matrix[0][0].min_sq = 1.5 * 1.5;
  matrix[0][0].max_sq = 3.0 * 3.0;
  // Si-O (0, 1): [1.2, 2.2]
  matrix[0][1].min_sq = 1.2 * 1.2;
  matrix[0][1].max_sq = 2.2 * 2.2;
  // O-O (1, 1): [1.8, 3.2]
  matrix[1][1].min_sq = 1.8 * 1.8;
  matrix[1][1].max_sq = 3.2 * 3.2;

  const auto report = CutoffSanityChecker::checkCutoffs(matrix, elements_, &cell_);
  EXPECT_EQ(report.severity, SanitySeverity::Optimal);
  EXPECT_FALSE(report.has_fix);
  EXPECT_EQ(report.title, "Optimal Cutoffs");
}

TEST_F(CutoffSanityCheckerTests, DetectsRangeInversion) {
  BondCutoffMatrix matrix(2);
  // Inversion for Si-O: min (2.5) > max (2.0)
  matrix[0][0].min_sq = 1.5 * 1.5;
  matrix[0][0].max_sq = 3.0 * 3.0;
  matrix[0][1].min_sq = 2.5 * 2.5;
  matrix[0][1].max_sq = 2.0 * 2.0;
  matrix[1][1].min_sq = 1.8 * 1.8;
  matrix[1][1].max_sq = 3.2 * 3.2;

  const auto report = CutoffSanityChecker::checkCutoffs(matrix, elements_, &cell_);
  EXPECT_EQ(report.severity, SanitySeverity::Error);
  EXPECT_TRUE(report.has_fix);
  EXPECT_EQ(report.fix_action, SanityFixAction::ResetToCovalent);
  EXPECT_EQ(report.title, "Inverted Cutoff Bounds");
  EXPECT_EQ(report.problematic_pair, "Si-O");
}

TEST_F(CutoffSanityCheckerTests, DetectsPBCViolation) {
  BondCutoffMatrix matrix(2);
  // PBC violation for Si-O: max (6.0) > L_min / 2 (5.0)
  matrix[0][0].min_sq = 1.5 * 1.5;
  matrix[0][0].max_sq = 3.0 * 3.0;
  matrix[0][1].min_sq = 1.2 * 1.2;
  matrix[0][1].max_sq = 6.0 * 6.0;
  matrix[1][1].min_sq = 1.8 * 1.8;
  matrix[1][1].max_sq = 3.2 * 3.2;

  const auto report = CutoffSanityChecker::checkCutoffs(matrix, elements_, &cell_);
  EXPECT_EQ(report.severity, SanitySeverity::Error);
  EXPECT_TRUE(report.has_fix);
  EXPECT_EQ(report.fix_action, SanityFixAction::ClampToBoxHalf);
  EXPECT_EQ(report.title, "PBC Violation (r_max > L/2)");
  EXPECT_EQ(report.problematic_pair, "Si-O");
  EXPECT_NEAR(report.suggested_value, 4.8, 1e-3);
}

TEST_F(CutoffSanityCheckerTests, DetectsMissingCutoff) {
  BondCutoffMatrix matrix(2);
  // Missing cutoff for O-O: max_sq == 0
  matrix[0][0].min_sq = 1.5 * 1.5;
  matrix[0][0].max_sq = 3.0 * 3.0;
  matrix[0][1].min_sq = 1.2 * 1.2;
  matrix[0][1].max_sq = 2.2 * 2.2;
  matrix[1][1].min_sq = 0.0;
  matrix[1][1].max_sq = 0.0;

  const auto report = CutoffSanityChecker::checkCutoffs(matrix, elements_, &cell_);
  EXPECT_EQ(report.severity, SanitySeverity::Warning);
  EXPECT_TRUE(report.has_fix);
  EXPECT_EQ(report.fix_action, SanityFixAction::AutoDetectFirstMin);
  EXPECT_EQ(report.title, "Unset Cutoff Distance");
}

TEST_F(CutoffSanityCheckerTests, DetectsPlausibilityWarnings) {
  BondCutoffMatrix matrix(2);
  // Unphysically small max cutoff (0.6 Å < MIN_PLAUSIBLE_CUTOFF 0.8 Å)
  matrix[0][0].min_sq = 0.1 * 0.1;
  matrix[0][0].max_sq = 0.6 * 0.6;
  matrix[0][1].min_sq = 1.2 * 1.2;
  matrix[0][1].max_sq = 2.2 * 2.2;
  matrix[1][1].min_sq = 1.8 * 1.8;
  matrix[1][1].max_sq = 3.2 * 3.2;

  const auto report = CutoffSanityChecker::checkCutoffs(matrix, elements_, &cell_);
  EXPECT_EQ(report.severity, SanitySeverity::Warning);
  EXPECT_TRUE(report.has_fix);
  EXPECT_EQ(report.fix_action, SanityFixAction::ResetToCovalent);
  EXPECT_EQ(report.title, "Unusually Short Cutoff");
}

TEST_F(CutoffSanityCheckerTests, AppliesFixClampToBoxHalf) {
  BondCutoffMatrix matrix(2);
  matrix[0][0].min_sq = 1.5 * 1.5;
  matrix[0][0].max_sq = 3.0 * 3.0;
  matrix[0][1].min_sq = 1.2 * 1.2;
  matrix[0][1].max_sq = 6.0 * 6.0; // Violates PBC (5.0 Å)
  matrix[1][1].min_sq = 1.8 * 1.8;
  matrix[1][1].max_sq = 3.2 * 3.2;

  const bool applied =
      CutoffSanityChecker::applyFix(matrix, elements_, &cell_, SanityFixAction::ClampToBoxHalf);
  EXPECT_TRUE(applied);

  const real_t new_max_d = std::sqrt(matrix[0][1].max_sq);
  EXPECT_LE(new_max_d, 4.81);
  EXPECT_GE(new_max_d, 4.79);

  // Re-check sanity should now be Optimal
  const auto new_report = CutoffSanityChecker::checkCutoffs(matrix, elements_, &cell_);
  EXPECT_EQ(new_report.severity, SanitySeverity::Optimal);
}

TEST_F(CutoffSanityCheckerTests, AppliesFixResetToCovalent) {
  BondCutoffMatrix matrix(2);
  // Range inversion
  matrix[0][0].min_sq = 3.0 * 3.0;
  matrix[0][0].max_sq = 1.0 * 1.0;
  matrix[0][1].min_sq = 2.5 * 2.5;
  matrix[0][1].max_sq = 1.5 * 1.5;
  matrix[1][1].min_sq = 2.0 * 2.0;
  matrix[1][1].max_sq = 1.0 * 1.0;

  const bool applied =
      CutoffSanityChecker::applyFix(matrix, elements_, &cell_, SanityFixAction::ResetToCovalent);
  EXPECT_TRUE(applied);

  // All pairs should now have valid positive non-inverted cutoff ranges
  bool all_bounds_valid = true;
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = i; j < 2; ++j) {
      if (matrix[i][j].min_sq >= matrix[i][j].max_sq || matrix[i][j].min_sq <= 0.0) {
        all_bounds_valid = false;
      }
    }
  }
  EXPECT_TRUE(all_bounds_valid);

  // Re-check report should be Optimal
  const auto new_report = CutoffSanityChecker::checkCutoffs(matrix, elements_, &cell_);
  EXPECT_EQ(new_report.severity, SanitySeverity::Optimal);
}

} // namespace
