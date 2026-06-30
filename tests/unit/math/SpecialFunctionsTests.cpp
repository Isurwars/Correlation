// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "math/SpecialFunctions.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {

using namespace correlation::math;
namespace {
class SpecialFunctionsTests : public ::testing::Test {};
} // namespace

TEST_F(SpecialFunctionsTests, FactorialCorrectness) {
  // Negative bounds
  EXPECT_DOUBLE_EQ(factorial(-1), 0.0);
  EXPECT_DOUBLE_EQ(factorial(-5), 0.0);

  // Table values
  EXPECT_DOUBLE_EQ(factorial(0), 1.0);
  EXPECT_DOUBLE_EQ(factorial(1), 1.0);
  EXPECT_DOUBLE_EQ(factorial(2), 2.0);
  EXPECT_DOUBLE_EQ(factorial(5), 120.0);
  EXPECT_DOUBLE_EQ(factorial(10), 3628800.0);
  EXPECT_DOUBLE_EQ(factorial(20), 2432902008176640000.0);

  // Gamma fallback for n > 20
  EXPECT_NEAR(factorial(21), 5.109094217170944e19, 1e10);
}

TEST_F(SpecialFunctionsTests, SphLegendreBoundaryCases) {
  // Out of bounds m
  EXPECT_DOUBLE_EQ(sph_legendre({.degree = 1, .order = -1}, 0.5), 0.0);
  EXPECT_DOUBLE_EQ(sph_legendre({.degree = 2, .order = 3}, 0.5), 0.0);
  EXPECT_DOUBLE_EQ(sph_legendre({.degree = 0, .order = 1}, 0.5), 0.0);
}

TEST_F(SpecialFunctionsTests, SphLegendreAnalyticValues) {
  const double theta = 0.5; // in radians

  // l=0, m=0: Y_0^0 = 1 / sqrt(4 * pi)
  double const y00_expected = 1.0 / std::sqrt(4.0 * M_PI);
  EXPECT_NEAR(sph_legendre({.degree = 0, .order = 0}, theta), y00_expected, 1e-9);

  // l=1, m=0: Y_1^0 = sqrt(3 / (4 * pi)) * cos(theta)
  double const y10_expected = std::sqrt(3.0 / (4.0 * M_PI)) * std::cos(theta);
  EXPECT_NEAR(sph_legendre({.degree = 1, .order = 0}, theta), y10_expected, 1e-9);

  // l=1, m=1: Y_1^1 = sqrt(3 / (8 * pi)) * sin(theta) (without Condon-Shortley phase)
  double const y11_expected = std::sqrt(3.0 / (8.0 * M_PI)) * std::sin(theta);
  EXPECT_NEAR(sph_legendre({.degree = 1, .order = 1}, theta), y11_expected, 1e-9);
}

TEST_F(SpecialFunctionsTests, SphLegendreBatchEquivalence) {
  std::vector<double> angles = {0.0, 0.2, 0.5, M_PI / 2.0, 2.1, M_PI - 0.1, M_PI};
  size_t const count = angles.size();
  std::vector<double> batch_results(count);

  // Test different combinations of l and m
  std::vector<std::pair<int, int>> const l_m_pairs = {{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {2, 2}, {3, 1}, {4, 2}};

  for (const auto &[degree, order] : l_m_pairs) {
    sph_legendre_batch({.degree = degree, .order = order}, angles.data(), batch_results.data(), count);
    for (size_t i = 0; i < count; ++i) {
      double const expected = sph_legendre({.degree = degree, .order = order}, angles[i]);
      EXPECT_NEAR(batch_results[i], expected, 1e-9)
          << "Mismatch at index " << i << " for l=" << degree << ", m=" << order;
    }
  }
}

TEST_F(SpecialFunctionsTests, SphLegendreBatchOutOfBoundsmPopulatesZero) {
  std::vector<double> angles = {0.1, 0.2, 0.3};
  std::vector<double> results(3, 1.23); // Prefill with dummy data

  sph_legendre_batch({.degree = 2, .order = 3}, angles.data(), results.data(), 3);
  EXPECT_DOUBLE_EQ(results[0], 0.0);
  EXPECT_DOUBLE_EQ(results[1], 0.0);
  EXPECT_DOUBLE_EQ(results[2], 0.0);
}

} // namespace correlation::testing
