// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/ml/ChargeCorrelationCalculator.hpp"
#include "calculators/ml/SpinCorrelationCalculator.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>

using correlation::real_t;

namespace correlation::calculators {
namespace {

TEST(SpinChargeCorrelationTests, SpinCalculatorFerroVsAntiferro) {
  correlation::core::Cell cell(std::array<real_t, 6>{10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Fe", {0.0, 0.0, 0.0});
  cell.addAtom("Fe", {2.5, 0.0, 0.0});

  std::vector<correlation::math::Vector3<real_t>> ferro_m{
      {0.0, 0.0, 2.0},
      {0.0, 0.0, 2.0},
  };
  auto const ferro_res = SpinCorrelationCalculator::calculate(cell, ferro_m, 5.0, 0.5);

  const auto bin_idx = static_cast<size_t>(2.5 / 0.5);
  EXPECT_EQ(ferro_res.pair_counts[bin_idx], 1U);
  EXPECT_DOUBLE_EQ(ferro_res.spin_correlation[bin_idx], 1.0);
  EXPECT_DOUBLE_EQ(ferro_res.magnetic_moment_corr[bin_idx], 4.0);

  std::vector<correlation::math::Vector3<real_t>> antiferro_m{
      {0.0, 0.0, 2.0},
      {0.0, 0.0, -2.0},
  };
  auto const anti_res = SpinCorrelationCalculator::calculate(cell, antiferro_m, 5.0, 0.5);

  EXPECT_EQ(anti_res.pair_counts[bin_idx], 1U);
  EXPECT_DOUBLE_EQ(anti_res.spin_correlation[bin_idx], -1.0);
  EXPECT_DOUBLE_EQ(anti_res.magnetic_moment_corr[bin_idx], -4.0);
}

TEST(SpinChargeCorrelationTests, ChargeCalculatorAttraction) {
  correlation::core::Cell cell(std::array<real_t, 6>{10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Na", {0.0, 0.0, 0.0});
  cell.addAtom("Cl", {2.8, 0.0, 0.0});

  std::vector<real_t> charges{1.0, -1.0};
  auto const res = ChargeCorrelationCalculator::calculate(cell, charges, 5.0, 0.5);

  const auto bin_idx = static_cast<size_t>(2.8 / 0.5);
  EXPECT_EQ(res.pair_counts[bin_idx], 1U);
  EXPECT_LT(res.g_qq[bin_idx], 0.0);
  EXPECT_LT(res.charge_integral[bin_idx], 0.0);
}

TEST(SpinChargeCorrelationTests, InvalidInputValidation) {
  correlation::core::Cell cell(std::array<real_t, 6>{10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  std::vector<real_t> charges{1.0};
  EXPECT_THROW((void)ChargeCorrelationCalculator::calculate(cell, charges, -1.0, 0.5),
               std::invalid_argument);
  EXPECT_THROW((void)ChargeCorrelationCalculator::calculate(cell, charges, 5.0, -0.1),
               std::invalid_argument);

  std::vector<real_t> empty_charges;
  EXPECT_THROW((void)ChargeCorrelationCalculator::calculate(cell, empty_charges, 5.0, 0.5),
               std::invalid_argument);
}

} // namespace
} // namespace correlation::calculators
