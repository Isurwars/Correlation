/**
 * @file ConstantsTests.cpp
 * @brief Unit tests for correlation::math constants and physical conversion factors.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "math/Constants.hpp"

#include <gtest/gtest.h>

namespace correlation::math::testing {

TEST(ConstantsTests, PiValuesAndDerivedMultiples) {
  EXPECT_NEAR(pi, 3.14159265358979323846, 1e-6);
  EXPECT_NEAR(two_pi, 2.0 * pi, 1e-6);
  EXPECT_NEAR(four_pi, 4.0 * pi, 1e-6);
}

TEST(ConstantsTests, AngularConversionsRoundTrip) {
  constexpr real_t angle_deg = static_cast<real_t>(45.0);
  real_t angle_rad = angle_deg * deg_to_rad;
  EXPECT_NEAR(angle_rad, pi / 4.0, 1e-6);

  real_t back_to_deg = angle_rad * rad_to_deg;
  EXPECT_NEAR(back_to_deg, angle_deg, 1e-5);
}

TEST(ConstantsTests, LengthUnitConversions) {
  constexpr real_t length_angstrom = static_cast<real_t>(1.0);
  real_t length_bohr = length_angstrom * angstrom_to_bohr;
  real_t back_to_angstrom = length_bohr * bohr_to_angstrom;

  EXPECT_NEAR(back_to_angstrom, length_angstrom, 1e-6);
  EXPECT_NEAR(bohr_to_angstrom * angstrom_to_bohr, static_cast<real_t>(1.0), 1e-6);
}

TEST(ConstantsTests, FrequencyEnergyConversions) {
  EXPECT_GT(thz_to_cminv, static_cast<real_t>(33.0));
  EXPECT_LT(thz_to_cminv, static_cast<real_t>(34.0));

  EXPECT_GT(thz_to_mev, static_cast<real_t>(4.0));
  EXPECT_LT(thz_to_mev, static_cast<real_t>(4.2));

  EXPECT_GT(kb_ev_per_k, static_cast<real_t>(0.0));
  EXPECT_GT(hbar_ev_ps, static_cast<real_t>(0.0));
}

} // namespace correlation::math::testing
