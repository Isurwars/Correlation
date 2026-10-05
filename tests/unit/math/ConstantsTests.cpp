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
  EXPECT_NEAR(PI, 3.14159265358979323846, 1e-6);
  EXPECT_NEAR(TWO_PI, 2.0 * PI, 1e-6);
  EXPECT_NEAR(FOUR_PI, 4.0 * PI, 1e-6);
}

TEST(ConstantsTests, AngularConversionsRoundTrip) {
  constexpr auto ANGLE_DEG = static_cast<real_t>(45.0);
  const real_t angle_rad = ANGLE_DEG * DEG_TO_RAD;
  EXPECT_NEAR(angle_rad, PI / 4.0, 1e-6);

  const real_t back_to_deg = angle_rad * RAD_TO_DEG;
  EXPECT_NEAR(back_to_deg, ANGLE_DEG, 1e-5);
}

TEST(ConstantsTests, LengthUnitConversions) {
  constexpr auto LENGTH_ANGSTROM = static_cast<real_t>(1.0);
  const real_t length_bohr = LENGTH_ANGSTROM * ANGSTROM_TO_BOHR;
  const real_t back_to_angstrom = length_bohr * BOHR_TO_ANGSTROM;

  EXPECT_NEAR(back_to_angstrom, LENGTH_ANGSTROM, 1e-6);
  EXPECT_NEAR(BOHR_TO_ANGSTROM * ANGSTROM_TO_BOHR, static_cast<real_t>(1.0), 1e-6);
}

TEST(ConstantsTests, FrequencyEnergyConversions) {
  EXPECT_GT(THZ_TO_CMINV, static_cast<real_t>(33.0));
  EXPECT_LT(THZ_TO_CMINV, static_cast<real_t>(34.0));

  EXPECT_GT(THZ_TO_MEV, static_cast<real_t>(4.0));
  EXPECT_LT(THZ_TO_MEV, static_cast<real_t>(4.2));

  EXPECT_GT(KB_EV_PER_K, static_cast<real_t>(0.0));
  EXPECT_GT(HBAR_EV_PS, static_cast<real_t>(0.0));
}

} // namespace correlation::math::testing
