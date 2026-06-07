// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "physics/PhysicalData.hpp"
#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::physics;

TEST(PhysicalDataTests, GetCovalentRadiusCorrectly) {
  // Test valid elements
  EXPECT_DOUBLE_EQ(getCovalentRadius("C"), 0.75);
  EXPECT_DOUBLE_EQ(getCovalentRadius("Si"), 1.16);
  EXPECT_DOUBLE_EQ(getCovalentRadius("O"), 0.63);

  // Test invalid element throws
  EXPECT_THROW(getCovalentRadius("Xx"), std::out_of_range);
}

TEST(PhysicalDataTests, GetAtomicMassCorrectly) {
  // Test valid elements
  EXPECT_DOUBLE_EQ(getAtomicMass("H"), 1.008);
  EXPECT_DOUBLE_EQ(getAtomicMass("Si"), 28.085);
  EXPECT_DOUBLE_EQ(getAtomicMass("O"), 15.999);

  // Test invalid element throws
  EXPECT_THROW(getAtomicMass("Xx"), std::out_of_range);
}

TEST(PhysicalDataTests, GetAtomicFormFactorsCorrectly) {
  // Test valid element (Silicon)
  // Silicon: {5.275329, 2.631338, 3.191038, 33.730728, 1.511514, 0.081119, 1.356849, 86.288643, 0.145073}
  auto silicon_ff = getAtomicFormFactors("Si");
  EXPECT_DOUBLE_EQ(silicon_ff[0], 5.275329);
  EXPECT_DOUBLE_EQ(silicon_ff[8], 0.145073);

  // Test invalid element throws
  EXPECT_THROW(getAtomicFormFactors("Xx"), std::out_of_range);
}

} // namespace correlation::testing
