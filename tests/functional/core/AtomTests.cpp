// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Atom.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"

#include <cmath>
#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::core;

namespace {
class AtomFunctionalTests : public ::testing::Test {};

TEST_F(AtomFunctionalTests, VerifyMethaneTetrahedralGeometry) {
  // Carbon atom at the center of the tetrahedron
  const Element element_c{.symbol = "C", .id = {6}};
  Atom carbon(element_c, {0.0, 0.0, 0.0}, 0);

  // Four Hydrogen atoms at the vertices of the tetrahedron
  // C-H bond length is approx 1.09 Angstroms
  const Element element_h{.symbol = "H", .id = {1}};
  const double bond_len = 1.09;
  const double factor = bond_len / std::numbers::sqrt3;

  Atom h_1(element_h, {factor, factor, factor}, 1);
  Atom h_2(element_h, {-factor, -factor, factor}, 2);
  Atom h_3(element_h, {-factor, factor, -factor}, 3);
  Atom h_4(element_h, {factor, -factor, -factor}, 4);

  // Verify bond distances from Carbon to each Hydrogen
  EXPECT_NEAR(distance(carbon, h_1), bond_len, 1e-6);
  EXPECT_NEAR(distance(carbon, h_2), bond_len, 1e-6);
  EXPECT_NEAR(distance(carbon, h_3), bond_len, 1e-6);
  EXPECT_NEAR(distance(carbon, h_4), bond_len, 1e-6);

  // Verify tetrahedral angle: acos(-1/3) which is ~ 109.47 degrees (1.910633 radians)
  const double expected_angle_rad = std::acos(-1.0 / 3.0);
  const double expected_angle_deg = 109.4712;

  double const angle_1_c_2 = angle(carbon, h_1, h_2);
  double const angle_1_c_3 = angle(carbon, h_1, h_3);
  double const angle_1_c_4 = angle(carbon, h_1, h_4);
  double const angle_2_c_3 = angle(carbon, h_2, h_3);
  double const angle_2_c_4 = angle(carbon, h_2, h_4);
  double const angle_3_c_4 = angle(carbon, h_3, h_4);

  EXPECT_NEAR(angle_1_c_2, expected_angle_rad, 1e-6);
  EXPECT_NEAR(angle_1_c_3, expected_angle_rad, 1e-6);
  EXPECT_NEAR(angle_1_c_4, expected_angle_rad, 1e-6);
  EXPECT_NEAR(angle_2_c_3, expected_angle_rad, 1e-6);
  EXPECT_NEAR(angle_2_c_4, expected_angle_rad, 1e-6);
  EXPECT_NEAR(angle_3_c_4, expected_angle_rad, 1e-6);

  // Convert to degrees and check
  EXPECT_NEAR(angle_1_c_2 * 180.0 / correlation::math::pi, expected_angle_deg, 1e-3);
}

TEST_F(AtomFunctionalTests, VerifyDynamicAttributeModifications) {
  // Simulate an atom participating in a reaction or phase transition
  Element element = {.symbol = "Al", .id = {13}};
  Atom atom(element, {1.0, 2.0, 3.0}, 1);

  // Change type from Al to O (e.g. oxidation simulation update)
  const Element new_element = {.symbol = "O", .id = {8}};
  atom.setElement(new_element);
  atom.setID(101);
  atom.setPosition({1.5, 2.5, 3.5});
  atom.setVelocity({0.1, -0.2, 0.5});

  EXPECT_EQ(atom.element().symbol, "O");
  EXPECT_EQ(atom.element_id(), 8);
  EXPECT_EQ(atom.id(), 101);

  EXPECT_NEAR(atom.position().x(), 1.5, 1e-6);
  EXPECT_NEAR(atom.position().y(), 2.5, 1e-6);
  EXPECT_NEAR(atom.position().z(), 3.5, 1e-6);

  EXPECT_NEAR(atom.velocity().x(), 0.1, 1e-6);
  EXPECT_NEAR(atom.velocity().y(), -0.2, 1e-6);
  EXPECT_NEAR(atom.velocity().z(), 0.5, 1e-6);
}

TEST_F(AtomFunctionalTests, VerifyAngleCollinearAndOverlapping) {
  const Element element_c{.symbol = "C", .id = {6}};
  const Element element_h{.symbol = "H", .id = {1}};

  Atom center(element_c, {0.0, 0.0, 0.0}, 0);

  // Collinear
  Atom atom_a(element_h, {1.0, 0.0, 0.0}, 1);
  Atom atom_b(element_h, {-1.0, 0.0, 0.0}, 2);

  EXPECT_NEAR(angle(center, atom_a, atom_b), correlation::math::pi, 1e-6);

  // Overlapping outer atom with center
  Atom overlapping(element_h, {0.0, 0.0, 0.0}, 3);
  EXPECT_NEAR(angle(center, overlapping, atom_b), 0.0, 1e-6);

  // Overlapping center atom itself
  EXPECT_NEAR(angle(center, atom_a, center), 0.0, 1e-6);
}

} // namespace
} // namespace correlation::testing
