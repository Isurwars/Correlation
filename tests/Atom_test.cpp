// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/Atom.hpp"
#include "../include/Constants.hpp"
#include "../include/LinearAlgebra.hpp"
#include <gtest/gtest.h>

// Use a namespace to encapsulate tests and avoid polluting the global scope.
namespace correlation::testing {

// Use a test fixture for all tests related to the Atom class.
class AtomTest : public ::testing::Test {};

//----------------------------------------------------------------------------//
//------------------------------ Constructor Tests ---------------------------//
//----------------------------------------------------------------------------//

TEST_F(AtomTest, DefaultConstructorInitializesToDefaultState) {
  // Arrange & Act
  Atom atom;
  Atom another_atom;

  // Assert
  EXPECT_DOUBLE_EQ(atom.position()[0], 0.0);
  EXPECT_DOUBLE_EQ(atom.position()[1], 0.0);
  EXPECT_DOUBLE_EQ(atom.position()[2], 0.0);
  EXPECT_EQ(atom.id(),
            another_atom.id()); // Verifies deterministic ID generation
}

TEST_F(AtomTest, ParameterizedConstructorSetsPropertiesCorrectly) {
  // Arrange
  const linalg::Vector3<double> expected_pos{1.0, 2.5, -3.0};

  // Act
  const Atom atom("O", expected_pos);

  // Assert
  EXPECT_EQ(atom.element(), "O");
  EXPECT_NEAR(atom.position()[0], expected_pos[0], 1E-8);
  EXPECT_NEAR(atom.position()[1], expected_pos[1], 1E-8);
  EXPECT_NEAR(atom.position()[2], expected_pos[2], 1E-8);
}

//----------------------------------------------------------------------------//
//-------------------------------- Method Tests ------------------------------//
//----------------------------------------------------------------------------//

TEST_F(AtomTest, DistanceCalculatesCorrectEuclideanDistance) {
  // Arrange: A classic 3-4-5 right triangle for an easy-to-verify distance.
  const Atom atom1("H", {0.0, 0.0, 0.0});
  const Atom atom2("H", {3.0, 4.0, 0.0});

  // Act
  const double distance = atom1.distance(atom2);

  // Assert
  EXPECT_NEAR(distance, 5.0, 1e-6);
}

TEST_F(AtomTest, AngleCalculatesNinetyDegreeAngle) {
  // Arrange
  const Atom center_atom("C", {0.0, 0.0, 0.0});
  const Atom neighbor_a("H", {1.0, 0.0, 0.0}); // Vector along X-axis
  const Atom neighbor_b("H", {0.0, 1.0, 0.0}); // Vector along Y-axis

  // Act
  const double angle = center_atom.angle(neighbor_a, neighbor_b);

  // Assert
  EXPECT_NEAR(angle, constants::pi / 2.0, 1e-6);
}

TEST_F(AtomTest, AddBondedAtomUpdatesBondedList) {
  // Arrange
  Atom atom_a("H", {0.0, 0.0, 0.0});
  const Atom atom_b("O", {1.0, 0.0, 0.0}, 123); // Give it a distinct ID

  // Act
  atom_a.addBondedAtom(atom_b);
  const auto bonded_atoms = atom_a.bonded_atoms();

  // Assert
  ASSERT_EQ(bonded_atoms.size(), 1);
  EXPECT_EQ(bonded_atoms[0].id(), 123);
  EXPECT_EQ(bonded_atoms[0].element(), "O");
}

TEST_F(AtomTest, ResetPositionAndElementUpdatesProperties) {
  // Arrange
  Atom atom; // Starts with default properties
  const linalg::Vector3<double> new_pos{2.0, 3.0, 4.0};

  // Act
  atom.resetPositionAndElement("Fe", new_pos);

  // Assert
  EXPECT_EQ(atom.element(), "Fe");
  EXPECT_NEAR(atom.position()[0], new_pos[0], 1E-8);
  EXPECT_NEAR(atom.position()[1], new_pos[1], 1E-8);
  EXPECT_NEAR(atom.position()[2], new_pos[2], 1E-8);
}

TEST_F(AtomTest, AngleHandlesCoincidentAtomsGracefully) {
  // Arrange: All atoms are at the origin, creating zero vectors.
  const Atom center_atom("C", {0.0, 0.0, 0.0});
  const Atom neighbor_a("H", {0.0, 0.0, 0.0});
  const Atom neighbor_b("H", {0.0, 0.0, 0.0});

  // Act
  const double angle = center_atom.angle(neighbor_a, neighbor_b);

  // Assert: The angle for zero vectors is undefined; returning 0.0 is a
  // reasonable way to handle this edge case.
  EXPECT_DOUBLE_EQ(angle, 0.0);
}

} // namespace correlation::testing
