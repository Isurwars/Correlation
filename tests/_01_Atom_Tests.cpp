// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Atom.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"

#include <gtest/gtest.h>

// Test fixture for the correlation::core::Atom class and related free
// functions.
class _01_Atom_Tests : public ::testing::Test {};

TEST_F(_01_Atom_Tests, DefaultConstructorInitializesCorrectly) {
  // Arrange & Act
  const correlation::core::Atom atom{};

  // Assert
  EXPECT_EQ(atom.id(), 0);
  EXPECT_EQ(atom.element().symbol, "");
  EXPECT_EQ(atom.element().id.value, -1);
  EXPECT_DOUBLE_EQ(correlation::math::norm(atom.position()), 0.0);
}

TEST_F(_01_Atom_Tests, AccessorsModifyStateCorrectly) {
  // Arrange
  correlation::core::Atom atom;
  const correlation::core::Element new_element = {"Fe", {26}};
  const correlation::math::Vector3<double> new_pos = {1.0, 2.0, 3.0};

  // Act
  atom.setID(42);
  atom.setElement(new_element);
  atom.setPosition(new_pos);

  // Assert
  EXPECT_EQ(atom.id(), 42);
  EXPECT_EQ(atom.element().symbol, "Fe");
  EXPECT_EQ(atom.element().id.value, 26);
  EXPECT_EQ(atom.element_id(), 26);
  EXPECT_NEAR(atom.position().x(), 1.0, 1e-9);
  EXPECT_NEAR(atom.position().y(), 2.0, 1e-9);
  EXPECT_NEAR(atom.position().z(), 3.0, 1e-9);
}

TEST_F(_01_Atom_Tests, ParameterizedConstructorSetsProperties) {
  // Arrange
  const correlation::core::Element element = {"O", {1}};
  const correlation::math::Vector3<double> expected_pos = {1.0, 2.5, -3.0};
  const correlation::core::AtomID expected_id = 123;

  // Act
  const correlation::core::Atom atom(element, expected_pos, expected_id);

  // Assert
  EXPECT_EQ(atom.element().symbol, "O");
  EXPECT_EQ(atom.element().id.value, 1);
  EXPECT_EQ(atom.id(), 123);
  EXPECT_NEAR(atom.position().x(), expected_pos.x(), 1e-9);
  EXPECT_NEAR(atom.position().y(), expected_pos.y(), 1e-9);
  EXPECT_NEAR(atom.position().z(), expected_pos.z(), 1e-9);
}

TEST_F(_01_Atom_Tests, DistanceFunctionCalculatesCorrectly) {
  // Arrange: A classic 3-4-5 right triangle.
  const correlation::core::Element element = {"H", {0}};
  const correlation::core::Atom atom1(element, {0.0, 0.0, 0.0}, 0);
  const correlation::core::Atom atom2(element, {3.0, 4.0, 0.0}, 1);

  // Act
  const double d = distance(atom1, atom2);

  // Assert
  EXPECT_NEAR(d, 5.0, 1e-9);
}

TEST_F(_01_Atom_Tests, AngleFunctionCalculatesNinetyDegrees) {
  // Arrange
  const correlation::core::Element element_c = {"C", {0}};
  const correlation::core::Element element_h = {"H", {1}};
  const correlation::core::Atom center(element_c, {0.0, 0.0, 0.0}, 0);
  const correlation::core::Atom neighbor_a(element_h, {1.0, 0.0, 0.0},
                                           1); // Along X-axis
  const correlation::core::Atom neighbor_b(element_h, {0.0, 1.0, 0.0},
                                           2); // Along Y-axis

  // Act
  const double calculated_angle = angle(center, neighbor_a, neighbor_b);

  // Assert
  EXPECT_NEAR(calculated_angle, correlation::math::pi / 2.0, 1e-9);
}

TEST_F(_01_Atom_Tests, AngleFunctionHandlesCollinearAtoms) {
  // Arrange
  const correlation::core::Element element = {"O", {0}};
  const correlation::core::Atom center(element, {0.0, 0.0, 0.0}, 0);
  const correlation::core::Atom neighbor_a(element, {1.0, 0.0, 0.0}, 1);
  const correlation::core::Atom neighbor_b(element, {-1.0, 0.0, 0.0},
                                           2); // 180 degrees

  // Act
  const double calculated_angle = angle(center, neighbor_a, neighbor_b);

  // Assert
  EXPECT_NEAR(calculated_angle, correlation::math::pi, 1e-9);
}

TEST_F(_01_Atom_Tests, AngleFunctionHandlesCoincidentAtoms) {
  // Arrange: All atoms at the same position.
  const correlation::core::Element element = {"C", {0}};
  const correlation::core::Atom center(element, {0.0, 0.0, 0.0}, 0);
  const correlation::core::Atom neighbor_a(element, {0.0, 0.0, 0.0}, 1);

  // Act
  const double calculated_angle = angle(center, neighbor_a, center);

  // Assert: Angle with a zero-length vector should gracefully return 0.
  EXPECT_DOUBLE_EQ(calculated_angle, 0.0);
}

TEST_F(_01_Atom_Tests, ElementIDEqualityWorks) {
  correlation::core::ElementID id1{5};
  correlation::core::ElementID id2{5};
  correlation::core::ElementID id3{10};

  EXPECT_TRUE(id1 == id2);
  EXPECT_FALSE(id1 == id3);
}

TEST_F(_01_Atom_Tests, AngleFunctionClampsFloatingPointInaccuracies) {
  const correlation::core::Element element = {"O", {0}};
  const correlation::core::Atom center(element, {0.0, 0.0, 0.0}, 0);

  // Using vectors that are collinear to trigger theta=0 calculations
  // which tests the clamp if floating point gets > 1.0 slightly
  const correlation::core::Atom neighbor_a(element, {1e-8, 1e-8, 1e-8}, 1);
  const correlation::core::Atom neighbor_b(element, {1e-8, 1e-8, 1e-8}, 2);

  const double calculated_angle = angle(center, neighbor_a, neighbor_b);
  EXPECT_DOUBLE_EQ(calculated_angle, 0.0);
}
