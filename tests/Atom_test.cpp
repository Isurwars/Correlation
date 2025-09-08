// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>

#include "../include/Atom.hpp"
#include "../include/Constants.hpp"

// Test fixture for the Atom class and related free functions.
class AtomTest : public ::testing::Test {};

TEST_F(AtomTest, DefaultConstructorInitializesCorrectly) {
  // Arrange & Act
  const Atom atom{};

  // Assert
  EXPECT_EQ(atom.id().value, 0);
  EXPECT_EQ(atom.element().symbol, "");
  EXPECT_EQ(atom.element().id.value, -1);
  EXPECT_DOUBLE_EQ(linalg::norm(atom.position()), 0.0);
}

TEST_F(AtomTest, ParameterizedConstructorSetsProperties) {
  // Arrange
  const Element element = {"O", {1}};
  const linalg::Vector3<double> expected_pos = {1.0, 2.5, -3.0};
  const AtomID expected_id = {123};

  // Act
  const Atom atom(element, expected_pos, expected_id);

  // Assert
  EXPECT_EQ(atom.element().symbol, "O");
  EXPECT_EQ(atom.element().id.value, 1);
  EXPECT_EQ(atom.id().value, 123);
  EXPECT_NEAR(atom.position().x(), expected_pos.x(), 1e-9);
  EXPECT_NEAR(atom.position().y(), expected_pos.y(), 1e-9);
  EXPECT_NEAR(atom.position().z(), expected_pos.z(), 1e-9);
}

TEST_F(AtomTest, DistanceFunctionCalculatesCorrectly) {
  // Arrange: A classic 3-4-5 right triangle.
  const Element element = {"H", {0}};
  const Atom atom1(element, {0.0, 0.0, 0.0}, {0});
  const Atom atom2(element, {3.0, 4.0, 0.0}, {1});

  // Act
  const double d = distance(atom1, atom2);

  // Assert
  EXPECT_NEAR(d, 5.0, 1e-9);
}

TEST_F(AtomTest, AngleFunctionCalculatesNinetyDegrees) {
  // Arrange
  const Element element_c = {"C", {0}};
  const Element element_h = {"H", {1}};
  const Atom center(element_c, {0.0, 0.0, 0.0}, {0});
  const Atom neighbor_a(element_h, {1.0, 0.0, 0.0}, {1}); // Along X-axis
  const Atom neighbor_b(element_h, {0.0, 1.0, 0.0}, {2}); // Along Y-axis

  // Act
  const double calculated_angle = angle(center, neighbor_a, neighbor_b);

  // Assert
  EXPECT_NEAR(calculated_angle, constants::pi / 2.0, 1e-9);
}

TEST_F(AtomTest, AngleFunctionHandlesCollinearAtoms) {
  // Arrange
  const Element element = {"O", {0}};
  const Atom center(element, {0.0, 0.0, 0.0}, {0});
  const Atom neighbor_a(element, {1.0, 0.0, 0.0}, {1});
  const Atom neighbor_b(element, {-1.0, 0.0, 0.0}, {2}); // 180 degrees

  // Act
  const double calculated_angle = angle(center, neighbor_a, neighbor_b);

  // Assert
  EXPECT_NEAR(calculated_angle, constants::pi, 1e-9);
}

TEST_F(AtomTest, AngleFunctionHandlesCoincidentAtoms) {
  // Arrange: All atoms at the same position.
  const Element element = {"C", {0}};
  const Atom center(element, {0.0, 0.0, 0.0}, {0});
  const Atom neighbor_a(element, {0.0, 0.0, 0.0}, {1});

  // Act
  const double calculated_angle = angle(center, neighbor_a, center);

  // Assert: Angle with a zero-length vector should gracefully return 0.
  EXPECT_DOUBLE_EQ(calculated_angle, 0.0);
}
