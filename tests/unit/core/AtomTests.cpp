// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Atom.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"

#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::core;

class AtomTests : public ::testing::Test {};

// --- Unitary Tests: Constructors & Accessors ---

TEST_F(AtomTests, DefaultConstructorInitializesCorrectly) {
  const Atom atom{};
  EXPECT_EQ(atom.id(), 0);
  EXPECT_EQ(atom.element().symbol, "");
  EXPECT_EQ(atom.element().id.value, -1);
  EXPECT_DOUBLE_EQ(correlation::math::norm(atom.position()), 0.0);
}

TEST_F(AtomTests, ParameterizedConstructorSetsProperties) {
  const Element element = {"O", {1}};
  const correlation::math::Vector3<double> expected_pos = {1.0, 2.5, -3.0};
  const AtomID expected_id = 123;

  const Atom atom(element, expected_pos, expected_id);

  EXPECT_EQ(atom.element().symbol, "O");
  EXPECT_EQ(atom.element().id.value, 1);
  EXPECT_EQ(atom.id(), 123);
  EXPECT_NEAR(atom.position().x(), expected_pos.x(), 1e-9);
  EXPECT_NEAR(atom.position().y(), expected_pos.y(), 1e-9);
  EXPECT_NEAR(atom.position().z(), expected_pos.z(), 1e-9);
}

TEST_F(AtomTests, AccessorsModifyStateCorrectly) {
  Atom atom;
  const Element new_element = {"Fe", {26}};
  const correlation::math::Vector3<double> new_pos = {1.0, 2.0, 3.0};

  atom.setID(42);
  atom.setElement(new_element);
  atom.setPosition(new_pos);

  EXPECT_EQ(atom.id(), 42);
  EXPECT_EQ(atom.element().symbol, "Fe");
  EXPECT_EQ(atom.element().id.value, 26);
  EXPECT_EQ(atom.element_id(), 26);
  EXPECT_NEAR(atom.position().x(), 1.0, 1e-9);
  EXPECT_NEAR(atom.position().y(), 2.0, 1e-9);
  EXPECT_NEAR(atom.position().z(), 3.0, 1e-9);
}

// --- Unitary Tests: Geometric Functions ---

TEST_F(AtomTests, DistanceFunctionCalculatesCorrectly) {
  const Element element = {"H", {0}};
  const Atom atom1(element, {0.0, 0.0, 0.0}, 0);
  const Atom atom2(element, {3.0, 4.0, 0.0}, 1);

  EXPECT_NEAR(distance(atom1, atom2), 5.0, 1e-9);
}

TEST_F(AtomTests, AngleFunctionCalculatesNinetyDegrees) {
  const Element element = {"C", {0}};
  const Atom center(element, {0.0, 0.0, 0.0}, 0);
  const Atom neighbor_a(element, {1.0, 0.0, 0.0}, 1);
  const Atom neighbor_b(element, {0.0, 1.0, 0.0}, 2);

  EXPECT_NEAR(angle(center, neighbor_a, neighbor_b),
              correlation::math::pi / 2.0, 1e-9);
}

// --- Limit Cases ---

TEST_F(AtomTests, AngleFunctionHandlesCoincidentAtoms) {
  const Element element = {"C", {0}};
  const Atom center(element, {0.0, 0.0, 0.0}, 0);
  const Atom neighbor_a(element, {0.0, 0.0, 0.0}, 1);

  // Angle with a zero-length vector should return 0.0
  EXPECT_DOUBLE_EQ(angle(center, neighbor_a, center), 0.0);
}

TEST_F(AtomTests, AngleFunctionHandlesCollinearAtoms) {
  const Element element = {"O", {0}};
  const Atom center(element, {0.0, 0.0, 0.0}, 0);
  const Atom neighbor_a(element, {1.0, 0.0, 0.0}, 1);
  const Atom neighbor_b(element, {-1.0, 0.0, 0.0}, 2); // 180 degrees

  EXPECT_NEAR(angle(center, neighbor_a, neighbor_b), correlation::math::pi,
              1e-9);
}

TEST_F(AtomTests, AngleFunctionClampsFloatingPointInaccuracies) {
  const Element element = {"O", {0}};
  const Atom center(element, {0.0, 0.0, 0.0}, 0);

  // Vectors that are slightly offset due to precision but logically collinear
  const Atom neighbor_a(element, {1e-15, 1e-15, 1e-15}, 1);
  const Atom neighbor_b(element, {2e-15, 2e-15, 2e-15}, 2);

  EXPECT_DOUBLE_EQ(angle(center, neighbor_a, neighbor_b), 0.0);
}

TEST_F(AtomTests, DistanceHandlesLargeCoordinates) {
  const Element element = {"H", {0}};
  const double large = 1e10;
  const Atom atom1(element, {large, large, large}, 0);
  const Atom atom2(element, {large + 3.0, large + 4.0, large}, 1);

  EXPECT_NEAR(distance(atom1, atom2), 5.0, 1e-7);
}

TEST_F(AtomTests, HandlesEmptyElementSymbol) {
  const Element element = {"", {-1}};
  const Atom atom(element, {0.0, 0.0, 0.0}, 0);

  EXPECT_EQ(atom.element().symbol, "");
  EXPECT_EQ(atom.element_id(), -1);
}

TEST_F(AtomTests, ElementIDEqualityWorks) {
  ElementID id1{5};
  ElementID id2{5};
  ElementID id3{10};

  EXPECT_TRUE(id1 == id2);
  EXPECT_FALSE(id1 == id3);
}

} // namespace correlation::testing
