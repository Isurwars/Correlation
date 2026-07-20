// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Atom.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"

#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::core;
namespace {
class AtomTests : public ::testing::Test {};
} // namespace
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
  const correlation::math::Vector3<real_t> expected_pos = {1.0, 2.5, -3.0};
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
  const correlation::math::Vector3<real_t> new_pos = {1.0, 2.0, 3.0};

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

  EXPECT_NEAR(angle(center, neighbor_a, neighbor_b), correlation::math::pi / 2.0, correlation::is_single_precision ? 1e-6 : 1e-9);
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

  EXPECT_NEAR(angle(center, neighbor_a, neighbor_b), correlation::math::pi, correlation::is_single_precision ? 1e-6 : 1e-9);
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
  const real_t large = correlation::is_single_precision ? 1e4 : 1e10;
  const Atom atom1(element, {large, large, large}, 0);
  const Atom atom2(element, {static_cast<real_t>(large + 3.0), static_cast<real_t>(large + 4.0), large}, 1);

  EXPECT_NEAR(distance(atom1, atom2), 5.0, correlation::is_single_precision ? 1e-3 : 1e-7);
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

TEST_F(AtomTests, DistanceBetweenIdenticalAtomsIsZero) {
  const Element element = {"H", {0}};
  const Atom atom(element, {1.2, -3.4, 5.6}, 1);
  EXPECT_DOUBLE_EQ(distance(atom, atom), 0.0);
}

TEST_F(AtomTests, DistanceHandlesSubnormalCoordinates) {
  const Element element = {"H", {0}};
  const real_t subnormal = 1e-308;
  const Atom atom1(element, {0.0, 0.0, 0.0}, 0);
  const Atom atom2(element, {subnormal, subnormal, subnormal}, 1);
  real_t dist = distance(atom1, atom2);
  EXPECT_TRUE(dist == 0.0 || std::abs(dist - std::sqrt(3.0) * subnormal) < 1e-310);
}

TEST_F(AtomTests, DistanceHandlesInfinityAndNaN) {
  const Element element = {"H", {0}};
  const Atom atom1(element, {0.0, 0.0, 0.0}, 0);
  const Atom atom_inf(element, {std::numeric_limits<real_t>::infinity(), 0.0, 0.0}, 1);
  const Atom atom_nan(element, {std::numeric_limits<real_t>::quiet_NaN(), 0.0, 0.0}, 2);

  EXPECT_TRUE(std::isinf(distance(atom1, atom_inf)));
  EXPECT_TRUE(std::isnan(distance(atom1, atom_nan)));
}

TEST_F(AtomTests, AngleFunctionHandlesNaNCoordinates) {
  const Element element = {.symbol = "C", .id = {0}};
  const Atom center(element, {0.0, 0.0, 0.0}, 0);
  const Atom atom_a(element, {std::numeric_limits<real_t>::quiet_NaN(), 1.0, 1.0}, 1);
  const Atom atom_b(element, {1.0, 1.0, 1.0}, 2);

  // Since NaN coordinate makes dot/norm_sq NaN or invalid, we expect standard clamp/acos behavior or nan
  // Let's assert it safely returns 0 or NaN, without crashing.
  EXPECT_TRUE(std::isnan(angle(center, atom_a, atom_b)) || angle(center, atom_a, atom_b) == 0.0);
}

TEST_F(AtomTests, HandlesLongAndSpecialElementSymbols) {
  const Element element = {.symbol = "Uun-110_LongSymbolTest!@#", .id = {110}};
  const Atom atom(element, {1.0, 2.0, 3.0}, 999999);

  EXPECT_EQ(atom.element().symbol, "Uun-110_LongSymbolTest!@#");
  EXPECT_EQ(atom.element_id(), 110);
  EXPECT_EQ(atom.id(), 999999);
}

} // namespace correlation::testing
