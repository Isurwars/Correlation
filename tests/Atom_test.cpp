#ifndef GTEST_ATOM_CPP
#define GTEST_ATOM_CPP
/* ----------------------------------------------------------------------------
 * Correlation: An Analysis Tool for Liquids and for Amorphous Solids
 * Copyright (c) 2013-2025 Isaías Rodríguez <isurwars@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the MIT License version as published in:
 * https://github.com/Isurwars/Correlation/blob/main/LICENSE
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 * ----------------------------------------------------------------------------
 */
#include <gtest/gtest.h>
#include "../include/Atom.hpp"
#include "../include/Constants.hpp"

//----------------------------------------------------------------------------//
//------------------------------ Constructor Tests ---------------------------//
//----------------------------------------------------------------------------//

TEST(AtomTests, DefaultConstructorSetsBasicProperties) {
  Atom atom;

  std::array<double, 3> expected_pos{0.0, 0.0, 0.0};
  EXPECT_EQ(atom.position(), expected_pos);

  Atom another_atom;
  EXPECT_EQ(atom.getID(), another_atom.getID());
}

TEST(AtomTests, ParameterizedConstructorSetsProperties) {
  Atom atom("O", {1.0, 2.5, -3.0});

  std::array<double, 3> expected_pos{1.0, 2.5, -3.0};
  EXPECT_EQ(atom.position(), expected_pos);
}

//----------------------------------------------------------------------------//
//-------------------------------- Method Tests ------------------------------//
//----------------------------------------------------------------------------//

TEST(AtomTests, DistanceCalculationIsAccurate) {
  Atom a1("H", {0.0, 0.0, 0.0}, 0);
  Atom a2("H", {3.0, 4.0, 0.0}, 1);

  EXPECT_NEAR(a1.distance(a2), 5.0, 1e-6);
}

TEST(AtomTests, BondAngleCalculationHandlesBasicCases) {
  Atom center("C", {0.0, 0.0, 0.0}, 0);
  Atom a("H", {1.0, 0.0, 0.0}, 1);
  Atom b("H", {0.0, 1.0, 0.0}, 2);

  double angle = center.getAngle(a, b);
  EXPECT_NEAR(angle, constants::pi / 2, 1e-6);
}

TEST(AtomTests, AddBondedAtomUpdatesList) {
  Atom atom_a("H", {0.2, 0.0, 0.0}, 0);
  Atom atom_b("H", {1.0, 0.0, 0.0}, 123);

  atom_a.addBondedAtom(atom_b);
  std::vector<int> bonded_ids = atom_a.getBondedAtomsID();

  ASSERT_EQ(bonded_ids.size(), 1);
  EXPECT_EQ(bonded_ids[0], 123);
}

TEST(AtomTests, SetAllUpdatesProperties) {
  Atom atom;
  atom.setAll("Fe", {2.0, 3.0, 4.0});

  std::array<double, 3> expected_pos{2.0, 3.0, 4.0};
  EXPECT_EQ(atom.position(), expected_pos);
}

TEST(AtomTests, BondAngleHandlesZeroVectorsGracefully) {
  Atom center("C", {0.0, 0.0, 0.0});
  Atom a("H", {0.0, 0.0, 0.0});
  Atom b("H", {0.0, 0.0, 0.0});

  double angle = center.getAngle(a, b);
  EXPECT_DOUBLE_EQ(angle, 0.0);
}

#endif // GTEST_ATOM_CPP
