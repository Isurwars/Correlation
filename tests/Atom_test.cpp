#ifndef GTEST_ATOM_CPP
#define GTEST_ATOM_CPP
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/Atom.hpp"
#include "../include/Constants.hpp"
#include "../include/LinearAlgebra.hpp"
#include <gtest/gtest.h>

//----------------------------------------------------------------------------//
//------------------------------ Constructor Tests ---------------------------//
//----------------------------------------------------------------------------//

TEST(AtomConstructorTest, DefaultConstructorSetsBasicProperties) {
  Atom atom;

  Vector3D pos{0.0, 0.0, 0.0};
  for (std::size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(atom.position()[i], pos[i], 1E-8);
  }
  Atom another_atom;
  EXPECT_EQ(atom.id(), another_atom.id());
}

TEST(AtomConstructorTest, ParameterizedConstructorSetsProperties) {
  Vector3D pos{1.0, 2.5, -3.0};
  Atom atom("O", pos);

  for (std::size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(atom.position()[i], pos[i], 1E-8);
  }
}

//----------------------------------------------------------------------------//
//-------------------------------- Method Tests ------------------------------//
//----------------------------------------------------------------------------//

TEST(AtomMethodTests, DistanceCalculationIsAccurate) {
  Atom a1("H", {0.0, 0.0, 0.0}, 0);
  Atom a2("H", {3.0, 4.0, 0.0}, 1);

  EXPECT_NEAR(a1.distance(a2), 5.0, 1e-6);
}

TEST(AtomMethodTests, BondAngleCalculationHandlesBasicCases) {
  Atom center("C", {0.0, 0.0, 0.0}, 0);
  Atom a("H", {1.0, 0.0, 0.0}, 1);
  Atom b("H", {0.0, 1.0, 0.0}, 2);

  double angle = center.getAngle(a, b);
  EXPECT_NEAR(angle, constants::pi / 2, 1e-6);
}

TEST(AtomMethodTests, AddBondedAtomUpdatesList) {
  Atom atom_a("H", {0.2, 0.0, 0.0}, 0);
  Atom atom_b("H", {1.0, 0.0, 0.0}, 123);

  atom_a.addBondedAtom(atom_b);
  std::vector<Atom> bonded = atom_a.bonded_atoms();

  ASSERT_EQ(bonded.size(), 1);
  EXPECT_EQ(bonded[0].id(), 123);
}

TEST(AtomMethodTests, SetAllUpdatesProperties) {
  Atom atom;
  atom.setAll("Fe", {2.0, 3.0, 4.0});

  Vector3D pos{2.0, 3.0, 4.0};
  for (std::size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(atom.position()[i], pos[i], 1E-8);
  }
}

TEST(AtomMethodTests, BondAngleHandlesZeroVectorsGracefully) {
  Atom center("C", {0.0, 0.0, 0.0});
  Atom a("H", {0.0, 0.0, 0.0});
  Atom b("H", {0.0, 0.0, 0.0});

  double angle = center.getAngle(a, b);
  EXPECT_DOUBLE_EQ(angle, 0.0);
}

#endif // GTEST_ATOM_CPP
