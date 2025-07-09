#ifndef GTEST_CELL_CPP
#define GTEST_CELL_CPP
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/Atom.hpp"
#include "../include/Cell.hpp"
#include <array>
#include <gtest/gtest.h>
#include <vector>

//----------------------------------------------------------------------------//
//----------------------------- Assert Array Near ----------------------------//
//----------------------------------------------------------------------------//
void AssertArrayNear(const std::array<double, 3> &vec1,
                     const std::array<double, 3> &vec2,
                     double abs_error = 1e-6) {
  for (size_t i = 0; i < vec1.size(); ++i) {
    EXPECT_NEAR(vec1[i], vec2[i], abs_error) << "at index " << i;
  }
}

//----------------------------------------------------------------------------//
//-------------------------- Cell Constructors Tests -------------------------//
//----------------------------------------------------------------------------//
TEST(CellConstructorTest, DefaultConstructor) {
  std::array<double, 6> expected{1.0, 1.0, 1.0, 90, 90, 90};
  Cell cell;
  ASSERT_EQ(cell.lattice_parameters(), expected);
  ASSERT_NEAR(cell.volume(), 1.0, 1e-6);
}

TEST(CellConstructorTest, LatticeParameterConstructor) {
  Cell cell({4.0, 4.0, 4.0, 90, 90, 90});

  AssertArrayNear(cell.v_a(), (std::array<double, 3>{4.0, 0.0, 0.0}));
  AssertArrayNear(cell.v_b(), (std::array<double, 3>{0.0, 4.0, 0.0}));
  AssertArrayNear(cell.v_c(), (std::array<double, 3>{0.0, 0.0, 4.0}));
  ASSERT_NEAR(cell.volume(), 64.0, 1e-6);
}

TEST(CellConstructorTest, LatticeVectorConstructor) {
  Cell cell({2.0, 0.0, 0.0}, {2.0, 2.0, 0.0}, {0.0, 0.0, 2.0});
  auto lat_param = cell.lattice_parameters();
  ASSERT_NEAR(lat_param[0], 2.0, 1e-6);
  ASSERT_NEAR(lat_param[1], 2.0 * 1.414213506, 1e-6);
  ASSERT_NEAR(lat_param[2], 2.0, 1e-6);
  ASSERT_NEAR(lat_param[3], 90.0, 1e-6);
  ASSERT_NEAR(lat_param[4], 90.0, 1e-6);
  ASSERT_NEAR(lat_param[5], 45.0, 1e-6);
  ASSERT_NEAR(cell.volume(), 8.0, 1e-6);
}

TEST(CellConstructorTest, LatticeVectorConstructor_2) {
  Vector3D v1 = {4.0, 0.0, 0.0};
  Vector3D v2 = {0.0, 4.0, 0.0};
  Vector3D v3 = {0.0, 0.0, 4.0};

  Cell cell(v1, v2, v3);

  ASSERT_NEAR(cell.volume(), 64.0, 1e-6);
  AssertArrayNear(cell.v_a(), v1);
  AssertArrayNear(cell.v_b(), v2);
  AssertArrayNear(cell.v_c(), v3);
}

//----------------------------------------------------------------------------//
//------------------------- Basic Functionality Tests ------------------------//
//----------------------------------------------------------------------------//
TEST(CellVolumeTest, VolumeCalculationNonOrthogonal) {
  Cell cell;
  cell.setLatticeParameters({5.0, 6.0, 7.0, 80, 90, 100});
  cell.calculateLatticeVectors();
  double expected_vol = cell.volume();
  ASSERT_GT(expected_vol, 0);
  ASSERT_TRUE(std::abs(expected_vol - 5 * 6 * 7) > 1.0);
}

TEST(CellVolumeTest, InvalidLatticeVector) {
  Vector3D zero = {0.0, 0.0, 0.0};
  Vector3D valid = {1.0, 0.0, 0.0};
  EXPECT_THROW(Cell(zero, valid, valid), std::invalid_argument);
  EXPECT_THROW(Cell(valid, zero, valid), std::invalid_argument);
  EXPECT_THROW(Cell(valid, valid, zero), std::invalid_argument);
}

TEST(CellVolumeTest, ThrowsOnInvalidVolume) {
  // Test various invalid volume scenarios

  // Negative volume
  EXPECT_THROW(Cell cell({-5.0, 1.0, 1.0, 90.0, 90.0, 90.0}), std::logic_error);

  // Zero volume
  EXPECT_THROW(Cell cell({0.0, 0.0, 0.0, 90.0, 90.0, 90.0}), std::logic_error);
}

//----------------------------------------------------------------------------//
//------------------------- Atom & Element Management Tests ------------------//
//----------------------------------------------------------------------------//
TEST(CellElementTest, AddElementAndAtoms) {
  Cell cell;
  cell.addElement("Si");
  cell.addElement("O");
  cell.addAtom(Atom("Si", {0.0, 0.0, 0.0}));
  cell.addAtom(Atom("O", {1.0, 1.0, 1.0}));
  cell.addAtom(Atom("Si", {2.0, 2.0, 2.0}));

  cell.populateElementID();
  cell.calculateElementNumbers();

  ASSERT_EQ(cell.elements().size(), 2);
  ASSERT_EQ(cell.element_numbers()[0], 1);
  ASSERT_EQ(cell.element_numbers()[1], 2);
  ASSERT_EQ(cell.atoms()[0].element_id(), 1);
  ASSERT_EQ(cell.atoms()[1].element_id(), 0);
}

TEST(CellElementTest, DuplicateElementHandling) {
  Cell cell;
  cell.addElement("H");
  cell.addElement("H");
  cell.addElement("C");
  ASSERT_EQ(cell.elements().size(), 2);
  ASSERT_EQ(cell.elements()[0], "C");
  ASSERT_EQ(cell.elements()[1], "H");
}

//----------------------------------------------------------------------------//
//-------------------------- Position Correction Tests -----------------------//
//----------------------------------------------------------------------------//
TEST(CellPositionTest, PositionWrapping) {
  Cell cell({10.0, 10.0, 10.0, 90, 90, 90});
  cell.calculateLatticeVectors();
  std::vector<Atom> atoms{Atom("H", {12.0, -3.0, 5.0}, 0)};
  cell.setAtoms(atoms);
  cell.populateBondLength(1.2);
  cell.correctPositions();

  auto pos = cell.atoms()[0].position();
  EXPECT_NEAR(pos[0], 2.0, 1e-6);
  EXPECT_NEAR(pos[1], 7.0, 1e-6);
  EXPECT_NEAR(pos[2], 5.0, 1e-6);
}

TEST(CellPositionTest, PositionWrappingNonOrthogonal) {
  Cell cell({10.0, 10.0, 10.0, 80, 90, 100});
  cell.calculateLatticeVectors();

  // Add atom near boundary
  Atom atom("O", {10.5, 10.5, 10.5}, 0);
  cell.addAtom(atom);
  cell.populateBondLength(1.2);
  cell.correctPositions();

  auto pos = cell.atoms()[0].position();
  // Should wrap back into [0, 10) in fractional coordinates
  EXPECT_LT(pos[0], 10.0);
  EXPECT_LT(pos[1], 10.0);
  EXPECT_LT(pos[2], 10.0);
}

TEST(CellPositionTest, EdgePositionWrapping) {
  const double edge = 10.0;
  Cell cell({edge, edge, edge, 90, 90, 90});
  Atom edgeAtom("O", {edge + 1e-7, edge + 1e-7, edge + 1e-7});
  cell.addAtom(edgeAtom);
  cell.correctPositions();

  auto pos = cell.atoms()[0].position();
  EXPECT_NEAR(pos[0], 0.0, 1e-6);
  EXPECT_NEAR(pos[1], 0.0, 1e-6);
  EXPECT_NEAR(pos[2], 0.0, 1e-6);
}

TEST(CellPositionTest, FractionalConversion) {
  Cell cell({10.0, 10.0, 10.0, 90, 90, 90});
  cell.calculateLatticeVectors();
  std::vector<Atom> atoms{Atom("O", {0.5, 0.5, 0.5}, 0)};
  cell.setAtoms(atoms);
  cell.populateBondLength(1.2);
  cell.correctFracPositions();

  auto pos = cell.atoms()[0].position();
  EXPECT_NEAR(pos[0], 5.0, 1e-6);
  EXPECT_NEAR(pos[1], 5.0, 1e-6);
  EXPECT_NEAR(pos[2], 5.0, 1e-6);
}

//----------------------------------------------------------------------------//
//-------------------------- Collision Detection Tests -----------------------//
//----------------------------------------------------------------------------//

TEST(CollisionTest, AtomicCollisionDetection) {
  Cell cell({10, 10, 10, 90, 90, 90});
  cell.addAtom(Atom("O", {0, 0, 0}));
  cell.addAtom(Atom("O", {0.05, 0, 0})); // Colisión
  cell.populateBondLength(1.2);
  cell.correctPositions();
  EXPECT_THROW(cell.distancePopulation(5.0, true), std::runtime_error);
}

TEST(CollisionTest, AtomicWarningDetection) {
  Cell cell({10, 10, 10, 90, 90, 90});
  cell.addAtom(Atom("O", {0, 0, 0}));
  cell.addAtom(Atom("O", {0.45, 0, 0})); // Warning
  cell.populateBondLength(1.2);
  cell.correctPositions();
  EXPECT_NO_FATAL_FAILURE(cell.distancePopulation(5.0, true));
}

//----------------------------------------------------------------------------//
//---------------------------- Bond Calculation Tests ------------------------//
//----------------------------------------------------------------------------//
TEST(CellBondTest, BondLengthPopulation) {
  Cell BondCell({14.0, 14.0, 14.0, 90, 90, 90});
  BondCell.setElements({"O", "Si"});
  BondCell.populateBondLength(1.2);

  const auto &bondMatrix = BondCell.bond_length();
  EXPECT_NEAR(bondMatrix[0][1], 1.79 * 1.2, 1e-6);
  EXPECT_NEAR(bondMatrix[1][0], 1.79 * 1.2, 1e-6);
}

TEST(CellBondTest, DistancePopulation) {
  Cell DistCell({14.0, 14.0, 14.0, 90, 90, 90});
  std::vector<Atom> atoms{Atom("Si", {1.0, 1.0, 1.0}, 0),
                          Atom("O", {2.5, 1.0, 1.0}, 1)};
  DistCell.setAtoms(atoms);
  DistCell.populateBondLength(1.2);
  DistCell.correctPositions();
  DistCell.distancePopulation(5.0, true);

  EXPECT_GT(DistCell.atoms()[0].bonded_atoms().size(), 0);
  EXPECT_GT(DistCell.atoms()[1].bonded_atoms().size(), 0);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

#endif // GTEST_CELL_CPP
