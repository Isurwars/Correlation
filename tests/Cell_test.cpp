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

// Fixture for common test setup
class CellTest : public ::testing::Test {
protected:
  void SetUp() override {
    cubicCell.setLatticeParameters({4.0, 4.0, 4.0, 90, 90, 90});
    cubicCell.calculateLatticeVectors();
  }

  Cell cubicCell;
  Cell defaultCell;
};

//----------------------------------------------------------------------------//
//--------------------------- Basic Functionality Tests ----------------------//
//----------------------------------------------------------------------------//
TEST_F(CellTest, DefaultConstructor) {
  std::array<double, 6> expected{1.0, 1.0, 1.0, 90, 90, 90};
  ASSERT_EQ(defaultCell.lattice_parameters(), expected);
}

void AssertArrayNear(const std::array<double, 3> &vec1,
                     const std::array<double, 3> &vec2,
                     double abs_error = 1e-6) {
  for (size_t i = 0; i < vec1.size(); ++i) {
    EXPECT_NEAR(vec1[i], vec2[i], abs_error) << "at index " << i;
  }
}

TEST_F(CellTest, LatticeVectorCalculation) {
  ASSERT_NEAR(cubicCell.volume(), 64.0, 1e-6);
  AssertArrayNear(cubicCell.v_a(), (std::array<double, 3>{4.0, 0.0, 0.0}));
  AssertArrayNear(cubicCell.v_b(), (std::array<double, 3>{0.0, 4.0, 0.0}));
  AssertArrayNear(cubicCell.v_c(), (std::array<double, 3>{0.0, 0.0, 4.0}));
}

//----------------------------------------------------------------------------//
//-------------------------- Position Correction Tests -----------------------//
//----------------------------------------------------------------------------//
TEST(CellPositionTest, PositionWrapping) {
  Cell PosWarCell({10.0, 10.0, 10.0, 90, 90, 90});
  PosWarCell.calculateLatticeVectors();
  std::vector<Atom> atoms{Atom("H", {12.0, -3.0, 5.0}, 0)};
  PosWarCell.setAtoms(atoms);
  PosWarCell.populateBondLength(1.2);
  PosWarCell.correctPositions();

  auto pos = PosWarCell.atoms()[0].position();
  EXPECT_NEAR(pos[0], 2.0, 1e-6);
  EXPECT_NEAR(pos[1], 7.0, 1e-6);
  EXPECT_NEAR(pos[2], 5.0, 1e-6);
}

TEST(CellPositionTest, FractionalConversion) {
  Cell FracCell({10.0, 10.0, 10.0, 90, 90, 90});
  FracCell.calculateLatticeVectors();
  std::vector<Atom> atoms{Atom("O", {0.5, 0.5, 0.5}, 0)};
  FracCell.setAtoms(atoms);
  FracCell.populateBondLength(1.2);
  FracCell.correctFracPositions();

  auto pos = FracCell.atoms()[0].position();
  EXPECT_NEAR(pos[0], 5.0, 1e-6);
  EXPECT_NEAR(pos[1], 5.0, 1e-6);
  EXPECT_NEAR(pos[2], 5.0, 1e-6);
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

//----------------------------------------------------------------------------//
//------------------------- Coordination Number Tests ------------------------//
//----------------------------------------------------------------------------/
TEST(CoordinationErrorTest, ThrowsOnEmptyCell) {
  Cell empty_cell;
  // No atoms added
  EXPECT_THROW(empty_cell.coordinationNumber(), std::logic_error);
}

TEST(CellCoordinationTest, CoordinationCounting) {
  Cell CoordCell({25.0, 20.0, 20.0, 90, 90, 90});
  CoordCell.calculateLatticeVectors();
  std::vector<Atom> atoms{Atom("C", {5.0, 2.5, 2.5}, 0),
                          Atom("H", {4.5, 2.5, 2.5}, 1),
                          Atom("H", {5.5, 2.5, 2.5}, 2)};
  CoordCell.setAtoms(atoms);
  CoordCell.populateBondLength(1.2);
  CoordCell.correctPositions();
  CoordCell.distancePopulation(9.0, true);
  CoordCell.coordinationNumber();

  const auto coordHist = CoordCell.Z();

  EXPECT_EQ(coordHist[2][2], 1);
  EXPECT_EQ(coordHist[3][1], 2);
}

//----------------------------------------------------------------------------//
//--------------------------- RDF Calculation Tests --------------------------//
//----------------------------------------------------------------------------//
TEST(RDFErrorTest, ThrowsOnInvalidBinWidth) {
  Cell validCell({3.0, 3.0, 3.0, 90.0, 90.0, 90.0});
  // Test various invalid bin widths
  EXPECT_THROW(validCell.calculateRDF(5.0, 0.0, false), std::invalid_argument);

  EXPECT_THROW(validCell.calculateRDF(5.0, -0.1, false), std::invalid_argument);
}

TEST(RDFErrorTest, ThrowsOnInvalidRCut) {
  Cell validCell2({3.0, 3.0, 3.0, 90.0, 90.0, 90.0});
  // Test various invalid cutoff radii
  EXPECT_THROW(validCell2.calculateRDF(0.0, 0.1, false), std::invalid_argument);

  EXPECT_THROW(validCell2.calculateRDF(-2.5, 0.1, false),
               std::invalid_argument);
}

TEST(RDFErrorTest, ThrowsOnInvalidVolume) {
  // Test various invalid volume scenarios
  Cell cell({-5.0, 1.0, 1.0, 90.0, 90.0, 90.0});
  // Negative volume
  EXPECT_THROW(cell.calculateRDF(3.0, 0.1, false), std::logic_error);

  Cell cell2({0.0, 0.0, 0.0, 90.0, 90.0, 90.0});
  // Zero volume
  EXPECT_THROW(cell.calculateRDF(3.0, 0.1, false), std::logic_error);
}

TEST(RDFErrorTest, ThrowsOnEmptyCell) {
  Cell empty_cell;
  // No atoms added
  EXPECT_THROW(empty_cell.calculateRDF(2.0, 0.1, false), std::logic_error);
}
TEST(CellRDFTest, BasicRDFCalculation) {
  Cell RDFCell({4.0, 15.0, 15.0, 90, 90, 90});
  std::vector<Atom> atoms{Atom("Ar", {1.5, 2.5, 2.5}, 0),
                          Atom("Ar", {2.5, 2.5, 2.5}, 1)};
  RDFCell.setAtoms(atoms);
  RDFCell.populateBondLength(1.2);
  RDFCell.correctPositions();
  RDFCell.distancePopulation(5.0, true);
  RDFCell.calculateRDF(6.0, 0.2);

  const auto &rdf = RDFCell.g();

  EXPECT_GT(rdf[1][5], 1.0);
  EXPECT_GT(rdf[1][15], 1.0);
  EXPECT_GT(rdf[1][20], 1.0);
  EXPECT_GT(rdf[1][25], 1.0);
}

//----------------------------------------------------------------------------//
//---------------------------- Plane Angles Tests ----------------------------//
//----------------------------------------------------------------------------//
TEST(PADErrorTest, ThrowsOnInvalidBinWidth) {
  Cell validCell({3.0, 3.0, 3.0, 90.0, 90.0, 90.0});
  // Test various invalid bin widths
  EXPECT_THROW(validCell.calculatePAD(20.0, -1.0), std::invalid_argument);

  EXPECT_THROW(validCell.calculatePAD(20.0, 0.0), std::invalid_argument);
}

TEST(PADErrorTest, ThrowsOnInvalidThetaCut) {
  Cell validCell2({3.0, 3.0, 3.0, 90.0, 90.0, 90.0});
  // Test various invalid bin widths
  EXPECT_THROW(validCell2.calculatePAD(-20.0, 1.0), std::invalid_argument);

  EXPECT_THROW(validCell2.calculatePAD(0.0, 1.0), std::invalid_argument);
}

TEST(PADTest, BasicAngleCalculation) {
  Cell waterCell({25.0, 20.0, 20.0, 90, 90, 90});

  // Create water molecule-like structure
  waterCell.addAtom(Atom("O", {0, 0, 0}, 0));
  waterCell.addAtom(Atom("H", {1, 0, 0}, 1));
  waterCell.addAtom(Atom("H", {-0.5, 0.866, 0}, 2));

  // Calculate bonds first
  waterCell.populateBondLength(1.2);
  waterCell.distancePopulation(2.0, true);
  // Calculate plane angles
  waterCell.planeAnglePopulation(true);
  waterCell.calculatePAD();

  const auto &angles = waterCell.F();

  ASSERT_GT(angles[angles.size() - 1][120], 0);
}

TEST(calculateSQTest, KnownCrystalStructure) {
  Cell fccCell({4.0, 4.0, 4.0, 90.0, 90.0, 90.0});
  std::vector<Atom> atoms{
      Atom("Pd", {0.0, 0.0, 0.0}, 0), Atom("Pd", {2.0, 2.0, 0.0}, 1),
      Atom("Pd", {0.0, 2.0, 2.0}, 2), Atom("Pd", {2.0, 0.0, 2.0}, 3)};
  fccCell.setAtoms(atoms);
  fccCell.populateBondLength(1.2);
  fccCell.distancePopulation(8.0, true);
  fccCell.calculateRDF(8.0, 0.01);
  fccCell.calculateSQ(10.0, 0.01);

  const auto &g = fccCell.g();

  //  for (size_t col = 0; col < g[0].size(); ++col) {
  //    for (size_t row = 0; row < g.size(); ++row) {
  //      std::cout << g[row][col] << " ";
  //    }
  //    std::cout << '\n';
  //  }

  // Verify first peak position matches expectation
  // EXPECT_NEAR(cell.S_q()[0][first_peak_bin], expected_value, tolerance);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

#endif // GTEST_CELL_CPP
