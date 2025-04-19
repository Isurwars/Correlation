#ifndef GTEST_CELL_CPP
#define GTEST_CELL_CPP
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
#include "../include/Atom.hpp"
#include "../include/Cell.hpp"
#include <array>
#include <gtest/gtest.h>
#include <vector>
#include <iostream>

// Fixture for common test setup
class CellTest : public ::testing::Test {
protected:
  void SetUp() override {
    cubicCell.setLatticeParameters({4.0, 4.0, 4.0, 90, 90, 90});
    cubicCell.setLatticeVectors();
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

void AssertArrayNear(const std::array<double, 3>& vec1,
		     const std::array<double, 3>& vec2,
		      double abs_error = 1e-6) {
    for (size_t i = 0; i < vec1.size(); ++i) {
	EXPECT_NEAR(vec1[i], vec2[i], abs_error) << "at index " << i;
    }
}

TEST_F(CellTest, LatticeVectorCalculation) {
  ASSERT_NEAR(cubicCell.getVolume(), 64.0, 1e-6);
  AssertArrayNear(cubicCell.v_a(), (std::array<double, 3>{4.0, 0.0, 0.0}));
  AssertArrayNear(cubicCell.v_b(), (std::array<double, 3>{0.0, 4.0, 0.0}));
  AssertArrayNear(cubicCell.v_c(), (std::array<double, 3>{0.0, 0.0, 4.0}));
}

//----------------------------------------------------------------------------//
//-------------------------- Position Correction Tests -----------------------//
//----------------------------------------------------------------------------//
TEST(CellPositionTest, PositionWrapping) {
  Cell PosWarCell({10.0, 10.0, 10.0, 90, 90, 90});
  PosWarCell.setLatticeVectors();
  std::vector<Atom> atoms{Atom("H", {12.0, -3.0, 5.0})};
  PosWarCell.setAtoms(atoms);
  PosWarCell.setElements({"H"});
  PosWarCell.correctPositions();

  auto pos = PosWarCell.atoms()[0].position();
  EXPECT_NEAR(pos[0], 2.0, 1e-6);
  EXPECT_NEAR(pos[1], 7.0, 1e-6);
  EXPECT_NEAR(pos[2], 5.0, 1e-6);
}

TEST(CellPositionTest, FractionalConversion) {
  Cell FracCell({10.0, 10.0, 10.0, 90, 90, 90});
  FracCell.setLatticeVectors();
  std::vector<Atom> atoms{Atom("O", {0.5, 0.5, 0.5})};
  FracCell.setAtoms(atoms);
  FracCell.setElements({"O"});
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
  BondCell.setElements({"Si", "O"});
  BondCell.populateBondLength(1.2);

  const auto &bondMatrix = BondCell.bond_length();
  EXPECT_NEAR(bondMatrix[0][1], 1.79 * 1.2, 1e-6);
  EXPECT_NEAR(bondMatrix[1][0], 1.79 * 1.2, 1e-6);
}

TEST(CellBondTest, DistancePopulation) {
  Cell DistCell({14.0, 14.0, 14.0, 90, 90, 90});
  std::vector<Atom> atoms{Atom("Si", {1.0, 1.0, 1.0}),
			  Atom("O", {2.5, 1.0, 1.0})};
  DistCell.setAtoms(atoms);
  DistCell.setElements({"Si", "O"});
  DistCell.populateBondLength(1.2);
  DistCell.distancePopulation(5.0, true);

  EXPECT_GT(DistCell.atoms()[0].bonded_atoms().size(), 0);
  EXPECT_GT(DistCell.atoms()[1].bonded_atoms().size(), 0);
}

//----------------------------------------------------------------------------//
//------------------------- Coordination Number Tests ------------------------//
//----------------------------------------------------------------------------//
TEST(CellCoordinationTest, CoordinationCounting) {
  Cell CoordCell({25.0, 20.0, 20.0, 90, 90, 90});
  std::vector<Atom> atoms{Atom("C", {5.0, 2.5, 2.5}),
			  Atom("H", {4.5, 2.5, 2.5}),
			  Atom("H", {5.5, 2.5, 2.5})};
  CoordCell.setAtoms(atoms);
  CoordCell.setElements({"C", "H"});
  CoordCell.populateBondLength(1.2);
  CoordCell.distancePopulation(9.0, true);
  CoordCell.coordinationNumber();

  const auto coordHist = CoordCell.Z();

  EXPECT_GE(coordHist[2][2], 1);
}

//----------------------------------------------------------------------------//
//--------------------------- RDF Calculation Tests --------------------------//
//----------------------------------------------------------------------------//
TEST(CellRDFTest, BasicRDFCalculation) {
  Cell RDFCell({15.0, 11.0, 13.0, 90, 90, 90});
  std::vector<Atom> atoms{Atom("Ar", {5.5, 2.5, 2.5}),
			  Atom("Ar", {6.5, 2.5, 2.5})};
  RDFCell.setAtoms(atoms);
  RDFCell.setElements({"Ar"});
  RDFCell.populateBondLength(1.2);
  RDFCell.distancePopulation(5.0, true);
  RDFCell.radialDistributionFunctions(5.0, 0.1);

  const auto &rdf = RDFCell.g();

  EXPECT_GT(rdf[1][10], rdf[1][0]);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

#endif // GTEST_CELL_CPP
