#ifndef GTEST_DISTRIBUTION_FUNCTIONS_CPP
#define GTEST_DISTRIBUTION_FUNCTIONS_CPP
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include <gtest/gtest.h>
#include <vector>

#include "../include/Atom.hpp"
#include "../include/Cell.hpp"
#include "../include/DistributionFunctions.hpp"
//----------------------------------------------------------------------------//
//------------------------- Coordination Number Tests ------------------------//
//----------------------------------------------------------------------------/
TEST(CoordinationErrorTest, ThrowsOnEmptyCell) {
  Cell empty_cell;
  DistributionFunctions empty_actions(empty_cell);
  // No atoms added
  EXPECT_THROW(empty_actions.coordinationNumber(), std::logic_error);
}

TEST(CellCoordinationTest, CoordinationCounting) {
  Cell coord_cell({25.0, 20.0, 20.0, 90, 90, 90});
  coord_cell.calculateLatticeVectors();
  std::vector<Atom> atoms{Atom("C", {5.0, 2.5, 2.5}, 0),
                          Atom("H", {4.5, 2.5, 2.5}, 1),
                          Atom("H", {5.5, 2.5, 2.5}, 2)};
  coord_cell.setAtoms(atoms);
  coord_cell.populateBondLength(1.2);
  coord_cell.correctPositions();
  coord_cell.distancePopulation(9.0, true);
  DistributionFunctions coord_actions(coord_cell);
  coord_actions.coordinationNumber();

  const auto coordHist = coord_actions.Z();

  EXPECT_EQ(coordHist[2][2], 1);
  EXPECT_EQ(coordHist[3][1], 2);
}

//----------------------------------------------------------------------------//
//--------------------------- RDF Calculation Tests --------------------------//
//----------------------------------------------------------------------------//
TEST(RDFErrorTest, ThrowsOnInvalidBinWidth) {
  Cell valid_cell({3.0, 3.0, 3.0, 90.0, 90.0, 90.0});
  DistributionFunctions valid(valid_cell);
  // Test various invalid bin widths
  EXPECT_THROW(valid.calculateRDF(5.0, 0.0, false), std::invalid_argument);

  EXPECT_THROW(valid.calculateRDF(5.0, -0.1, false), std::invalid_argument);
}

TEST(RDFErrorTest, ThrowsOnInvalidRCut) {
  Cell valid_cell_2({3.0, 3.0, 3.0, 90.0, 90.0, 90.0});
  DistributionFunctions valid2(valid_cell_2);
  // Test various invalid cutoff radii
  EXPECT_THROW(valid2.calculateRDF(0.0, 0.1, false), std::invalid_argument);

  EXPECT_THROW(valid2.calculateRDF(-2.5, 0.1, false), std::invalid_argument);
}

TEST(RDFErrorTest, ThrowsOnEmptyCell) {
  Cell empty_cell;
  DistributionFunctions error_actions(empty_cell);
  // No atoms added
  EXPECT_THROW(error_actions.calculateRDF(2.0, 0.1, false), std::logic_error);
}

TEST(CellRDFTest, BasicRDFCalculation) {
  Cell RDF_cell({4.0, 15.0, 15.0, 90, 90, 90});
  std::vector<Atom> atoms{Atom("Ar", {1.5, 2.5, 2.5}, 0),
                          Atom("Ar", {2.5, 2.5, 2.5}, 1)};
  RDF_cell.setAtoms(atoms);
  RDF_cell.populateBondLength(1.2);
  RDF_cell.correctPositions();
  RDF_cell.distancePopulation(5.0, true);
  DistributionFunctions RDF_actions(RDF_cell);
  RDF_actions.calculateRDF(6.0, 0.2);

  const auto &rdf = RDF_actions.g();

  EXPECT_GT(rdf[1][5], 1.0);
  EXPECT_GT(rdf[1][15], 1.0);
  EXPECT_GT(rdf[1][20], 1.0);
  EXPECT_GT(rdf[1][25], 1.0);
}

//----------------------------------------------------------------------------//
//---------------------------- Plane Angles Tests ----------------------------//
//----------------------------------------------------------------------------//
TEST(PADErrorTest, ThrowsOnInvalidBinWidth) {
  Cell valid_cell_3({3.0, 3.0, 3.0, 90.0, 90.0, 90.0});
  DistributionFunctions valid(valid_cell_3);
  // Test various invalid bin widths
  EXPECT_THROW(valid.calculatePAD(20.0, -1.0), std::invalid_argument);

  EXPECT_THROW(valid.calculatePAD(20.0, 0.0), std::invalid_argument);
}

TEST(PADErrorTest, ThrowsOnInvalidThetaCut) {
  Cell valid_cell_4({3.0, 3.0, 3.0, 90.0, 90.0, 90.0});
  DistributionFunctions valid_2(valid_cell_4);
  // Test various invalid bin widths
  EXPECT_THROW(valid_2.calculatePAD(-20.0, 1.0), std::invalid_argument);

  EXPECT_THROW(valid_2.calculatePAD(0.0, 1.0), std::invalid_argument);
}

TEST(PADTest, BasicAngleCalculation) {
  Cell water_cell({25.0, 20.0, 20.0, 90, 90, 90});

  // Create water molecule-like structure
  water_cell.addAtom(Atom("O", {0, 0, 0}, 0));
  water_cell.addAtom(Atom("H", {1, 0, 0}, 1));
  water_cell.addAtom(Atom("H", {-0.5, 0.866, 0}, 2));

  // Calculate bonds first
  water_cell.populateBondLength(1.2);
  water_cell.distancePopulation(2.0, true);
  water_cell.planeAnglePopulation(true);

  // Calculate plane angles
  DistributionFunctions water_actions(water_cell);
  water_actions.calculatePAD();

  const auto &angles = water_actions.F();

  ASSERT_GT(angles[angles.size() - 1][120], 0);
}

TEST(calculateSQTest, KnownCrystalStructure) {
  Cell fcc_cell({4.0, 4.0, 4.0, 90.0, 90.0, 90.0});
  std::vector<Atom> atoms{
      Atom("Pd", {0.0, 0.0, 0.0}, 0), Atom("Pd", {2.0, 2.0, 0.0}, 1),
      Atom("Pd", {0.0, 2.0, 2.0}, 2), Atom("Pd", {2.0, 0.0, 2.0}, 3)};
  fcc_cell.setAtoms(atoms);
  fcc_cell.populateBondLength(1.2);
  fcc_cell.distancePopulation(20.0, true);

  // Cell ACTIONS
  DistributionFunctions fcc_actions(fcc_cell);
  fcc_actions.calculateRDF(20.0, 0.05);
  fcc_actions.calculateSQ(25.0, 0.05, 9.0);

  const auto &g = fcc_actions.g();

  for (size_t col = 0; col < g[0].size(); ++col) {
    for (size_t row = 0; row < g.size(); ++row) {
      std::cout << g[row][col] << " ";
    }
    std::cout << '\n';
  }

  const auto &s = fcc_actions.S();

  for (size_t col = 0; col < s[0].size(); ++col) {
    for (size_t row = 0; row < s.size(); ++row) {
      std::cout << s[row][col] << " ";
    }
    std::cout << '\n';
  }

  // Verify first peak position matches expectation
  // EXPECT_NEAR(cell.S_q()[0][first_peak_bin], expected_value, tolerance);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

#endif // GTEST_DISTRIBUTION_FUNCTIONS_CPP
