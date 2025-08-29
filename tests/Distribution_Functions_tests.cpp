// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/Atom.hpp"
#include "../include/Cell.hpp"
#include "../include/DistributionFunctions.hpp"
#include <algorithm> // For std::max_element
#include <gtest/gtest.h>
#include <iterator> // For std::distance
#include <vector>

namespace correlation::testing {

// A single test fixture for all DistributionFunctions tests.
class DistributionFunctionsTest : public ::testing::Test {};

//----------------------------------------------------------------------------//
//------------------------- Coordination Number Tests
//------------------------//
//----------------------------------------------------------------------------//
TEST_F(DistributionFunctionsTest, CoordinationNumberThrowsOnEmptyCell) {
  // Arrange
  Cell empty_cell;
  DistributionFunctions actions(empty_cell);

  // Act & Assert
  EXPECT_THROW(actions.coordinationNumber(), std::logic_error);
}

TEST_F(DistributionFunctionsTest, CoordinationNumberCountsCorrectly) {
  // Arrange: A central Carbon with two Hydrogen neighbors.
  Cell cell({25.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  cell.calculateLatticeVectors();
  std::vector<Atom> atoms{Atom("C", {5.0, 2.5, 2.5}, 0),
                          Atom("H", {4.5, 2.5, 2.5}, 1),
                          Atom("H", {5.5, 2.5, 2.5}, 2)};
  cell.setAtoms(atoms);
  cell.populateBondLength(1.2);
  cell.correctPositions();
  cell.distancePopulation(5.0, true);
  DistributionFunctions actions(cell);

  // Act
  actions.coordinationNumber();
  const auto &coord_hist = actions.Z();

  // Assert:
  EXPECT_EQ(coord_hist[2][2], 1);
  EXPECT_EQ(coord_hist[3][1], 2);
}

//----------------------------------------------------------------------------//
//--------------------------- RDF Calculation Tests
//--------------------------//
//----------------------------------------------------------------------------//
TEST_F(DistributionFunctionsTest, CalculateRDFThrowsOnInvalidParameters) {
  // Arrange
  Cell cell({3.0, 3.0, 3.0, 90.0, 90.0, 90.0});
  DistributionFunctions actions(cell);

  // Act & Assert for various invalid inputs
  EXPECT_THROW(actions.calculateRDF(5.0, 0.0, false), std::invalid_argument);
  EXPECT_THROW(actions.calculateRDF(5.0, -0.1, false), std::invalid_argument);
  EXPECT_THROW(actions.calculateRDF(0.0, 0.1, false), std::invalid_argument);
  EXPECT_THROW(actions.calculateRDF(-2.5, 0.1, false), std::invalid_argument);
}

TEST_F(DistributionFunctionsTest, CalculateRDFThrowsOnEmptyCell) {
  // Arrange
  Cell empty_cell;
  DistributionFunctions actions(empty_cell);

  // Act & Assert
  EXPECT_THROW(actions.calculateRDF(2.0, 0.1, false), std::logic_error);
}

TEST_F(DistributionFunctionsTest, RDFPeakPositionIsCorrect) {
  // Arrange: Two atoms exactly 1.0 unit apart.
  Cell cell({15.0, 15.0, 15.0, 90.0, 90.0, 90.0});
  std::vector<Atom> atoms{Atom("Ar", {5.0, 5.0, 5.0}, 0),
                          Atom("Ar", {6.0, 5.0, 5.0}, 1)};
  cell.setAtoms(atoms);
  cell.distancePopulation(5.0, true);
  DistributionFunctions actions(cell);
  const double bin_width = 0.2;

  // Act
  actions.calculateRDF(6.0, bin_width);
  const auto &rdf = actions.g();
  const auto &total_rdf =
      rdf.back(); // Assuming the last entry is the total g(r)

  // Assert: Find the peak of the RDF and verify its position.
  // The distance is 1.0, so the peak should be in bin 1.0 / 0.2 = 5.
  auto max_it = std::max_element(total_rdf.begin(), total_rdf.end());
  size_t peak_index = std::distance(total_rdf.begin(), max_it);

  EXPECT_EQ(peak_index, 5);
}

//----------------------------------------------------------------------------//
//---------------------------- PAD Calculation Tests
//-------------------------//
//----------------------------------------------------------------------------//
TEST_F(DistributionFunctionsTest, CalculatePADThrowsOnInvalidParameters) {
  // Arrange
  Cell cell({3.0, 3.0, 3.0, 90.0, 90.0, 90.0});
  DistributionFunctions actions(cell);

  // Act & Assert
  EXPECT_THROW(actions.calculatePAD(20.0, -1.0), std::invalid_argument);
  EXPECT_THROW(actions.calculatePAD(20.0, 0.0), std::invalid_argument);
  EXPECT_THROW(actions.calculatePAD(-20.0, 1.0), std::invalid_argument);
  EXPECT_THROW(actions.calculatePAD(0.0, 1.0), std::invalid_argument);
}

TEST_F(DistributionFunctionsTest, PADPeakPositionIsCorrectForWaterMolecule) {
  // Arrange: A water-like structure with a known 120-degree angle.
  Cell cell({25.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  cell.addAtom(Atom("O", {0.0, 0.0, 0.0}, 0));
  cell.addAtom(Atom("H", {1.0, 0.0, 0.0}, 1)); // H-O-H angle is 120 deg
  cell.addAtom(Atom("H", {-0.5, 0.866, 0.0}, 2));
  cell.distancePopulation(2.0, true);
  cell.planeAnglePopulation(true);
  DistributionFunctions actions(cell);

  // Act
  actions.calculatePAD(180.0, 1.0); // Use 1-degree bins
  const auto &pad = actions.F();
  const auto &h_o_h_pad =
      pad.back(); // Assuming last is total or relevant partial

  // Assert: The peak of the angle distribution should be at 120 degrees.
  auto max_it = std::max_element(h_o_h_pad.begin(), h_o_h_pad.end());
  size_t peak_index = std::distance(h_o_h_pad.begin(), max_it);

  EXPECT_EQ(peak_index, 120);
}

//----------------------------------------------------------------------------//
//---------------------------- SQ Calculation Tests
//--------------------------//
//----------------------------------------------------------------------------//
TEST_F(DistributionFunctionsTest, CalculateSQMatchesKnownFCCStructure) {
  // Arrange: A simple FCC cell of Palladium.
  Cell fcc_cell({4.0, 4.0, 4.0, 90.0, 90.0, 90.0});
  std::vector<Atom> atoms{
      Atom("Pd", {0.0, 0.0, 0.0}, 0), Atom("Pd", {2.0, 2.0, 0.0}, 1),
      Atom("Pd", {0.0, 2.0, 2.0}, 2), Atom("Pd", {2.0, 0.0, 2.0}, 3)};
  fcc_cell.setAtoms(atoms);
  fcc_cell.distancePopulation(20.0, true);
  DistributionFunctions actions(fcc_cell);

  // Act
  actions.calculateRDF(20.0, 0.05);
  actions.calculateSQ(25.0, 0.05, 9.0);

  // Assert: Check against known, pre-calculated values for this structure.
  // This serves as a regression test.
  EXPECT_NEAR(actions.g().back()[56], 38.29, 0.1);
  EXPECT_NEAR(actions.S().back()[87], 2.46184, 0.1);
}

} // namespace correlation::testing
