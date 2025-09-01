// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/Atom.hpp"
#include "../include/Cell.hpp"
#include "../include/Constants.hpp"
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {

void AssertVectorNear(const linalg::Vector3<double> &vec1,
                      const linalg::Vector3<double> &vec2,
                      double abs_error = 1e-6) {
  EXPECT_NEAR(vec1[0], vec2[0], abs_error) << "at component x";
  EXPECT_NEAR(vec1[1], vec2[1], abs_error) << "at component y";
  EXPECT_NEAR(vec1[2], vec2[2], abs_error) << "at component z";
}

// A single test fixture for all Cell-related tests.
class CellTest : public ::testing::Test {
protected:
  // Use SetUp to create common objects for multiple tests.
  void SetUp() override {
    // A standard 10x10x10 orthogonal cell used in many position tests.
    orthogonal_cell_ = Cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
    orthogonal_cell_.calculateLatticeVectors();
  }

  Cell orthogonal_cell_;
};

//----------------------------------------------------------------------------//
//-------------------------- Constructor and Volume Tests --------------------//
//----------------------------------------------------------------------------//
TEST_F(CellTest, DefaultConstructorCreatesUnitCell) {
  // Arrange & Act
  Cell cell;
  const std::array<double, 6> expected_params{1.0, 1.0, 1.0, 90.0, 90.0, 90.0};

  // Assert
  EXPECT_EQ(cell.lattice_parameters(), expected_params);
  EXPECT_NEAR(cell.volume(), 1.0, 1e-6);
}

TEST_F(CellTest, LatticeParameterConstructorForCubicCell) {
  // Arrange & Act
  Cell cell({4.0, 4.0, 4.0, 90.0, 90.0, 90.0});

  // Assert
  AssertVectorNear(cell.v_a(), {4.0, 0.0, 0.0});
  AssertVectorNear(cell.v_b(), {0.0, 4.0, 0.0});
  AssertVectorNear(cell.v_c(), {0.0, 0.0, 4.0});
  EXPECT_NEAR(cell.volume(), 64.0, 1e-6);
}

TEST_F(CellTest, LatticeVectorConstructorCalculatesParameters) {
  // Arrange & Act
  Cell cell({2.0, 0.0, 0.0}, {2.0, 2.0, 0.0}, {0.0, 0.0, 2.0});
  const auto lat_param = cell.lattice_parameters();

  // Assert
  EXPECT_NEAR(lat_param[0], 2.0, 1e-6);
  EXPECT_NEAR(lat_param[1], 2.0 * std::sqrt(2.0), 1e-6); // Avoid magic numbers
  EXPECT_NEAR(lat_param[2], 2.0, 1e-6);
  EXPECT_NEAR(lat_param[3], 90.0, 1e-6); // alpha
  EXPECT_NEAR(lat_param[4], 90.0, 1e-6); // beta
  EXPECT_NEAR(lat_param[5], 45.0, 1e-6); // gamma
  EXPECT_NEAR(cell.volume(), 8.0, 1e-6);
}

TEST_F(CellTest, VolumeCalculationForNonOrthogonalCellIsCorrect) {
  // Arrange
  const double a = 5.0, b = 6.0, c = 7.0;
  const double alpha_rad = 80.0 * constants::pi / 180.0;
  const double beta_rad = 90.0 * constants::pi / 180.0;
  const double gamma_rad = 100.0 * constants::pi / 180.0;
  Cell cell;
  cell.setLatticeParameters({a, b, c, 80.0, 90.0, 100.0});
  cell.calculateLatticeVectors();

  // Act: Calculate the expected volume using the formal equation.
  const double cos_a = std::cos(alpha_rad);
  const double cos_b = std::cos(beta_rad);
  const double cos_g = std::cos(gamma_rad);
  const double volume_sqrt = 1.0 - cos_a * cos_a - cos_b * cos_b -
                             cos_g * cos_g + 2.0 * cos_a * cos_b * cos_g;
  const double expected_volume = a * b * c * std::sqrt(volume_sqrt);

  // Assert
  EXPECT_NEAR(cell.volume(), expected_volume, 1e-6);
}

TEST_F(CellTest, ConstructorThrowsOnInvalidLatticeVector) {
  // Arrange
  const linalg::Vector3<double> zero_vec = {0.0, 0.0, 0.0};
  const linalg::Vector3<double> valid_vec = {1.0, 0.0, 0.0};

  // Act & Assert
  EXPECT_THROW(Cell(zero_vec, valid_vec, valid_vec), std::invalid_argument);
  EXPECT_THROW(Cell(valid_vec, zero_vec, valid_vec), std::invalid_argument);
  EXPECT_THROW(Cell(valid_vec, valid_vec, zero_vec), std::invalid_argument);
}

TEST_F(CellTest, ConstructorThrowsOnInvalidLatticeParameters) {
  // Act & Assert
  // NOTE: Fixed incorrect EXPECT_THROW syntax by removing variable declaration.
  EXPECT_THROW(Cell({-5.0, 1.0, 1.0, 90.0, 90.0, 90.0}), std::logic_error);
  EXPECT_THROW(Cell({0.0, 1.0, 1.0, 90.0, 90.0, 90.0}), std::logic_error);
}

//----------------------------------------------------------------------------//
//------------------------- Atom & Element Management Tests ------------------//
//----------------------------------------------------------------------------//
TEST_F(CellTest, ElementAndAtomManagement) {
  // Arrange
  Cell cell;
  cell.addElement("Si");
  cell.addElement("O");
  cell.addAtom(Atom("Si", {0.0, 0.0, 0.0}));
  cell.addAtom(Atom("O", {1.0, 1.0, 1.0}));
  cell.addAtom(Atom("Si", {2.0, 2.0, 2.0}));

  // Act
  cell.populateElementID();
  cell.calculateElementNumbers();

  // Assert
  ASSERT_EQ(cell.elements().size(), 2);
  // Assuming elements are sorted alphabetically: O, Si
  ASSERT_EQ(cell.element_numbers()[0], 1);    // Count of O
  ASSERT_EQ(cell.element_numbers()[1], 2);    // Count of Si
  ASSERT_EQ(cell.atoms()[0].element_id(), 1); // First Si
  ASSERT_EQ(cell.atoms()[1].element_id(), 0); // O
}

TEST_F(CellTest, DuplicateElementIsNotAdded) {
  // Arrange
  Cell cell;
  cell.addElement("H");
  cell.addElement("H"); // Duplicate
  cell.addElement("C");

  // Act
  const auto &elements = cell.elements();

  // Assert
  ASSERT_EQ(elements.size(), 2);
  // Assuming elements are sorted:
  EXPECT_EQ(elements[0], "C");
  EXPECT_EQ(elements[1], "H");
}

//----------------------------------------------------------------------------//
//-------------------------- Position Correction Tests -----------------------//
//----------------------------------------------------------------------------//
TEST_F(CellTest, CorrectPositionsWrapsAtomIntoOrthogonalCell) {
  // Arrange: Uses orthogonal_cell_ from SetUp
  Atom atom_outside("H", {12.0, -3.0, 5.0});
  orthogonal_cell_.addAtom(atom_outside);

  // Act
  orthogonal_cell_.correctPositions();
  const auto final_pos = orthogonal_cell_.atoms()[0].position();

  // Assert
  EXPECT_NEAR(final_pos[0], 2.0, 1e-6);
  EXPECT_NEAR(final_pos[1], 7.0, 1e-6);
  EXPECT_NEAR(final_pos[2], 5.0, 1e-6);
}

TEST_F(CellTest, FractionalToCartesianConversionIsCorrect) {
  // Arrange: Uses orthogonal_cell_ from SetUp
  Atom frac_atom("O", {0.5, 0.25, 0.1}); // Fractional coordinates
  orthogonal_cell_.addAtom(frac_atom);

  // Act
  orthogonal_cell_.correctFracPositions();
  const auto final_pos = orthogonal_cell_.atoms()[0].position();

  // Assert
  EXPECT_NEAR(final_pos[0], 5.0, 1e-6); // 0.5 * 10.0
  EXPECT_NEAR(final_pos[1], 2.5, 1e-6); // 0.25 * 10.0
  EXPECT_NEAR(final_pos[2], 1.0, 1e-6); // 0.1 * 10.0
}

//----------------------------------------------------------------------------//
//---------------------------- Collision & Bond Tests ------------------------//
//----------------------------------------------------------------------------//
TEST_F(CellTest, DistancePopulationThrowsOnCollision) {
  // Arrange: Uses orthogonal_cell_ from SetUp
  orthogonal_cell_.addAtom(Atom("O", {5.0, 5.0, 5.0}));
  orthogonal_cell_.addAtom(Atom("O", {5.05, 5.0, 5.0})); // Too close
  orthogonal_cell_.populateBondLength(1.2);

  // Act & Assert
  EXPECT_THROW(orthogonal_cell_.distancePopulation(5.0, true),
               std::runtime_error);
}

TEST_F(CellTest, DistancePopulationFindsNeighbors) {
  // Arrange: Uses orthogonal_cell_ from SetUp
  Atom atom1("Si", {1.0, 1.0, 1.0}, 0);
  Atom atom2("O", {2.5, 1.0, 1.0}, 1);
  orthogonal_cell_.addAtom(atom1);
  orthogonal_cell_.addAtom(atom2);
  orthogonal_cell_.populateBondLength(1.2);

  // Act
  orthogonal_cell_.distancePopulation(5.0, true);

  // Assert
  EXPECT_FALSE(orthogonal_cell_.atoms()[0].bonded_atoms().empty());
  EXPECT_FALSE(orthogonal_cell_.atoms()[1].bonded_atoms().empty());
  EXPECT_EQ(orthogonal_cell_.atoms()[0].bonded_atoms()[0].id(), 1);
}

} // namespace correlation::testing
