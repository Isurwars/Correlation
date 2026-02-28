// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <cmath>
#include <gtest/gtest.h>

#include "../include/Cell.hpp"
#include "../include/PhysicalData.hpp"

// Test fixture for the Cell class.
class _02_Cell_Tests : public ::testing::Test {
protected:
  // A standard 10x10x10 orthogonal cell created once for all tests in this
  // fixture.
  Cell orthogonal_cell_{{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}};
};

// --- Constructor and Lattice Tests ---

TEST_F(_02_Cell_Tests, DefaultConstructorInitializesCorrectly) {
  // Arrange & Act
  Cell cell;

  // Assert
  EXPECT_TRUE(cell.isEmpty());
  EXPECT_EQ(cell.atomCount(), 0);
  EXPECT_EQ(cell.volume(), 0.0);
  EXPECT_DOUBLE_EQ(cell.getEnergy(), 0.0);
}

TEST_F(_02_Cell_Tests, ParameterConstructorForCubicCell) {
  // Arrange
  const std::array<double, 6> params = {4.0, 4.0, 4.0, 90.0, 90.0, 90.0};

  // Act
  Cell cell(params);
  const auto &vectors = cell.latticeVectors();

  // Assert
  EXPECT_NEAR(linalg::norm(vectors[0]), 4.0, 1e-9);
  EXPECT_NEAR(linalg::norm(vectors[1]), 4.0, 1e-9);
  EXPECT_NEAR(linalg::norm(vectors[2]), 4.0, 1e-9);
  EXPECT_NEAR(cell.volume(), 64.0, 1e-9);
}

TEST_F(_02_Cell_Tests, VectorConstructorCalculatesParameters) {
  // Arrange & Act
  Cell cell({2.0, 0.0, 0.0}, {0.0, 3.0, 0.0}, {0.0, 0.0, 4.0});
  const auto &params = cell.lattice_parameters();

  // Assert
  EXPECT_NEAR(params[0], 2.0, 1e-9);  // a
  EXPECT_NEAR(params[1], 3.0, 1e-9);  // b
  EXPECT_NEAR(params[2], 4.0, 1e-9);  // c
  EXPECT_NEAR(params[3], 90.0, 1e-9); // alpha
  EXPECT_NEAR(params[4], 90.0, 1e-9); // beta
  EXPECT_NEAR(params[5], 90.0, 1e-9); // gamma
  EXPECT_NEAR(cell.volume(), 24.0, 1e-9);
}

TEST_F(_02_Cell_Tests, VolumeForNonOrthogonalCellIsCorrect) {
  // Arrange
  const std::array<double, 6> params = {5.0, 6.0, 7.0, 80.0, 90.0, 100.0};
  Cell cell(params);

  // Act: Expected volume from standard formula
  const double cos_a = std::cos(80.0 * constants::deg2rad);
  const double cos_b = std::cos(90.0 * constants::deg2rad);
  const double cos_g = std::cos(100.0 * constants::deg2rad);
  const double vol_sqrt = 1.0 - cos_a * cos_a - cos_b * cos_b - cos_g * cos_g +
                          2 * cos_a * cos_b * cos_g;
  const double expected_volume = 5.0 * 6.0 * 7.0 * std::sqrt(vol_sqrt);

  // Assert
  EXPECT_NEAR(cell.volume(), expected_volume, 1e-9);
}

TEST_F(_02_Cell_Tests, ThrowsOnInvalidLatticeParameters) {
  EXPECT_THROW(Cell({-5.0, 1.0, 1.0, 90.0, 90.0, 90.0}), std::invalid_argument);
  EXPECT_THROW(Cell({0.0, 1.0, 1.0, 90.0, 90.0, 90.0}), std::invalid_argument);
  EXPECT_THROW(Cell({0.0, 1.0, 1.0, 90.0, 90.0, 90.0}), std::invalid_argument);
}

TEST_F(_02_Cell_Tests, RuleOfFiveWorksCorrectly) {
  // Arrange
  const std::array<double, 6> params = {4.0, 4.0, 4.0, 90.0, 90.0, 90.0};
  Cell cell(params);
  cell.addAtom("H", {0.5, 0.5, 0.5});
  cell.setEnergy(1.23);

  // Copy Constructor
  Cell cell_copy(cell);
  EXPECT_EQ(cell_copy.volume(), cell.volume());
  EXPECT_EQ(cell_copy.atomCount(), cell.atomCount());
  EXPECT_DOUBLE_EQ(cell_copy.getEnergy(), 1.23);

  // Copy Assignment
  Cell cell_assigned;
  cell_assigned = cell;
  EXPECT_EQ(cell_assigned.volume(), cell.volume());
  EXPECT_EQ(cell_assigned.atomCount(), cell.atomCount());
  EXPECT_DOUBLE_EQ(cell_assigned.getEnergy(), 1.23);

  // Move Constructor
  Cell cell_moved(std::move(cell_copy));
  EXPECT_EQ(cell_moved.volume(), cell.volume());
  EXPECT_EQ(cell_moved.atomCount(), cell.atomCount());
  EXPECT_DOUBLE_EQ(cell_moved.getEnergy(), 1.23);

  // Move Assignment
  Cell cell_move_assigned;
  cell_move_assigned = std::move(cell_assigned);
  EXPECT_EQ(cell_move_assigned.volume(), cell.volume());
  EXPECT_EQ(cell_move_assigned.atomCount(), cell.atomCount());
  EXPECT_DOUBLE_EQ(cell_move_assigned.getEnergy(), 1.23);
}

TEST_F(_02_Cell_Tests, AccessorsWorkCorrectly) {
  // Arrange
  Cell cell;
  // Test setLatticeParameters
  const std::array<double, 6> params = {10.0, 10.0, 10.0, 90.0, 90.0, 90.0};
  cell.setLatticeParameters(params);

  // Assert Lattice Parameters and Volume
  EXPECT_NEAR(cell.volume(), 1000.0, 1e-9);

  // Test Inverse Lattice Vectors (Identity for 10x10x10 cube should be close to
  // diag(0.1, 0.1, 0.1))
  const auto &inv_vecs = cell.inverseLatticeVectors();
  EXPECT_NEAR(inv_vecs[0][0], 0.1, 1e-9);

  // Test Energy
  cell.setEnergy(-13.6);
  EXPECT_DOUBLE_EQ(cell.getEnergy(), -13.6);
}

// --- Atom & Element Management Tests ---

TEST_F(_02_Cell_Tests, AddAtomManagesElementsAndAtomsCorrectly) {
  // Arrange
  Cell cell;

  // Act
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("O", {1.0, 1.0, 1.0});
  cell.addAtom("Si", {2.0, 2.0, 2.0});

  // Assert
  const auto &elements = cell.elements();
  const auto &atoms = cell.atoms();

  ASSERT_EQ(elements.size(), 2);
  EXPECT_EQ(elements[0].symbol, "Si");
  EXPECT_EQ(elements[0].id.value, 0);
  EXPECT_EQ(elements[1].symbol, "O");
  EXPECT_EQ(elements[1].id.value, 1);

  ASSERT_EQ(atoms.size(), 3);
  EXPECT_EQ(atoms[0].element().symbol, "Si");
  EXPECT_EQ(atoms[0].element().id.value, 0);
  EXPECT_EQ(atoms[1].element().symbol, "O");
  EXPECT_EQ(atoms[1].element().id.value, 1);
  EXPECT_EQ(atoms[2].element().symbol, "Si");
  EXPECT_EQ(atoms[2].element().id.value, 0);
  EXPECT_EQ(atoms[2].element().id.value, 0);
}

TEST_F(_02_Cell_Tests, FindElementWorksCorrectly) {
  // Arrange
  Cell cell;
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  // Act
  auto element_opt = cell.findElement("Si");
  auto missing_opt = cell.findElement("Au");

  // Assert
  ASSERT_TRUE(element_opt.has_value());
  EXPECT_EQ(element_opt->symbol, "Si");
  EXPECT_FALSE(missing_opt.has_value());
}

// --- Position Manipulation Tests ---

TEST_F(_02_Cell_Tests, WrapPositionsWrapsAtomIntoCell) {
  // Arrange: Uses orthogonal_cell_ from fixture
  orthogonal_cell_.addAtom("H", {12.0, -3.0, 5.0});

  // Act
  orthogonal_cell_.wrapPositions();
  const auto &final_pos = orthogonal_cell_.atoms().front().position();

  // Assert
  EXPECT_NEAR(final_pos.x(), 2.0, 1e-9); // 12.0 mod 10.0
  EXPECT_NEAR(final_pos.y(), 7.0, 1e-9); // -3.0 mod 10.0
  EXPECT_NEAR(final_pos.z(), 5.0, 1e-9); // 5.0 mod 10.0
}

TEST_F(_02_Cell_Tests, MoveSemanticsLeavesMovedFromStateEmpty) {
  // Arrange
  const std::array<double, 6> params = {4.0, 4.0, 4.0, 90.0, 90.0, 90.0};
  Cell cell(params);
  cell.addAtom("H", {0.5, 0.5, 0.5});

  // Act
  Cell cell_moved(std::move(cell));

  // Assert
  EXPECT_TRUE(cell.isEmpty());
  EXPECT_EQ(cell.atomCount(), 0);
}

TEST_F(_02_Cell_Tests, InverseLatticeVectorsAreCorrect) {
  const std::array<double, 6> params = {2.0, 4.0, 5.0, 90.0, 90.0, 90.0};
  Cell cell(params);

  const auto &inv = cell.inverseLatticeVectors();

  // 1/2, 1/4, 1/5 for diagonal elements of inverse matrix of orthogonal cell
  EXPECT_NEAR(inv[0][0], 0.5, 1e-9);
  EXPECT_NEAR(inv[1][1], 0.25, 1e-9);
  EXPECT_NEAR(inv[2][2], 0.2, 1e-9);
  EXPECT_NEAR(inv[0][1], 0.0, 1e-9);
}

TEST_F(_02_Cell_Tests, WrapPositionsHandlesNegativeCoordinates) {
  // Arrange
  Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  cell.addAtom("H", {-15.0, -1.0, -25.0});

  // Act
  cell.wrapPositions();
  const auto &final_pos = cell.atoms().front().position();

  // Assert (-15 mod 10 = 5, -1 mod 10 = 9, -25 mod 10 = 5)
  EXPECT_NEAR(final_pos.x(), 5.0, 1e-9);
  EXPECT_NEAR(final_pos.y(), 9.0, 1e-9);
  EXPECT_NEAR(final_pos.z(), 5.0, 1e-9);
}