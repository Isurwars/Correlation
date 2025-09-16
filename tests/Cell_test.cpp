// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <cmath>
#include <gtest/gtest.h>

#include "../include/Cell.hpp"
#include "../include/PhysicalData.hpp"

// Test fixture for the Cell class.
class CellTest : public ::testing::Test {
protected:
  // A standard 10x10x10 orthogonal cell created once for all tests in this
  // fixture.
  Cell orthogonal_cell_{{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}};
};

// --- Constructor and Lattice Tests ---

TEST_F(CellTest, ParameterConstructorForCubicCell) {
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

TEST_F(CellTest, VectorConstructorCalculatesParameters) {
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

TEST_F(CellTest, VolumeForNonOrthogonalCellIsCorrect) {
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

TEST_F(CellTest, ThrowsOnInvalidLatticeParameters) {
  EXPECT_THROW(Cell({-5.0, 1.0, 1.0, 90.0, 90.0, 90.0}), std::invalid_argument);
  EXPECT_THROW(Cell({0.0, 1.0, 1.0, 90.0, 90.0, 90.0}), std::invalid_argument);
}

// --- Atom & Element Management Tests ---

TEST_F(CellTest, AddAtomManagesElementsAndAtomsCorrectly) {
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
}

// --- Position Manipulation Tests ---

TEST_F(CellTest, WrapPositionsWrapsAtomIntoCell) {
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
