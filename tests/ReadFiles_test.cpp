// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/Cell.hpp"
#include "../include/ReadFiles.hpp"
#include <fstream>
#include <gtest/gtest.h>

// NOTE: Added namespace for consistency with other test modules.
namespace correlation::testing {

//----------------------------------------------------------------------------//
//--------------------------- Test Fixture for File I/O ----------------------//
//----------------------------------------------------------------------------//
class ReadFilesTest : public ::testing::Test {
protected:
  // This function runs before each test to create temporary files.
  void SetUp() override {
    // Create a temporary bond file
    std::ofstream bond_file("test_bonds.txt");
    ASSERT_TRUE(bond_file.is_open());
    bond_file << "Si O  1.61\n";
    bond_file << "O  O  1.48\n";
    bond_file << "H  H  0.74\n";
    bond_file.close();

    // Create a temporary CAR file
    std::ofstream car_file("test.car");
    ASSERT_TRUE(car_file.is_open());
    car_file << "!BIOSYM archive 3\n";
    car_file << "PBC=ON\n";
    car_file << "PBC 10.5 11.5 12.5 90.0 90.0 90.0\n";
    car_file << "C1      1.00  2.00  3.00 XXXX 1      xx      C    0.000\n";
    car_file << "Si2     4.50  5.50  6.50 XXXX 1      xx      Si   0.000\n";
    car_file << "end\n";
    car_file << "end\n";
    car_file.close();

    // Create a temporary CELL file
    std::ofstream cell_file("test.cell");
    ASSERT_TRUE(cell_file.is_open());
    cell_file << "%BLOCK LATTICE_ABC\n";
    cell_file << " 15.0 15.0 20.0\n";
    cell_file << " 90.0 90.0 120.0\n";
    cell_file << "%ENDBLOCK LATTICE_ABC\n\n";
    cell_file << "%BLOCK POSITIONS_ABS\n";
    cell_file << " C   1.1 2.2 3.3\n";
    cell_file << " O   4.4 5.5 6.6\n";
    cell_file << "%ENDBLOCK POSITIONS_ABS\n";
    cell_file.close();
  }

  // This function runs after each test to clean up temporary files.
  void TearDown() override {
    remove("test_bonds.txt");
    remove("test.car");
    remove("test.cell");
  }
};

//----------------------------------------------------------------------------//
//--------------------------------- Test Cases -------------------------------//
//----------------------------------------------------------------------------//
TEST_F(ReadFilesTest, ParsesBondFileCorrectly) {
  // Arrange
  Cell test_cell;
  // NOTE: Element indices depend on the order of addElement: H=0, O=1, Si=2
  test_cell.addElement("H");
  test_cell.addElement("O");
  test_cell.addElement("Si");
  test_cell.populateBondLength(1.2); // Populate with default values first

  // Act
  auto bonds = readBond("test_bonds.txt", test_cell);

  // Assert
  ASSERT_EQ(bonds.size(), 3);
  EXPECT_DOUBLE_EQ(bonds[2][1], 1.61); // Si-O
  EXPECT_DOUBLE_EQ(bonds[1][2], 1.61); // O-Si (symmetric)
  EXPECT_DOUBLE_EQ(bonds[1][1], 1.48); // O-O
  EXPECT_DOUBLE_EQ(bonds[0][0], 0.74); // H-H
  EXPECT_NE(bonds[0][2], 0.0);         // H-Si should remain a default value
}

TEST_F(ReadFilesTest, ParsesCarFileCorrectly) {
  // Arrange & Act
  Cell result_cell = readCar("test.car");

  // Assert: Check lattice parameters
  const auto &params = result_cell.lattice_parameters();
  EXPECT_DOUBLE_EQ(params[0], 10.5);
  EXPECT_DOUBLE_EQ(params[1], 11.5);
  EXPECT_DOUBLE_EQ(params[2], 12.5);
  EXPECT_DOUBLE_EQ(params[3], 90.0);
  EXPECT_DOUBLE_EQ(params[4], 90.0);
  EXPECT_DOUBLE_EQ(params[5], 90.0);

  // Assert: Check atoms
  const auto &atoms = result_cell.atoms();
  ASSERT_EQ(atoms.size(), 2);

  EXPECT_EQ(atoms[0].element(), "C");
  EXPECT_DOUBLE_EQ(atoms[0].position()[0], 1.0);
  EXPECT_DOUBLE_EQ(atoms[0].position()[1], 2.0);
  EXPECT_DOUBLE_EQ(atoms[0].position()[2], 3.0);

  EXPECT_EQ(atoms[1].element(), "Si");
  EXPECT_DOUBLE_EQ(atoms[1].position()[0], 4.5);
  EXPECT_DOUBLE_EQ(atoms[1].position()[1], 5.5);
  EXPECT_DOUBLE_EQ(atoms[1].position()[2], 6.5);
}

TEST_F(ReadFilesTest, ParsesCellFileCorrectly) {
  // Arrange & Act
  Cell result_cell = readCell("test.cell");

  // Assert: Check lattice parameters
  const auto &params = result_cell.lattice_parameters();
  EXPECT_DOUBLE_EQ(params[0], 15.0);  // a
  EXPECT_DOUBLE_EQ(params[1], 15.0);  // b
  EXPECT_DOUBLE_EQ(params[2], 20.0);  // c
  EXPECT_DOUBLE_EQ(params[3], 90.0);  // alpha
  EXPECT_DOUBLE_EQ(params[4], 90.0);  // beta
  EXPECT_DOUBLE_EQ(params[5], 120.0); // gamma

  // Assert: Check atoms
  const auto &atoms = result_cell.atoms();
  ASSERT_EQ(atoms.size(), 2);

  EXPECT_EQ(atoms[0].element(), "C");
  EXPECT_DOUBLE_EQ(atoms[0].position()[0], 1.1);
  EXPECT_DOUBLE_EQ(atoms[0].position()[1], 2.2);
  EXPECT_DOUBLE_EQ(atoms[0].position()[2], 3.3);

  EXPECT_EQ(atoms[1].element(), "O");
  EXPECT_DOUBLE_EQ(atoms[1].position()[0], 4.4);
  // NOTE: Completed assertions for the Oxygen atom's position.
  EXPECT_DOUBLE_EQ(atoms[1].position()[1], 5.5);
  EXPECT_DOUBLE_EQ(atoms[1].position()[2], 6.6);
}

} // namespace correlation::testing
