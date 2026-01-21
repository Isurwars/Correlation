// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <fstream>
#include <gtest/gtest.h>

#include "../include/Cell.hpp"
#include "../include/FileIO.hpp"
#include "../include/LinearAlgebra.hpp"

// A single test fixture for all File I/O related tests.
// This fixture handles the creation and cleanup of temporary files needed for
// tests, ensuring that tests are self-contained and do not rely on external
// data files.
class FileIOTest : public ::testing::Test {
protected:
  // This function runs before each test to create temporary files.
  void SetUp() override {
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

    // Create a temporary CIF file for a simple rock-salt structure
    std::ofstream cif_file("test.cif");
    ASSERT_TRUE(cif_file.is_open());
    cif_file << "data_NaCl\n";
    cif_file
        << "_cell_length_a 5.64\n_cell_length_b 5.64\n_cell_length_c 5.64\n";
    cif_file
        << "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n";
    cif_file << "loop_\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_"
                "fract_y\n_atom_site_fract_z\n";
    cif_file << " Na 0.0 0.0 0.0\n Cl 0.5 0.5 0.5\n";
    cif_file << "loop_\n_symmetry_equiv_pos_as_xyz\n 'x, y, z'\n";
    cif_file << "loop_\n";
    cif_file.close();
  }

  // This function runs after each test to clean up temporary files.
  void TearDown() override {
    remove("test.car");
    remove("test.cell");
    remove("test.cif");
  }
};

//----------------------------------------------------------------------------//
//--------------------------------- Test Cases -------------------------------//
//----------------------------------------------------------------------------//

TEST_F(FileIOTest, ReadCarFileCorrectly) {
  // Arrange & Act
  FileIO::FileType type = FileIO::determineFileType("test.car");
  Cell result_cell = FileIO::readStructure("test.car", type);

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

  EXPECT_EQ(atoms[0].element().symbol, "C");
  EXPECT_DOUBLE_EQ(atoms[0].position().x(), 1.0);
  EXPECT_DOUBLE_EQ(atoms[0].position().y(), 2.0);
  EXPECT_DOUBLE_EQ(atoms[0].position().z(), 3.0);

  EXPECT_EQ(atoms[1].element().symbol, "Si");
  EXPECT_DOUBLE_EQ(atoms[1].position().x(), 4.5);
  EXPECT_DOUBLE_EQ(atoms[1].position().y(), 5.5);
  EXPECT_DOUBLE_EQ(atoms[1].position().z(), 6.5);
}

TEST_F(FileIOTest, ReadCellFileCorrectly) {
  // Arrange & Act
  FileIO::FileType type = FileIO::determineFileType("test.cell");
  Cell result_cell = FileIO::readStructure("test.cell", type);

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

  EXPECT_EQ(atoms[0].element().symbol, "C");
  EXPECT_DOUBLE_EQ(atoms[0].position().x(), 1.1);
  EXPECT_DOUBLE_EQ(atoms[0].position().y(), 2.2);
  EXPECT_DOUBLE_EQ(atoms[0].position().z(), 3.3);

  EXPECT_EQ(atoms[1].element().symbol, "O");
  EXPECT_DOUBLE_EQ(atoms[1].position().x(), 4.4);
  EXPECT_DOUBLE_EQ(atoms[1].position().y(), 5.5);
  EXPECT_DOUBLE_EQ(atoms[1].position().z(), 6.6);
}

TEST_F(FileIOTest, ReadCifFileCorrectly) {
  // Arrange & Act
  FileIO::FileType type = FileIO::determineFileType("test.cif");
  Cell result_cell = FileIO::readStructure("test.cif", type);

  // Assert: Check lattice parameters
  const auto &params = result_cell.lattice_parameters();
  EXPECT_DOUBLE_EQ(params[0], 5.64); // a
  EXPECT_DOUBLE_EQ(params[1], 5.64); // b
  EXPECT_DOUBLE_EQ(params[2], 5.64); // c
  EXPECT_DOUBLE_EQ(params[5], 90.0); // gamma

  // Assert: Check total atom count (2 unique atoms * 1 symm ops = 2 atoms)
  const auto &atoms = result_cell.atoms();
  ASSERT_EQ(atoms.size(), 2);

  // Assert: Check that specific, expected atoms were generated by symmetry
  bool na_origin_found = false;
  bool cl_center_found = false;
  for (const auto &atom : atoms) {
    if (atom.element().symbol == "Na" && linalg::norm(atom.position()) < 1e-4) {
      na_origin_found = true;
    }
    linalg::Vector3<double> expected_cl_pos = {2.82, 2.82, 2.82};
    if (atom.element().symbol == "Cl" &&
        linalg::norm(atom.position() - expected_cl_pos) < 1e-4) {
      cl_center_found = true;
    }
  }
  EXPECT_TRUE(na_origin_found)
      << "Did not find the original Na atom at (0,0,0)";
  EXPECT_TRUE(cl_center_found)
      << "Did not find the original Cl atom at the cell center";
}
