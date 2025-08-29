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

namespace { // Use anonymous namespace to keep helpers local to this test file.
std::string cleanToken(std::string s) {
  size_t paren_pos = s.find('(');
  if (paren_pos != std::string::npos) {
    s = s.substr(0, paren_pos);
  }
  s.erase(0, s.find_first_not_of(" \t\n\r"));
  s.erase(s.find_last_not_of(" \t\n\r") + 1);
  return s;
}

void replaceAll(std::string &str, const std::string &from,
                const std::string &to) {
  size_t start_pos = 0;
  while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length();
  }
}

Vector3D applySymmetry(const std::string &op, const Vector3D &v) {
  Vector3D result = {0.0, 0.0, 0.0};
  std::stringstream ss(op);
  std::string component;
  int i = 0;
  while (std::getline(ss, component, ',')) {
    double val = 0.0;
    double sign = 1.0;
    std::string term;
    replaceAll(component, "x", std::to_string(v[0]));
    replaceAll(component, "y", std::to_string(v[1]));
    replaceAll(component, "z", std::to_string(v[2]));
    replaceAll(component, "+", " + ");
    replaceAll(component, "-", " - ");
    std::stringstream term_ss(component);
    while (term_ss >> term) {
      if (term == "+") {
        sign = 1.0;
      } else if (term == "-") {
        sign = -1.0;
      } else {
        size_t slash_pos = term.find('/');
        if (slash_pos != std::string::npos) {
          val += sign * (std::stod(term.substr(0, slash_pos)) /
                         std::stod(term.substr(slash_pos + 1)));
        } else {
          val += sign * std::stod(term);
        }
        sign = 1.0;
      }
    }
    result[i++] = val;
  }
  return result;
}
} // namespace

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
    cif_file << " 'x+1/2, y+1/2, z'\n 'x, y+1/2, z+1/2'\n 'x+1/2, y, z+1/2'\n";
    cif_file.close();
  }

  // This function runs after each test to clean up temporary files.
  void TearDown() override {
    remove("test_bonds.txt");
    remove("test.car");
    remove("test.cell");
    remove("test.cif");
  }
};

//----------------------------------------------------------------------------//
//--------------------------------- Test Cases -------------------------------//
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//--------------------------------- Bonds File -------------------------------//
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

//----------------------------------------------------------------------------//
//---------------------------------- CAR File --------------------------------//
//----------------------------------------------------------------------------//

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
//----------------------------------------------------------------------------//
//--------------------------------- CELL File --------------------------------//
//----------------------------------------------------------------------------//

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

//----------------------------------------------------------------------------//
//------------------------ CIF Helper Function Tests -------------------------//
//----------------------------------------------------------------------------//
TEST_F(ReadFilesTest, CleanTokenHandlesUncertaintyAndWhitespace) {
  EXPECT_EQ(cleanToken(" 1.234(5) "), "1.234");
  EXPECT_EQ(cleanToken("Si1"), "Si1");
  EXPECT_EQ(cleanToken(" 90.000 "), "90.000");
}

TEST_F(ReadFilesTest, ApplySymmetryPerformsCorrectTransformations) {
  Vector3D v = {0.1, 0.2, 0.3};
  Vector3D result;

  result = applySymmetry("x, y, z", v);
  EXPECT_NEAR(result[0], 0.1, 1e-6);
  EXPECT_NEAR(result[1], 0.2, 1e-6);

  result = applySymmetry("-x, y+1/2, -z+1.0", v);
  EXPECT_NEAR(result[0], -0.1, 1e-6);
  EXPECT_NEAR(result[1], 0.7, 1e-6);
  EXPECT_NEAR(result[2], 0.7, 1e-6);
}

//----------------------------------------------------------------------------//
//--------------------------------- CIF File ---------------------------------//
//----------------------------------------------------------------------------//
TEST_F(ReadFilesTest, ParsesCifFileCorrectly) {
  // Arrange & Act
  Cell result_cell = readCif("test.cif");

  // Assert: Check lattice parameters
  const auto &params = result_cell.lattice_parameters();
  EXPECT_DOUBLE_EQ(params[0], 5.64); // a
  EXPECT_DOUBLE_EQ(params[1], 5.64); // b
  EXPECT_DOUBLE_EQ(params[2], 5.64); // c
  EXPECT_DOUBLE_EQ(params[5], 90.0); // gamma

  // Assert: Check total atom count (2 unique atoms * 4 symm ops = 8 atoms)
  const auto &atoms = result_cell.atoms();
  ASSERT_EQ(atoms.size(), 8);

  // Assert: Check a few generated atom positions (in Cartesian)
  // Na at (0,0,0) -> (0.0, 0.0, 0.0)
  // Cl at (0.5, 0.5, 0.5) -> (2.82, 2.82, 2.82)
  bool na_origin_found = false;
  bool cl_center_found = false;
  for (const auto &atom : atoms) {
    if (atom.element() == "Na" && norm(atom.position()) < 1e-4) {
      na_origin_found = true;
    }
    Vector3D expected_cl_pos = {2.82, 2.82, 2.82};
    if (atom.element() == "Cl" &&
        norm(atom.position() - expected_cl_pos) < 1e-4) {
      cl_center_found = true;
    }
  }
  EXPECT_TRUE(na_origin_found);
  EXPECT_TRUE(cl_center_found);
}

} // namespace correlation::testing
