#include "readers/CifReader.hpp"
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>

using namespace correlation::readers;

TEST(CifReaderTests, Properties) {
  CifReader reader;
  EXPECT_EQ(reader.getName(), "CIF");
  EXPECT_FALSE(reader.isTrajectory());
  auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 1);
  EXPECT_EQ(exts[0], "cif");
}

TEST(CifReaderTests, ReadsStructureWithSymmetryAndUncertainty) {
  std::ofstream out("temp_cif.cif");
  out << "data_NaCl\n"
      << "_cell_length_a 5.64(2)\n" // test uncertainty parentheses
      << "_cell_length_b 5.64\n"
      << "_cell_length_c 5.64\n"
      << "_cell_angle_alpha 90(1)\n"
      << "_cell_angle_beta 90\n"
      << "_cell_angle_gamma 90\n"
      << "loop_\n"
      << "_space_group_symop_operation_xyz\n" // test alternative symmetry header
      << " 'x, y, z'\n"
      << " '-x, -y, -z'\n"
      << " 'x+1/2, y, z'\n"
      << "loop_\n"
      << "_atom_site_type_symbol\n"
      << "_atom_site_fract_x\n"
      << "_atom_site_fract_y\n"
      << "_atom_site_fract_z\n"
      << " Na 0.0 0.0 0.0\n"
      << " Cl 0.25 0.25 0.25\n"
      << " Na 0.0 0.0 0.0\n"; // duplicate atom - should be ignored
  out.close();

  CifReader reader;
  auto cell = reader.readStructure("temp_cif.cif");
  std::remove("temp_cif.cif");

  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[0], 5.64);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[3], 90.0);

  // Na 0.0 0.0 0.0 with symmetry:
  // (x, y, z) -> 0.0 0.0 0.0
  // (-x, -y, -z) -> 0.0 0.0 0.0 (duplicate, skipped)
  // (x+1/2, y, z) -> 0.5 0.0 0.0
  // Cl 0.25 0.25 0.25 with symmetry:
  // (x, y, z) -> 0.25 0.25 0.25
  // (-x, -y, -z) -> -0.25 -0.25 -0.25 -> 0.75 0.75 0.75
  // (x+1/2, y, z) -> 0.75 0.25 0.25
  // Duplicate Na input (0 0 0) should be skipped entirely.
  // Total atoms: 2 Na (0 0 0 and 0.5 0 0) + 3 Cl = 5 atoms.
  EXPECT_EQ(cell.atomCount(), 5);
}

TEST(CifReaderTests, ThrowsOnReadTrajectory) {
  CifReader reader;
  EXPECT_THROW(reader.readTrajectory("dummy.cif"), std::runtime_error);
}
