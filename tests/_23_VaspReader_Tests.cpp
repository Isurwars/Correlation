// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "readers/VaspReader.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <stdexcept>

namespace {

// Helper to find the test data directory relative to the executable
std::string getTestDataDir() {
  // Try common paths relative to the build directory
  std::vector<std::string> candidates = {
      "../../tests/data/",   // build/tests -> tests/data
      "../tests/data/",      // build -> tests/data
      "tests/data/",         // project root
      "data/",               // tests/data (if cwd is tests)
  };
  for (const auto &dir : candidates) {
    if (std::filesystem::exists(dir + "Si.poscar")) {
      return dir;
    }
  }
  // Try absolute path from source dir (CMake sets CWD to build/tests)
  return "../../tests/data/";
}

} // namespace

class _23_VaspReader_Tests : public ::testing::Test {
protected:
  std::string data_dir_;
  void SetUp() override { data_dir_ = getTestDataDir(); }
};

TEST_F(_23_VaspReader_Tests, ParseSiDiamondPoscar) {
  auto cell = correlation::readers::VaspReader::read(data_dir_ + "Si.poscar");

  // 8 Si atoms
  EXPECT_EQ(cell.atomCount(), 8);

  // 1 element type
  EXPECT_EQ(cell.elements().size(), 1);
  EXPECT_EQ(cell.elements()[0].symbol, "Si");

  // Lattice should be 5.43 x 5.43 x 5.43 (scaling factor * unit vectors)
  auto params = cell.lattice_parameters();
  EXPECT_NEAR(params[0], 5.43, 1e-6);
  EXPECT_NEAR(params[1], 5.43, 1e-6);
  EXPECT_NEAR(params[2], 5.43, 1e-6);
  EXPECT_NEAR(params[3], 90.0, 1e-6);
  EXPECT_NEAR(params[4], 90.0, 1e-6);
  EXPECT_NEAR(params[5], 90.0, 1e-6);
}

TEST_F(_23_VaspReader_Tests, ParseCartesianCoordinates) {
  auto cell =
      correlation::readers::VaspReader::read(data_dir_ + "Si_cartesian.poscar");

  EXPECT_EQ(cell.atomCount(), 8);
  EXPECT_EQ(cell.elements().size(), 1);
  EXPECT_EQ(cell.elements()[0].symbol, "Si");

  // First atom should be at origin
  auto pos = cell.atoms()[0].position();
  EXPECT_NEAR(pos[0], 0.0, 1e-6);
  EXPECT_NEAR(pos[1], 0.0, 1e-6);
  EXPECT_NEAR(pos[2], 0.0, 1e-6);
}

TEST_F(_23_VaspReader_Tests, ParseSelectiveDynamics) {
  auto cell =
      correlation::readers::VaspReader::read(data_dir_ + "Si_seldyn.poscar");

  // Should have 4 atoms despite the Selective Dynamics line
  EXPECT_EQ(cell.atomCount(), 4);
  EXPECT_EQ(cell.elements().size(), 1);
  EXPECT_EQ(cell.elements()[0].symbol, "Si");
}

TEST_F(_23_VaspReader_Tests, ParseMultiSpecies) {
  auto cell =
      correlation::readers::VaspReader::read(data_dir_ + "SiO_multi.poscar");

  // 4 Si + 4 O = 8 total atoms
  EXPECT_EQ(cell.atomCount(), 8);

  // 2 element types
  EXPECT_EQ(cell.elements().size(), 2);

  // Verify species assignment
  std::map<std::string, int> counts;
  for (const auto &atom : cell.atoms()) {
    counts[cell.elements()[atom.element_id()].symbol]++;
  }
  EXPECT_EQ(counts["Si"], 4);
  EXPECT_EQ(counts["O"], 4);
}

TEST_F(_23_VaspReader_Tests, ScalingFactorApplied) {
  auto cell = correlation::readers::VaspReader::read(data_dir_ + "Si.poscar");

  // The POSCAR has scaling factor 5.43 with unit vectors
  // So lattice should be 5.43 Å
  auto params = cell.lattice_parameters();
  EXPECT_NEAR(params[0], 5.43, 1e-6);
}

TEST_F(_23_VaspReader_Tests, ReaderIsRegisteredInFactory) {
  correlation::readers::VaspReader reader;
  EXPECT_EQ(reader.getName(), "VASP POSCAR");
  EXPECT_FALSE(reader.isTrajectory());

  auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 3);
  EXPECT_EQ(exts[0], "poscar");
  EXPECT_EQ(exts[1], "contcar");
  EXPECT_EQ(exts[2], "vasp");
}

TEST_F(_23_VaspReader_Tests, ReadTrajectoryThrows) {
  correlation::readers::VaspReader reader;
  EXPECT_THROW(reader.readTrajectory(data_dir_ + "Si.poscar"),
               std::runtime_error);
}

TEST_F(_23_VaspReader_Tests, NonExistentFileThrows) {
  EXPECT_THROW(
      correlation::readers::VaspReader::read("nonexistent_file.poscar"),
      std::runtime_error);
}
