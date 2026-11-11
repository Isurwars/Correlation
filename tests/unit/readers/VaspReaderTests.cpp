/**
 * @file VaspReaderTests.cpp
 * @brief Unit tests for VaspReader POSCAR and CONTCAR parsing.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "../../PrecisionTestHelper.hpp"
#include "core/Cell.hpp"
#include "readers/VaspReader.hpp"

#include <filesystem>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace correlation::testing {

using correlation::real_t;

namespace {

// Helper to find the test data directory relative to the executable
[[nodiscard]] std::string getTestDataDir() {
  std::vector<std::string> const candidates = {
      "../../tests/data/", // build/tests -> tests/data
      "../tests/data/",    // build -> tests/data
      "tests/data/",       // project root
      "data/",             // tests/data (if cwd is tests)
  };
  for (const auto &dir : candidates) {
    if (std::filesystem::exists(dir + "vasp/Si.poscar")) {
      return dir;
    }
  }
  return "../../tests/data/";
}

class VaspReaderTests : public ::testing::Test {
protected:
  [[nodiscard]] const std::string &dataDir() const noexcept { return data_dir_; }

  void SetUp() override { data_dir_ = getTestDataDir() + "vasp/"; }

private:
  std::string data_dir_;
};

} // namespace

TEST_F(VaspReaderTests, ParseSiDiamondPoscar) {
  const auto cell = correlation::readers::VaspReader::read(dataDir() + "Si.poscar");

  // 8 Si atoms
  EXPECT_EQ(cell.atomCount(), 8);

  // 1 element type
  EXPECT_EQ(cell.elements().size(), 1);
  EXPECT_EQ(cell.elements()[0].symbol, "Si");

  // Lattice should be 5.43 x 5.43 x 5.43 (scaling factor * unit vectors)
  const auto params = cell.latticeParameters();
  EXPECT_THAT(params[0], IsRealEq(static_cast<real_t>(5.43)));
  EXPECT_THAT(params[1], IsRealEq(static_cast<real_t>(5.43)));
  EXPECT_THAT(params[2], IsRealEq(static_cast<real_t>(5.43)));
  EXPECT_THAT(params[3], IsRealEq(static_cast<real_t>(90.0)));
  EXPECT_THAT(params[4], IsRealEq(static_cast<real_t>(90.0)));
  EXPECT_THAT(params[5], IsRealEq(static_cast<real_t>(90.0)));
}

TEST_F(VaspReaderTests, ParseCartesianCoordinates) {
  const auto cell = correlation::readers::VaspReader::read(dataDir() + "Si_cartesian.poscar");

  EXPECT_EQ(cell.atomCount(), 8);
  EXPECT_EQ(cell.elements().size(), 1);
  EXPECT_EQ(cell.elements()[0].symbol, "Si");

  // First atom should be at origin
  const auto pos = cell.atoms()[0].position();
  EXPECT_THAT(pos[0], IsRealEq(static_cast<real_t>(0.0)));
  EXPECT_THAT(pos[1], IsRealEq(static_cast<real_t>(0.0)));
  EXPECT_THAT(pos[2], IsRealEq(static_cast<real_t>(0.0)));
}

TEST_F(VaspReaderTests, ParseSelectiveDynamics) {
  const auto cell = correlation::readers::VaspReader::read(dataDir() + "Si_seldyn.poscar");

  // Should have 4 atoms despite the Selective Dynamics line
  EXPECT_EQ(cell.atomCount(), 4);
  EXPECT_EQ(cell.elements().size(), 1);
  EXPECT_EQ(cell.elements()[0].symbol, "Si");
}

TEST_F(VaspReaderTests, ParseMultiSpecies) {
  const auto cell = correlation::readers::VaspReader::read(dataDir() + "SiO_multi.poscar");

  // 4 Si + 4 O = 8 total atoms
  EXPECT_EQ(cell.atomCount(), 8);

  // 2 element types
  EXPECT_EQ(cell.elements().size(), 2);

  // Verify species assignment
  std::map<std::string, int> counts;
  for (const auto &atom : cell.atoms()) {
    counts[cell.elements()[atom.elementId()].symbol]++;
  }
  EXPECT_EQ(counts["Si"], 4);
  EXPECT_EQ(counts["O"], 4);
}

TEST_F(VaspReaderTests, ScalingFactorApplied) {
  const auto cell = correlation::readers::VaspReader::read(dataDir() + "Si.poscar");

  // The POSCAR has scaling factor 5.43 with unit vectors
  // So lattice should be 5.43 Å
  const auto params = cell.latticeParameters();
  EXPECT_THAT(params[0], IsRealEq(static_cast<real_t>(5.43)));
}

TEST_F(VaspReaderTests, ReaderIsRegisteredInFactory) {
  correlation::readers::VaspReader const reader;
  EXPECT_EQ(reader.getName(), "VASP POSCAR");
  EXPECT_FALSE(reader.isTrajectory());

  const auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 3);
  EXPECT_EQ(exts[0], "poscar");
  EXPECT_EQ(exts[1], "contcar");
  EXPECT_EQ(exts[2], "vasp");
}

TEST_F(VaspReaderTests, ReadTrajectoryThrows) {
  correlation::readers::VaspReader reader;
  EXPECT_THROW(reader.readTrajectory(dataDir() + "Si.poscar"), std::runtime_error);
}

TEST_F(VaspReaderTests, NonExistentFileThrows) {
  EXPECT_THROW(correlation::readers::VaspReader::read("nonexistent_file.poscar"),
               std::runtime_error);
}

} // namespace correlation::testing
