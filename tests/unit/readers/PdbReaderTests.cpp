/**
 * @file PdbReaderTests.cpp
 * @brief Unit tests for PdbReader structure and trajectory parsing.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "../../PrecisionTestHelper.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/PdbReader.hpp"

#include <filesystem>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace correlation::testing {

using correlation::real_t;

namespace {

[[nodiscard]] std::string getTestDataDir() {
  std::vector<std::string> const candidates = {
      "../../tests/data/",
      "../tests/data/",
      "tests/data/",
      "data/",
  };
  for (const auto &dir : candidates) {
    if (std::filesystem::exists(dir + "pdb/clean.pdb")) {
      return dir + "pdb/";
    }
  }
  return "../../tests/data/pdb/";
}

class PdbReaderTests : public ::testing::Test {
protected:
  [[nodiscard]] const std::string &dataDir() const noexcept { return data_dir_; }

  void SetUp() override { data_dir_ = getTestDataDir(); }

private:
  std::string data_dir_;
};

} // namespace

TEST_F(PdbReaderTests, ReadSingleStructure) {
  correlation::readers::PdbReader reader;
  const auto cell = reader.readStructure(dataDir() + "clean.pdb");

  // Check lattice
  const auto params = cell.latticeParameters();
  EXPECT_THAT(params[0], IsRealEq(static_cast<real_t>(20.0)));
  EXPECT_THAT(params[1], IsRealEq(static_cast<real_t>(25.0)));
  EXPECT_THAT(params[2], IsRealEq(static_cast<real_t>(30.0)));

  // Check atoms
  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "N");
  EXPECT_THAT(cell.atoms()[0].position().x(), IsRealEq(static_cast<real_t>(1.0)));
  EXPECT_THAT(cell.atoms()[0].position().y(), IsRealEq(static_cast<real_t>(2.0)));
  EXPECT_THAT(cell.atoms()[0].position().z(), IsRealEq(static_cast<real_t>(3.0)));

  EXPECT_EQ(cell.atoms()[1].element().symbol, "C");
  EXPECT_THAT(cell.atoms()[1].position().x(), IsRealEq(static_cast<real_t>(2.0)));
}

TEST_F(PdbReaderTests, ReadFallbackElements) {
  correlation::readers::PdbReader reader;
  const auto cell = reader.readStructure(dataDir() + "clean_fallback.pdb");

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "O");
  EXPECT_EQ(cell.atoms()[1].element().symbol, "H");
}

TEST_F(PdbReaderTests, ReadTrajectoryFrames) {
  correlation::readers::PdbReader reader;
  const auto traj = reader.readTrajectory(dataDir() + "clean_traj.pdb");

  EXPECT_EQ(traj.getFrameCount(), 2);

  // Frame 1
  const auto &cell1 = traj.getFrame(0);
  EXPECT_EQ(cell1.atomCount(), 1);
  EXPECT_THAT(cell1.atoms()[0].position().x(), IsRealEq(static_cast<real_t>(0.0)));

  // Frame 2
  const auto &cell2 = traj.getFrame(1);
  EXPECT_EQ(cell2.atomCount(), 1);
  EXPECT_THAT(cell2.atoms()[0].position().x(), IsRealEq(static_cast<real_t>(1.0)));
}

TEST_F(PdbReaderTests, ReadTrajectorySingleFrameFallback) {
  correlation::readers::PdbReader reader;
  const auto traj = reader.readTrajectory(dataDir() + "clean.pdb");

  EXPECT_EQ(traj.getFrameCount(), 1);
  EXPECT_EQ(traj.getFrame(0).atomCount(), 2);
}

TEST_F(PdbReaderTests, ThrowOnNonExistentFile) {
  correlation::readers::PdbReader reader;
  EXPECT_THROW(reader.readStructure("non_existent_file.pdb"), std::runtime_error);
  EXPECT_THROW(reader.readTrajectory("non_existent_file.pdb"), std::runtime_error);
}

} // namespace correlation::testing
