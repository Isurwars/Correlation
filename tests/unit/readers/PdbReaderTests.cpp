// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/PdbReader.hpp"

#include <filesystem>
#include <gtest/gtest.h>

namespace correlation::testing {
namespace {
std::string getTestDataDir() {
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
public:
  std::string data_dir_;

protected:
  void SetUp() override { data_dir_ = getTestDataDir(); }
};
} // namespace

TEST_F(PdbReaderTests, ReadSingleStructure) {
  correlation::readers::PdbReader reader;
  auto cell = reader.readStructure(data_dir_ + "clean.pdb");

  // Check lattice
  auto params = cell.lattice_parameters();
  EXPECT_DOUBLE_EQ(params[0], 20.0);
  EXPECT_DOUBLE_EQ(params[1], 25.0);
  EXPECT_DOUBLE_EQ(params[2], 30.0);

  // Check atoms
  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "N");
  EXPECT_DOUBLE_EQ(cell.atoms()[0].position().x(), 1.0);
  EXPECT_DOUBLE_EQ(cell.atoms()[0].position().y(), 2.0);
  EXPECT_DOUBLE_EQ(cell.atoms()[0].position().z(), 3.0);

  EXPECT_EQ(cell.atoms()[1].element().symbol, "C");
  EXPECT_DOUBLE_EQ(cell.atoms()[1].position().x(), 2.0);
}

TEST_F(PdbReaderTests, ReadFallbackElements) {
  correlation::readers::PdbReader reader;
  auto cell = reader.readStructure(data_dir_ + "clean_fallback.pdb");

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "O");
  EXPECT_EQ(cell.atoms()[1].element().symbol, "H");
}

TEST_F(PdbReaderTests, ReadTrajectoryFrames) {
  correlation::readers::PdbReader reader;
  auto traj = reader.readTrajectory(data_dir_ + "clean_traj.pdb");

  EXPECT_EQ(traj.getFrameCount(), 2);

  // Frame 1
  auto cell1 = traj.getFrame(0);
  EXPECT_EQ(cell1.atomCount(), 1);
  EXPECT_DOUBLE_EQ(cell1.atoms()[0].position().x(), 0.0);

  // Frame 2
  auto cell2 = traj.getFrame(1);
  EXPECT_EQ(cell2.atomCount(), 1);
  EXPECT_DOUBLE_EQ(cell2.atoms()[0].position().x(), 1.0);
}

TEST_F(PdbReaderTests, ReadTrajectorySingleFrameFallback) {
  correlation::readers::PdbReader reader;
  auto traj = reader.readTrajectory(data_dir_ + "clean.pdb");

  EXPECT_EQ(traj.getFrameCount(), 1);
  EXPECT_EQ(traj.getFrame(0).atomCount(), 2);
}

TEST_F(PdbReaderTests, ThrowOnNonExistentFile) {
  correlation::readers::PdbReader reader;
  EXPECT_THROW(reader.readStructure("non_existent_file.pdb"), std::runtime_error);
  EXPECT_THROW(reader.readTrajectory("non_existent_file.pdb"), std::runtime_error);
}

} // namespace correlation::testing
