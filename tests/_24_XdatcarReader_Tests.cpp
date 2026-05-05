// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "readers/XdatcarReader.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <stdexcept>

namespace {

std::string getTestDataDir() {
  std::vector<std::string> candidates = {
      "../../tests/data/",
      "../tests/data/",
      "tests/data/",
      "data/",
  };
  for (const auto &dir : candidates) {
    if (std::filesystem::exists(dir + "Si.xdatcar")) {
      return dir;
    }
  }
  return "../../tests/data/";
}

} // namespace

class _24_XdatcarReader_Tests : public ::testing::Test {
protected:
  std::string data_dir_;
  void SetUp() override { data_dir_ = getTestDataDir(); }
};

TEST_F(_24_XdatcarReader_Tests, ParseThreeFrameTrajectory) {
  auto frames =
      correlation::readers::XdatcarReader::read(data_dir_ + "Si.xdatcar");

  EXPECT_EQ(frames.size(), 3);
}

TEST_F(_24_XdatcarReader_Tests, FrameAtomCountConsistent) {
  auto frames =
      correlation::readers::XdatcarReader::read(data_dir_ + "Si.xdatcar");

  for (const auto &frame : frames) {
    EXPECT_EQ(frame.atomCount(), 4);
  }
}

TEST_F(_24_XdatcarReader_Tests, LatticeConsistentAcrossFrames) {
  auto frames =
      correlation::readers::XdatcarReader::read(data_dir_ + "Si.xdatcar");

  // All frames should share the same lattice (5.43 Å cubic)
  for (const auto &frame : frames) {
    auto params = frame.lattice_parameters();
    EXPECT_NEAR(params[0], 5.43, 1e-6);
    EXPECT_NEAR(params[1], 5.43, 1e-6);
    EXPECT_NEAR(params[2], 5.43, 1e-6);
  }
}

TEST_F(_24_XdatcarReader_Tests, SpeciesAreCorrect) {
  auto frames =
      correlation::readers::XdatcarReader::read(data_dir_ + "Si.xdatcar");

  ASSERT_FALSE(frames.empty());
  EXPECT_EQ(frames[0].elements().size(), 1);
  EXPECT_EQ(frames[0].elements()[0].symbol, "Si");
}

TEST_F(_24_XdatcarReader_Tests, PositionsDifferBetweenFrames) {
  auto frames =
      correlation::readers::XdatcarReader::read(data_dir_ + "Si.xdatcar");

  ASSERT_GE(frames.size(), 2);
  // First atom in frame 0 vs frame 1 should differ
  auto pos0 = frames[0].atoms()[0].position();
  auto pos1 = frames[1].atoms()[0].position();

  double dist_sq = (pos0[0] - pos1[0]) * (pos0[0] - pos1[0]) +
                   (pos0[1] - pos1[1]) * (pos0[1] - pos1[1]) +
                   (pos0[2] - pos1[2]) * (pos0[2] - pos1[2]);
  EXPECT_GT(dist_sq, 1e-12);
}

TEST_F(_24_XdatcarReader_Tests, ReadTrajectoryReturnsTrajectory) {
  correlation::readers::XdatcarReader reader;
  auto traj = reader.readTrajectory(data_dir_ + "Si.xdatcar");

  EXPECT_EQ(traj.getFrameCount(), 3);
}

TEST_F(_24_XdatcarReader_Tests, ReadStructureReturnsFirstFrame) {
  correlation::readers::XdatcarReader reader;
  auto cell = reader.readStructure(data_dir_ + "Si.xdatcar");

  EXPECT_EQ(cell.atomCount(), 4);
}

TEST_F(_24_XdatcarReader_Tests, ReaderIsRegisteredInFactory) {
  correlation::readers::XdatcarReader reader;
  EXPECT_EQ(reader.getName(), "VASP XDATCAR");
  EXPECT_TRUE(reader.isTrajectory());
  EXPECT_EQ(reader.getExtensions().size(), 1);
  EXPECT_EQ(reader.getExtensions()[0], "xdatcar");
}

TEST_F(_24_XdatcarReader_Tests, NonExistentFileThrows) {
  EXPECT_THROW(
      correlation::readers::XdatcarReader::read("nonexistent_file.xdatcar"),
      std::runtime_error);
}
