// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/XdatcarReader.hpp"

#include <filesystem>
#include <gtest/gtest.h>
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
    if (std::filesystem::exists(dir + "xdatcar/Si.xdatcar")) {
      return dir + "xdatcar/";
    }
  }
  return "../../tests/data/xdatcar/";
}

class XdatcarReaderTests : public ::testing::Test {
public:
  std::string data_dir_;

protected:
  void SetUp() override { data_dir_ = getTestDataDir(); }
};
} // namespace

TEST_F(XdatcarReaderTests, ParseThreeFrameTrajectory) {
  correlation::readers::XdatcarReader reader;
  auto traj = reader.readTrajectory(data_dir_ + "Si.xdatcar");

  EXPECT_EQ(traj.getFrameCount(), 3);
}

TEST_F(XdatcarReaderTests, FrameAtomCountConsistent) {
  correlation::readers::XdatcarReader reader;
  auto traj = reader.readTrajectory(data_dir_ + "Si.xdatcar");

  for (size_t i = 0; i < traj.getFrameCount(); ++i) {
    auto frame = traj.getFrame(i);
    EXPECT_EQ(frame.atomCount(), 4);
  }
}

TEST_F(XdatcarReaderTests, LatticeConsistentAcrossFrames) {
  correlation::readers::XdatcarReader reader;
  auto traj = reader.readTrajectory(data_dir_ + "Si.xdatcar");

  // All frames should share the same lattice (5.43 Å cubic)
  for (size_t i = 0; i < traj.getFrameCount(); ++i) {
    auto frame = traj.getFrame(i);
    auto params = frame.lattice_parameters();
    EXPECT_NEAR(params[0], 5.43, 1e-6);
    EXPECT_NEAR(params[1], 5.43, 1e-6);
    EXPECT_NEAR(params[2], 5.43, 1e-6);
  }
}

TEST_F(XdatcarReaderTests, SpeciesAreCorrect) {
  correlation::readers::XdatcarReader reader;
  auto traj = reader.readTrajectory(data_dir_ + "Si.xdatcar");

  ASSERT_GT(traj.getFrameCount(), 0U);
  auto frame0 = traj.getFrame(0);
  EXPECT_EQ(frame0.elements().size(), 1);
  EXPECT_EQ(frame0.elements()[0].symbol, "Si");
}

TEST_F(XdatcarReaderTests, PositionsDifferBetweenFrames) {
  correlation::readers::XdatcarReader reader;
  auto traj = reader.readTrajectory(data_dir_ + "Si.xdatcar");

  ASSERT_GE(traj.getFrameCount(), 2U);
  // First atom in frame 0 vs frame 1 should differ
  auto frame0 = traj.getFrame(0);
  auto frame1 = traj.getFrame(1);
  auto pos0 = frame0.atoms()[0].position();
  auto pos1 = frame1.atoms()[0].position();

  double dist_sq = (pos0[0] - pos1[0]) * (pos0[0] - pos1[0]) + (pos0[1] - pos1[1]) * (pos0[1] - pos1[1]) +
                   (pos0[2] - pos1[2]) * (pos0[2] - pos1[2]);
  EXPECT_GT(dist_sq, 1e-12);
}

TEST_F(XdatcarReaderTests, ReadTrajectoryReturnsTrajectory) {
  correlation::readers::XdatcarReader reader;
  auto traj = reader.readTrajectory(data_dir_ + "Si.xdatcar");

  EXPECT_EQ(traj.getFrameCount(), 3);
}

TEST_F(XdatcarReaderTests, ReadStructureReturnsFirstFrame) {
  correlation::readers::XdatcarReader reader;
  auto cell = reader.readStructure(data_dir_ + "Si.xdatcar");

  EXPECT_EQ(cell.atomCount(), 4);
}

TEST_F(XdatcarReaderTests, ReaderIsRegisteredInFactory) {
  correlation::readers::XdatcarReader reader;
  EXPECT_EQ(reader.getName(), "VASP XDATCAR");
  EXPECT_TRUE(reader.isTrajectory());
  EXPECT_EQ(reader.getExtensions().size(), 1);
  EXPECT_EQ(reader.getExtensions()[0], "xdatcar");
}

TEST_F(XdatcarReaderTests, NonExistentFileThrows) {
  correlation::readers::XdatcarReader reader;
  EXPECT_THROW(reader.readTrajectory("nonexistent_file.xdatcar"), std::runtime_error);
}
