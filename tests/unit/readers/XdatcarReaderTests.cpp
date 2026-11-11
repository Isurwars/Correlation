/**
 * @file XdatcarReaderTests.cpp
 * @brief Unit tests for XdatcarReader VASP MD trajectory parsing.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "../../PrecisionTestHelper.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/XdatcarReader.hpp"

#include <filesystem>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <stdexcept>
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
    if (std::filesystem::exists(dir + "xdatcar/Si.xdatcar")) {
      return dir + "xdatcar/";
    }
  }
  return "../../tests/data/xdatcar/";
}

class XdatcarReaderTests : public ::testing::Test {
protected:
  [[nodiscard]] const std::string &dataDir() const noexcept { return data_dir_; }

  void SetUp() override { data_dir_ = getTestDataDir(); }

private:
  std::string data_dir_;
};

} // namespace

TEST_F(XdatcarReaderTests, ParseThreeFrameTrajectory) {
  correlation::readers::XdatcarReader reader;
  const auto traj = reader.readTrajectory(dataDir() + "Si.xdatcar");

  EXPECT_EQ(traj.getFrameCount(), 3);
}

TEST_F(XdatcarReaderTests, FrameAtomCountConsistent) {
  correlation::readers::XdatcarReader reader;
  const auto traj = reader.readTrajectory(dataDir() + "Si.xdatcar");

  for (size_t idx = 0; idx < traj.getFrameCount(); ++idx) {
    const auto &frame = traj.getFrame(idx);
    EXPECT_EQ(frame.atomCount(), 4);
  }
}

TEST_F(XdatcarReaderTests, LatticeConsistentAcrossFrames) {
  correlation::readers::XdatcarReader reader;
  const auto traj = reader.readTrajectory(dataDir() + "Si.xdatcar");

  // All frames should share the same lattice (5.43 Å cubic)
  for (size_t idx = 0; idx < traj.getFrameCount(); ++idx) {
    const auto &frame = traj.getFrame(idx);
    const auto params = frame.latticeParameters();
    EXPECT_THAT(params[0], IsRealEq(static_cast<real_t>(5.43)));
    EXPECT_THAT(params[1], IsRealEq(static_cast<real_t>(5.43)));
    EXPECT_THAT(params[2], IsRealEq(static_cast<real_t>(5.43)));
  }
}

TEST_F(XdatcarReaderTests, SpeciesAreCorrect) {
  correlation::readers::XdatcarReader reader;
  const auto traj = reader.readTrajectory(dataDir() + "Si.xdatcar");

  ASSERT_GT(traj.getFrameCount(), 0U);
  const auto &frame0 = traj.getFrame(0);
  EXPECT_EQ(frame0.elements().size(), 1);
  EXPECT_EQ(frame0.elements()[0].symbol, "Si");
}

TEST_F(XdatcarReaderTests, PositionsDifferBetweenFrames) {
  correlation::readers::XdatcarReader reader;
  const auto traj = reader.readTrajectory(dataDir() + "Si.xdatcar");

  ASSERT_GE(traj.getFrameCount(), 2U);
  // First atom in frame 0 vs frame 1 should differ
  const auto &frame0 = traj.getFrame(0);
  const auto &frame1 = traj.getFrame(1);
  const auto &pos0 = frame0.atoms()[0].position();
  const auto &pos1 = frame1.atoms()[0].position();

  const real_t dist_sq = (pos0[0] - pos1[0]) * (pos0[0] - pos1[0]) +
                         (pos0[1] - pos1[1]) * (pos0[1] - pos1[1]) +
                         (pos0[2] - pos1[2]) * (pos0[2] - pos1[2]);
  EXPECT_GT(dist_sq, static_cast<real_t>(1e-12));
}

TEST_F(XdatcarReaderTests, ReadTrajectoryReturnsTrajectory) {
  correlation::readers::XdatcarReader reader;
  const auto traj = reader.readTrajectory(dataDir() + "Si.xdatcar");

  EXPECT_EQ(traj.getFrameCount(), 3);
}

TEST_F(XdatcarReaderTests, ReadStructureReturnsFirstFrame) {
  correlation::readers::XdatcarReader reader;
  const auto cell = reader.readStructure(dataDir() + "Si.xdatcar");

  EXPECT_EQ(cell.atomCount(), 4);
}

TEST_F(XdatcarReaderTests, ReaderIsRegisteredInFactory) {
  correlation::readers::XdatcarReader const reader;
  EXPECT_EQ(reader.getName(), "VASP XDATCAR");
  EXPECT_TRUE(reader.isTrajectory());
  const auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 1);
  EXPECT_EQ(exts[0], "xdatcar");
}

TEST_F(XdatcarReaderTests, NonExistentFileThrows) {
  correlation::readers::XdatcarReader reader;
  EXPECT_THROW(reader.readTrajectory("nonexistent_file.xdatcar"), std::runtime_error);
}

} // namespace correlation::testing
