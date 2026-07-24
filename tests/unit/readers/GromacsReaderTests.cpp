// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/GromacsReader.hpp"

#include <filesystem>
#include <gtest/gtest.h>

namespace {

std::string getTestDataDir() {
  std::vector<std::string> const candidates = {
      "../../tests/data/",
      "../tests/data/",
      "tests/data/",
      "data/",
  };
  for (const auto &dir : candidates) {
    if (std::filesystem::exists(dir + "gromacs/clean.gro")) {
      return dir + "gromacs/";
    }
  }
  return "../../tests/data/gromacs/";
}

class GromacsReaderTests : public ::testing::Test {
public:
  std::string data_dir_;

protected:
  void SetUp() override { data_dir_ = getTestDataDir(); }
};
} // namespace

TEST_F(GromacsReaderTests, ReadMultiFrameTrajectory) {
  correlation::readers::GromacsReader reader;
  EXPECT_TRUE(reader.isTrajectory());

  auto traj = reader.readTrajectory(data_dir_ + "clean.gro");
  EXPECT_EQ(traj.getFrameCount(), 2);

  // Frame 1
  auto frame1 = traj.getFrame(0);
  EXPECT_EQ(frame1.atomCount(), 3);
  EXPECT_EQ(frame1.atoms()[0].element().symbol, "O"); // "OW" -> "O"
  EXPECT_EQ(frame1.atoms()[1].element().symbol, "H"); // "HW1" -> "H"

  auto pos1 = frame1.atoms()[0].position();
  // GROMACS is in nm, we convert to Angstroms (* 10)
  EXPECT_THAT(pos1.x(), correlation::testing::IsRealEq(1.0));
  EXPECT_THAT(pos1.y(), correlation::testing::IsRealEq(2.0));
  EXPECT_THAT(pos1.z(), correlation::testing::IsRealEq(3.0));

  auto lat1 = frame1.lattice_parameters();
  EXPECT_THAT(lat1[0], correlation::testing::IsRealEq(10.0));
  EXPECT_THAT(lat1[1], correlation::testing::IsRealEq(20.0));
  EXPECT_THAT(lat1[2], correlation::testing::IsRealEq(30.0));

  // Frame 2
  auto frame2 = traj.getFrame(1);
  EXPECT_EQ(frame2.atomCount(), 3);
  auto pos2 = frame2.atoms()[0].position();
  EXPECT_THAT(pos2.x(), correlation::testing::IsRealEq(1.1));
  EXPECT_THAT(pos2.y(), correlation::testing::IsRealEq(2.1));
  EXPECT_THAT(pos2.z(), correlation::testing::IsRealEq(3.1));

  auto lat2 = frame2.lattice_parameters();
  EXPECT_THAT(lat2[0], correlation::testing::IsRealEq(11.0));
}

TEST_F(GromacsReaderTests, ReadStructureReturnsLastFrame) {
  correlation::readers::GromacsReader reader;
  auto cell = reader.readStructure(data_dir_ + "clean.gro");

  EXPECT_EQ(cell.atomCount(), 3);
  auto pos = cell.atoms()[0].position();
  // Should be Frame 2
  EXPECT_THAT(pos.x(), correlation::testing::IsRealEq(1.1));
}
