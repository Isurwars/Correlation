// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/XYZReader.hpp"

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
    if (std::filesystem::exists(dir + "xyz/clean.xyz")) {
      return dir + "xyz/";
    }
  }
  return "../../tests/data/xyz/";
}

class XYZReaderTests : public ::testing::Test {
public:
  std::string data_dir_;

protected:
  void SetUp() override { data_dir_ = getTestDataDir(); }
};
} // namespace

TEST_F(XYZReaderTests, ReadStandardXYZ) {
  correlation::readers::XYZReader reader;
  auto traj = reader.readTrajectory(data_dir_ + "clean.xyz");
  EXPECT_EQ(traj.getFrameCount(), 1);

  auto cell = traj.getFrame(0);
  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
  EXPECT_EQ(cell.atoms()[1].element().symbol, "O");
}

TEST_F(XYZReaderTests, ReadExtendedXYZ) {
  correlation::readers::XYZReader reader;
  auto traj = reader.readTrajectory(data_dir_ + "clean_extended.xyz");
  auto cell = traj.getFrame(0);

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");

  // Check energy
  EXPECT_THAT(cell.getEnergy(), correlation::testing::IsRealEq(-15.2));

  // Check lattice
  auto params = cell.lattice_parameters();
  EXPECT_THAT(params[0], correlation::testing::IsRealEq(10.0));
  EXPECT_THAT(params[1], correlation::testing::IsRealEq(10.0));
  EXPECT_THAT(params[2], correlation::testing::IsRealEq(10.0));
}

TEST_F(XYZReaderTests, ReadExtendedXYZCustomColumns) {
  correlation::readers::XYZReader reader;
  auto traj = reader.readTrajectory(data_dir_ + "clean_extended_cols.xyz");
  auto cell = traj.getFrame(0);

  EXPECT_EQ(cell.atomCount(), 2);
  // Column 4 is species ("type")
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
  EXPECT_EQ(cell.atoms()[1].element().symbol, "O");

  // Check energy
  EXPECT_THAT(cell.getEnergy(), correlation::testing::IsRealEq(-10.0));

  // Pos should be in cols 1,2,3
  EXPECT_THAT(cell.atoms()[1].position().x(), correlation::testing::IsRealEq(1.2));
}
