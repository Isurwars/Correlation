/**
 * @file XYZReaderTests.cpp
 * @brief Unit tests for XYZReader standard and extended XYZ trajectory parsing.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "../../PrecisionTestHelper.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/XYZReader.hpp"

#include <filesystem>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace correlation::testing {

using correlation::real_t;
using correlation::testing::IsRealEq;

namespace {

[[nodiscard]] std::string getTestDataDir() {
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
protected:
  [[nodiscard]] const std::string &dataDir() const noexcept { return data_dir_; }

  void SetUp() override { data_dir_ = getTestDataDir(); }

private:
  std::string data_dir_;
};

} // namespace

TEST_F(XYZReaderTests, ReadStandardXYZ) {
  correlation::readers::XYZReader reader;
  const auto traj = reader.readTrajectory(dataDir() + "clean.xyz");
  EXPECT_EQ(traj.getFrameCount(), 1);

  const auto cell = traj.getFrame(0);
  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
  EXPECT_EQ(cell.atoms()[1].element().symbol, "O");
}

TEST_F(XYZReaderTests, ReadExtendedXYZ) {
  correlation::readers::XYZReader reader;
  const auto traj = reader.readTrajectory(dataDir() + "clean_extended.xyz");
  const auto cell = traj.getFrame(0);

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");

  // Check energy
  EXPECT_THAT(cell.getEnergy(), IsRealEq(static_cast<real_t>(-15.2)));

  // Check lattice
  const auto params = cell.latticeParameters();
  EXPECT_THAT(params[0], IsRealEq(static_cast<real_t>(10.0)));
  EXPECT_THAT(params[1], IsRealEq(static_cast<real_t>(10.0)));
  EXPECT_THAT(params[2], IsRealEq(static_cast<real_t>(10.0)));
}

TEST_F(XYZReaderTests, ReadExtendedXYZCustomColumns) {
  correlation::readers::XYZReader reader;
  const auto traj = reader.readTrajectory(dataDir() + "clean_extended_cols.xyz");
  const auto cell = traj.getFrame(0);

  EXPECT_EQ(cell.atomCount(), 2);
  // Column 4 is species ("type")
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
  EXPECT_EQ(cell.atoms()[1].element().symbol, "O");

  // Check energy
  EXPECT_THAT(cell.getEnergy(), IsRealEq(static_cast<real_t>(-10.0)));

  // Pos should be in cols 1,2,3
  EXPECT_THAT(cell.atoms()[1].position().x(), IsRealEq(static_cast<real_t>(1.2)));
}

} // namespace correlation::testing
