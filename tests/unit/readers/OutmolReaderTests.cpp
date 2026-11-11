/**
 * @file OutmolReaderTests.cpp
 * @brief Unit tests for OutmolReader trajectory and structure parsing.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "math/Constants.hpp"
#include "readers/OutmolReader.hpp"

#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace correlation::readers {

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
    if (std::filesystem::exists(dir + "outmol/clean_fmt1.outmol")) {
      return dir + "outmol/";
    }
  }
  return "../../tests/data/outmol/";
}

} // namespace

TEST(OutmolReaderTests, Properties) {
  OutmolReader const reader;
  EXPECT_EQ(reader.getName(), "Outmol");
  EXPECT_TRUE(reader.isTrajectory());
  const auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 1);
  EXPECT_EQ(exts[0], "outmol");
}

TEST(OutmolReaderTests, ReadsTrajectoryFormat1) {
  // Format 1: $cell vectors + $coordinates
  std::string const data_dir = getTestDataDir();
  OutmolReader reader;
  const auto traj = reader.readTrajectory(data_dir + "clean_fmt1.outmol");
  const auto cell = reader.readStructure(data_dir + "clean_fmt1.outmol");

  EXPECT_EQ(traj.getFrameCount(), 1);
  const auto &frame = traj.getFrame(0);
  EXPECT_EQ(frame.atomCount(), 2);

  // Bohr to angstrom check
  EXPECT_NEAR(frame.latticeParameters()[0],
              static_cast<real_t>(10.0) * correlation::math::BOHR_TO_ANGSTROM,
              static_cast<real_t>(1e-5));
  EXPECT_NEAR(frame.atoms()[0].position().x(),
              static_cast<real_t>(1.0) * correlation::math::BOHR_TO_ANGSTROM,
              static_cast<real_t>(1e-5));
  EXPECT_EQ(frame.atoms()[0].element().symbol, "C");
  EXPECT_EQ(frame.atoms()[1].element().symbol, "H");

  EXPECT_EQ(cell.atomCount(), 2);
}

TEST(OutmolReaderTests, ReadsTrajectoryFormat2) {
  // Format 2: $cell vectors + ATOMIC COORDINATES (au)
  std::string const data_dir = getTestDataDir();
  OutmolReader reader;
  const auto traj = reader.readTrajectory(data_dir + "clean_fmt2.outmol");

  EXPECT_EQ(traj.getFrameCount(), 1);
  const auto &frame = traj.getFrame(0);
  EXPECT_EQ(frame.atomCount(), 2);

  EXPECT_NEAR(frame.latticeParameters()[0],
              static_cast<real_t>(12.0) * correlation::math::BOHR_TO_ANGSTROM,
              static_cast<real_t>(1e-5));
  EXPECT_NEAR(frame.atoms()[0].position().x(),
              static_cast<real_t>(2.0) * correlation::math::BOHR_TO_ANGSTROM,
              static_cast<real_t>(1e-5));
  EXPECT_EQ(frame.atoms()[0].element().symbol, "Si");
  EXPECT_EQ(frame.atoms()[1].element().symbol, "O");
}

TEST(OutmolReaderTests, ThrowsOnInvalidFile) {
  OutmolReader reader;
  EXPECT_THROW(reader.readTrajectory("nonexistent_outmol.outmol"), std::runtime_error);
}

} // namespace correlation::readers
