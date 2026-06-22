#include "math/Constants.hpp"
#include "readers/OutmolReader.hpp"
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>

using namespace correlation::readers;

TEST(OutmolReaderTests, Properties) {
  OutmolReader const reader;
  EXPECT_EQ(reader.getName(), "Outmol");
  EXPECT_TRUE(reader.isTrajectory());
  auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 1);
  EXPECT_EQ(exts[0], "outmol");
}

#include <filesystem>
#include <vector>

namespace {

std::string getTestDataDir() {
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

TEST(OutmolReaderTests, ReadsTrajectoryFormat1) {
  // Format 1: $cell vectors + $coordinates
  std::string const data_dir = getTestDataDir();
  OutmolReader reader;
  auto traj = reader.readTrajectory(data_dir + "clean_fmt1.outmol");
  auto cell = reader.readStructure(data_dir + "clean_fmt1.outmol");

  EXPECT_EQ(traj.getFrameCount(), 1);
  const auto &f = traj.getFrame(0);
  EXPECT_EQ(f.atomCount(), 2);

  // Bohr to angstrom check
  EXPECT_NEAR(f.lattice_parameters()[0], 10.0 * correlation::math::bohr_to_angstrom, 1e-5);
  EXPECT_NEAR(f.atoms()[0].position().x(), 1.0 * correlation::math::bohr_to_angstrom, 1e-5);
  EXPECT_EQ(f.atoms()[0].element().symbol, "C");
  EXPECT_EQ(f.atoms()[1].element().symbol, "H");

  EXPECT_EQ(cell.atomCount(), 2);
}

TEST(OutmolReaderTests, ReadsTrajectoryFormat2) {
  // Format 2: $cell vectors + ATOMIC COORDINATES (au)
  std::string const data_dir = getTestDataDir();
  OutmolReader reader;
  auto traj = reader.readTrajectory(data_dir + "clean_fmt2.outmol");

  EXPECT_EQ(traj.getFrameCount(), 1);
  const auto &f = traj.getFrame(0);
  EXPECT_EQ(f.atomCount(), 2);

  EXPECT_NEAR(f.lattice_parameters()[0], 12.0 * correlation::math::bohr_to_angstrom, 1e-5);
  EXPECT_NEAR(f.atoms()[0].position().x(), 2.0 * correlation::math::bohr_to_angstrom, 1e-5);
  EXPECT_EQ(f.atoms()[0].element().symbol, "Si");
  EXPECT_EQ(f.atoms()[1].element().symbol, "O");
}

TEST(OutmolReaderTests, ThrowsOnInvalidFile) {
  OutmolReader reader;
  EXPECT_THROW(reader.readTrajectory("nonexistent_outmol.outmol"), std::runtime_error);
}
