#include "readers/CarReader.hpp"
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>

using namespace correlation::readers;

TEST(CarReaderTests, Properties) {
  CarReader const reader;
  EXPECT_EQ(reader.getName(), "Accelrys CAR");
  EXPECT_FALSE(reader.isTrajectory());
  auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 1);
  EXPECT_EQ(exts[0], "car");
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
    if (std::filesystem::exists(dir + "car/clean.car")) {
      return dir + "car/";
    }
  }
  return "../../tests/data/car/";
}

} // namespace

TEST(CarReaderTests, ReadsStructure) {
  std::string const data_dir = getTestDataDir();
  CarReader reader;
  auto cell = reader.readStructure(data_dir + "clean.car");

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[0], 10.5);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[1], 11.5);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[2], 12.5);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
  EXPECT_EQ(cell.atoms()[1].element().symbol, "H");
}

TEST(CarReaderTests, ReadsStructurePbcOff) {
  std::string const data_dir = getTestDataDir();
  CarReader reader;
  auto cell = reader.readStructure(data_dir + "clean_off.car");

  EXPECT_EQ(cell.atomCount(), 1);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[0], 100.0); // PBC=OFF sets 100.0
  EXPECT_EQ(cell.atoms()[0].element().symbol, "O");
}

TEST(CarReaderTests, ThrowsOnReadTrajectory) {
  CarReader reader;
  EXPECT_THROW(reader.readTrajectory("dummy.car"), std::runtime_error);
}
