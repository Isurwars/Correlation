#include "readers/OnetepDatReader.hpp"
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>

using namespace correlation::readers;

TEST(OnetepDatReaderTests, Properties) {
  OnetepDatReader const reader;
  EXPECT_EQ(reader.getName(), "ONETEP DAT");
  EXPECT_FALSE(reader.isTrajectory());
  auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 1);
  EXPECT_EQ(exts[0], "dat");
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
    if (std::filesystem::exists(dir + "onetep_dat/clean.dat")) {
      return dir + "onetep_dat/";
    }
  }
  return "../../tests/data/onetep_dat/";
}

} // namespace

TEST(OnetepDatReaderTests, ReadsStructureCartesianAndBohr) {
  std::string const data_dir = getTestDataDir();
  OnetepDatReader reader;
  auto cell = reader.readStructure(data_dir + "clean.dat");

  EXPECT_EQ(cell.atomCount(), 1);
  // Since unit bohr was used for lattice: 10.0 bohr (meaning the cell has those parameters directly, wait, does
  // OnetepDatReader convert bohr to angstrom? No, OnetepDatReader doesn't do a bohr-to-angstrom conversion for
  // lattice/positions, it just skips the token "bohr" or "angstrom"). Let's verify. Yes: in OnetepDatReader.cpp: if
  // (lower_first == "angstrom" || lower_first == "bohr") { continue; } So it skips them. Thus cell parameters
  // remain 10.0.
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[0], 10.0);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
  EXPECT_DOUBLE_EQ(cell.atoms()[0].position().x(), 1.0);
}

TEST(OnetepDatReaderTests, ReadsStructureLatticeAbcAndPositionsFrac) {
  std::string const data_dir = getTestDataDir();
  OnetepDatReader reader;
  auto cell = reader.readStructure(data_dir + "clean_abc.dat");

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[0], 15.0);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "H");
  EXPECT_EQ(cell.atoms()[1].element().symbol, "O");

  const auto &p1 = cell.atoms()[0].position();
  EXPECT_NEAR(p1.x(), 1.5, 1e-5);
  EXPECT_NEAR(p1.y(), 1.5, 1e-5);
  EXPECT_NEAR(p1.z(), 1.5, 1e-5);

  const auto &p2 = cell.atoms()[1].position();
  EXPECT_NEAR(p2.x(), 3.0, 1e-5);
  EXPECT_NEAR(p2.y(), 3.0, 1e-5);
  EXPECT_NEAR(p2.z(), 13.5, 1e-5);
}

TEST(OnetepDatReaderTests, ThrowsOnReadTrajectory) {
  OnetepDatReader reader;
  EXPECT_THROW(reader.readTrajectory("dummy.dat"), std::runtime_error);
}
