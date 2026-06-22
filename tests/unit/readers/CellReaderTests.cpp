#include "readers/CellReader.hpp"
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>

using namespace correlation::readers;

TEST(CellReaderTests, Properties) {
  CellReader reader;
  EXPECT_EQ(reader.getName(), "CASTEP CELL");
  EXPECT_FALSE(reader.isTrajectory());
  auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 1);
  EXPECT_EQ(exts[0], "cell");
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
    if (std::filesystem::exists(dir + "cell/clean.cell")) {
      return dir + "cell/";
    }
  }
  return "../../tests/data/cell/";
}

} // namespace

TEST(CellReaderTests, ReadsStructureLatticeAbc) {
  std::string const data_dir = getTestDataDir();
  CellReader reader;
  auto cell = reader.readStructure(data_dir + "clean.cell");

  EXPECT_EQ(cell.atomCount(), 1);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[0], 15.0);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[5], 120.0);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
}

TEST(CellReaderTests, ReadsStructureLatticeCartAndFrac) {
  std::string const data_dir = getTestDataDir();
  CellReader reader;
  auto cell = reader.readStructure(data_dir + "clean_cart.cell");

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[0], 10.0);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[1], 10.0);
  EXPECT_DOUBLE_EQ(cell.lattice_parameters()[2], 10.0);

  EXPECT_EQ(cell.atoms()[0].element().symbol, "Si");
  EXPECT_EQ(cell.atoms()[1].element().symbol, "O");

  // Wrapped position check: Si should be (5.0, 5.0, 5.0)
  const auto &p1 = cell.atoms()[0].position();
  EXPECT_NEAR(p1.x(), 5.0, 1e-5);
  EXPECT_NEAR(p1.y(), 5.0, 1e-5);
  EXPECT_NEAR(p1.z(), 5.0, 1e-5);

  // Wrapped position check: O (1.2 0.1 -0.3) -> (0.2 0.1 0.7) -> (2.0, 1.0, 7.0)
  const auto &p2 = cell.atoms()[1].position();
  EXPECT_NEAR(p2.x(), 2.0, 1e-5);
  EXPECT_NEAR(p2.y(), 1.0, 1e-5);
  EXPECT_NEAR(p2.z(), 7.0, 1e-5);
}

TEST(CellReaderTests, ThrowsOnReadTrajectory) {
  CellReader reader;
  EXPECT_THROW(reader.readTrajectory("dummy.cell"), std::runtime_error);
}
