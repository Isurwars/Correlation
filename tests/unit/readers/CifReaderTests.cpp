#include "readers/CifReader.hpp"

#include <gtest/gtest.h>

using namespace correlation::readers;

TEST(CifReaderTests, Properties) {
  CifReader const reader;
  EXPECT_EQ(reader.getName(), "CIF");
  EXPECT_FALSE(reader.isTrajectory());
  auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 1);
  EXPECT_EQ(exts[0], "cif");
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
    if (std::filesystem::exists(dir + "cif/clean.cif")) {
      return dir + "cif/";
    }
  }
  return "../../tests/data/cif/";
}

} // namespace

TEST(CifReaderTests, ReadsStructureWithSymmetryAndUncertainty) {
  std::string const data_dir = getTestDataDir();
  CifReader reader;
  auto cell = reader.readStructure(data_dir + "clean.cif");

  EXPECT_THAT(cell.lattice_parameters()[0], correlation::testing::IsRealEq(5.64));
  EXPECT_THAT(cell.lattice_parameters()[3], correlation::testing::IsRealEq(90.0));

  // Na 0.0 0.0 0.0 with symmetry:
  // (x, y, z) -> 0.0 0.0 0.0
  // (-x, -y, -z) -> 0.0 0.0 0.0 (duplicate, skipped)
  // (x+1/2, y, z) -> 0.5 0.0 0.0
  // Cl 0.25 0.25 0.25 with symmetry:
  // (x, y, z) -> 0.25 0.25 0.25
  // (-x, -y, -z) -> -0.25 -0.25 -0.25 -> 0.75 0.75 0.75
  // (x+1/2, y, z) -> 0.75 0.25 0.25
  // Duplicate Na input (0 0 0) should be skipped entirely.
  // Total atoms: 2 Na (0 0 0 and 0.5 0 0) + 3 Cl = 5 atoms.
  EXPECT_EQ(cell.atomCount(), 5);
}

TEST(CifReaderTests, ThrowsOnReadTrajectory) {
  CifReader reader;
  EXPECT_THROW(reader.readTrajectory("dummy.cif"), std::runtime_error);
}
