#include "math/Constants.hpp"
#include "readers/CastepMdReader.hpp"

#include <gtest/gtest.h>

using namespace correlation::readers;

TEST(CastepMdReaderTests, Properties) {
  CastepMdReader reader;
  EXPECT_EQ(reader.getName(), "CASTEP MD");
  EXPECT_TRUE(reader.isTrajectory());
  auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 1);
  EXPECT_EQ(exts[0], "md");
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
    if (std::filesystem::exists(dir + "castep_md/clean.md")) {
      return dir + "castep_md/";
    }
  }
  return "../../tests/data/castep_md/";
}

} // namespace

TEST(CastepMdReaderTests, ReadsTrajectory) {
  std::string const data_dir = getTestDataDir();
  CastepMdReader reader;
  auto traj = reader.readTrajectory(data_dir + "clean.md");
  auto struct_cell = reader.readStructure(data_dir + "clean.md");

  // Trajectory checks
  EXPECT_EQ(traj.getFrameCount(), 2);

  // Frame 1 check
  const auto &file_1 = traj.getFrame(0);
  EXPECT_EQ(file_1.atomCount(), 1);
  EXPECT_THAT(file_1.getEnergy(), correlation::testing::IsRealEq(-31.8206146));
  // Bohr to angstrom conversion checking:
  // lattice a = 100.0 bohr = 100.0 * 0.529177210903...
  EXPECT_NEAR(file_1.lattice_parameters()[0], 100.0 * correlation::math::bohr_to_angstrom, 1e-5);
  EXPECT_NEAR(file_1.atoms()[0].position().x(), 1.0 * correlation::math::bohr_to_angstrom, 1e-5);

  // Frame 2 check (reuses lattice but has different energy and positions)
  const auto &file_2 = traj.getFrame(1);
  EXPECT_EQ(file_2.atomCount(), 1);
  EXPECT_THAT(file_2.getEnergy(), correlation::testing::IsRealEq(-40.0));
  EXPECT_NEAR(file_2.lattice_parameters()[0], 100.0 * correlation::math::bohr_to_angstrom, 1e-5);
  EXPECT_NEAR(file_2.atoms()[0].position().x(), 1.1 * correlation::math::bohr_to_angstrom, 1e-5);

  // readStructure returns the first frame
  EXPECT_EQ(struct_cell.atomCount(), 1);
  EXPECT_THAT(struct_cell.getEnergy(), correlation::testing::IsRealEq(-31.8206146));
}
