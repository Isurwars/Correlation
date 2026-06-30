#include "readers/ArcReader.hpp"

#include <gtest/gtest.h>

using namespace correlation::readers;

TEST(ArcReaderTests, Properties) {
  ArcReader const reader;
  EXPECT_EQ(reader.getName(), "Accelrys ARC");
  EXPECT_TRUE(reader.isTrajectory());
  auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 1);
  EXPECT_EQ(exts[0], "arc");
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
    if (std::filesystem::exists(dir + "arc/clean.arc")) {
      return dir + "arc/";
    }
  }
  return "../../tests/data/arc/";
}

} // namespace

TEST(ArcReaderTests, ReadsTrajectory) {
  std::string const data_dir = getTestDataDir();
  ArcReader reader;
  auto traj = reader.readTrajectory(data_dir + "clean.arc");

  EXPECT_EQ(traj.getFrameCount(), 2);

  // Frame 1 check
  const auto &frame_1 = traj.getFrame(0);
  EXPECT_EQ(frame_1.atomCount(), 1);
  EXPECT_DOUBLE_EQ(frame_1.lattice_parameters()[0], 10.0);
  EXPECT_DOUBLE_EQ(frame_1.getEnergy(), -12.345);
  EXPECT_EQ(frame_1.atoms()[0].element().symbol, "C");

  // Frame 2 check
  const auto &frame_2 = traj.getFrame(1);
  EXPECT_EQ(frame_2.atomCount(), 1);
  EXPECT_DOUBLE_EQ(frame_2.lattice_parameters()[0], 100.0); // PBC=OFF sets 100.0
  EXPECT_DOUBLE_EQ(frame_2.getEnergy(), 42.0);
  EXPECT_EQ(frame_2.atoms()[0].element().symbol, "C");
}

TEST(ArcReaderTests, ThrowsOnReadStructure) {
  ArcReader reader;
  EXPECT_THROW(reader.readStructure("dummy.arc"), std::runtime_error);
}
