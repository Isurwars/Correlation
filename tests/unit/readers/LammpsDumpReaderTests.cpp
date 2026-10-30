#include "readers/LammpsDumpReader.hpp"

#include <gtest/gtest.h>

using namespace correlation::readers;

TEST(LammpsDumpReaderTests, Properties) {
  LammpsDumpReader reader;
  EXPECT_EQ(reader.getName(), "LAMMPS Dump");
  EXPECT_TRUE(reader.isTrajectory());
  auto exts = reader.getExtensions();
  EXPECT_EQ(exts.size(), 2);
  EXPECT_EQ(exts[0], "dump");
  EXPECT_EQ(exts[1], "lammpstrj");
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
    if (std::filesystem::exists(dir + "lammps/clean.dump")) {
      return dir + "lammps/";
    }
  }
  return "../../tests/data/lammps/";
}

} // namespace

TEST(LammpsDumpReaderTests, ReadsTrajectoryOrtho) {
  std::string const data_dir = getTestDataDir();
  LammpsDumpReader reader;
  auto traj = reader.readTrajectory(data_dir + "clean.dump");

  EXPECT_EQ(traj.getFrameCount(), 1);
  const auto &frame = traj.getFrame(0);
  EXPECT_EQ(frame.atomCount(), 1);
  EXPECT_DOUBLE_EQ(frame.latticeParameters()[0], 10.0);
  EXPECT_DOUBLE_EQ(frame.latticeParameters()[1], 11.0);
  EXPECT_DOUBLE_EQ(frame.latticeParameters()[2], 12.0);
  // Standard atom mapping should fall back to type since no element was specified
  EXPECT_EQ(frame.atoms()[0].element().symbol, "2");
}

TEST(LammpsDumpReaderTests, ReadsTrajectoryTriclinicScaledAndElement) {
  std::string const data_dir = getTestDataDir();
  LammpsDumpReader reader;
  auto traj = reader.readTrajectory(data_dir + "clean_tri.dump");

  EXPECT_EQ(traj.getFrameCount(), 1);
  const auto &frame = traj.getFrame(0);
  EXPECT_EQ(frame.atomCount(), 2);

  // Check elements
  EXPECT_EQ(frame.atoms()[0].element().symbol, "C");
  EXPECT_EQ(frame.atoms()[1].element().symbol, "H");

  // Triclinic lattice:
  // lx = 10, ly = 10, lz = 10
  // a = (lx, 0, 0) = (10, 0, 0)
  // b = (xy, ly, 0) = (1, 10, 0)
  // c = (xz, yz, lz) = (2, 3, 10)
  // For Atom 1 (xs=0.5, ys=0.5, zs=0.5):
  // pos = 0.5*a + 0.5*b + 0.5*c = (5.0, 0.0, 0.0) + (0.5, 5.0, 0.0) + (1.0, 1.5, 5.0) =
  // (6.5, 6.5, 5.0)
  const auto &pos1 = frame.atoms()[0].position();
  EXPECT_NEAR(pos1.x(), 6.5, 1e-5);
  EXPECT_NEAR(pos1.y(), 6.5, 1e-5);
  EXPECT_NEAR(pos1.z(), 5.0, 1e-5);
}

TEST(LammpsDumpReaderTests, ThrowsOnInvalidFile) {
  LammpsDumpReader reader;
  EXPECT_THROW(reader.readTrajectory("nonexistent_file.dump"), std::runtime_error);
}

#include "readers/NaturalSort.hpp"
#include <fstream>

TEST(LammpsDumpReaderTests, NaturalSortOrdering) {
  std::vector<std::string> files = {
      "dump.10.lammpstrj",
      "dump.2.lammpstrj",
      "dump.1.lammpstrj",
      "dump.20.lammpstrj",
  };
  naturalSort(files);

  EXPECT_EQ(files[0], "dump.1.lammpstrj");
  EXPECT_EQ(files[1], "dump.2.lammpstrj");
  EXPECT_EQ(files[2], "dump.10.lammpstrj");
  EXPECT_EQ(files[3], "dump.20.lammpstrj");
}

namespace {

void writeDumpFrame(const std::string &path, int64_t timestep, float x_coord) {
  std::ofstream out(path);
  out << "ITEM: TIMESTEP\n"
      << timestep << "\n"
      << "ITEM: NUMBER OF ATOMS\n"
      << 1 << "\n"
      << "ITEM: BOX BOUNDS pp pp pp\n"
      << "0.0 10.0\n0.0 10.0\n0.0 10.0\n"
      << "ITEM: ATOMS id type x y z\n"
      << "1 1 " << x_coord << " 0.0 0.0\n";
}

} // namespace

TEST(LammpsDumpReaderTests, MultiFileReadAndNaturalOrder) {
  const std::filesystem::path temp_dir =
      std::filesystem::temp_directory_path() / "lammps_multi_test";
  std::filesystem::create_directories(temp_dir);

  const std::string file1 = (temp_dir / "dump.1.lammpstrj").string();
  const std::string file2 = (temp_dir / "dump.2.lammpstrj").string();
  const std::string file10 = (temp_dir / "dump.10.lammpstrj").string();

  writeDumpFrame(file1, 100, 1.0F);
  writeDumpFrame(file2, 200, 2.0F);
  writeDumpFrame(file10, 1000, 10.0F);

  LammpsDumpReader reader;
  // Pass in unsorted order to verify natural sorting
  auto traj = reader.readTrajectory(std::vector<std::string>{file10, file1, file2});

  EXPECT_EQ(traj.getFrameCount(), 3);
  EXPECT_NEAR(traj.getFrame(0).atoms()[0].position().x(), 1.0, 1e-4);
  EXPECT_NEAR(traj.getFrame(1).atoms()[0].position().x(), 2.0, 1e-4);
  EXPECT_NEAR(traj.getFrame(2).atoms()[0].position().x(), 10.0, 1e-4);

  std::filesystem::remove_all(temp_dir);
}

TEST(LammpsDumpReaderTests, WildcardPatternExpansion) {
  const std::filesystem::path temp_dir =
      std::filesystem::temp_directory_path() / "lammps_glob_test";
  std::filesystem::create_directories(temp_dir);

  writeDumpFrame((temp_dir / "chunk_01.dump").string(), 10, 1.0F);
  writeDumpFrame((temp_dir / "chunk_02.dump").string(), 20, 2.0F);

  LammpsDumpReader reader;
  auto traj = reader.readTrajectory((temp_dir / "chunk_*.dump").string());

  EXPECT_EQ(traj.getFrameCount(), 2);
  EXPECT_NEAR(traj.getFrame(0).atoms()[0].position().x(), 1.0, 1e-4);
  EXPECT_NEAR(traj.getFrame(1).atoms()[0].position().x(), 2.0, 1e-4);

  std::filesystem::remove_all(temp_dir);
}

TEST(LammpsDumpReaderTests, DeduplicatesConsecutiveRestartTimesteps) {
  const std::filesystem::path temp_dir =
      std::filesystem::temp_directory_path() / "lammps_restart_test";
  std::filesystem::create_directories(temp_dir);

  const std::string file1 = (temp_dir / "run1.dump").string();
  const std::string file2 = (temp_dir / "run2.dump").string();

  // File 1 has steps 100, 200
  {
    std::ofstream out(file1);
    out << "ITEM: TIMESTEP\n100\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n"
        << "0.0 10.0\n0.0 10.0\n0.0 10.0\nITEM: ATOMS id type x y z\n1 1 1.0 0.0 0.0\n"
        << "ITEM: TIMESTEP\n200\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n"
        << "0.0 10.0\n0.0 10.0\n0.0 10.0\nITEM: ATOMS id type x y z\n1 1 2.0 0.0 0.0\n";
  }

  // File 2 restarts at 200 (duplicate), continues to 300
  {
    std::ofstream out(file2);
    out << "ITEM: TIMESTEP\n200\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n"
        << "0.0 10.0\n0.0 10.0\n0.0 10.0\nITEM: ATOMS id type x y z\n1 1 2.0 0.0 0.0\n"
        << "ITEM: TIMESTEP\n300\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n"
        << "0.0 10.0\n0.0 10.0\n0.0 10.0\nITEM: ATOMS id type x y z\n1 1 3.0 0.0 0.0\n";
  }

  LammpsDumpReader reader;
  auto traj = reader.readTrajectory(std::vector<std::string>{file1, file2});

  // Step 200 in file 2 should be deduplicated
  EXPECT_EQ(traj.getFrameCount(), 3);
  EXPECT_NEAR(traj.getFrame(0).atoms()[0].position().x(), 1.0, 1e-4);
  EXPECT_NEAR(traj.getFrame(1).atoms()[0].position().x(), 2.0, 1e-4);
  EXPECT_NEAR(traj.getFrame(2).atoms()[0].position().x(), 3.0, 1e-4);

  std::filesystem::remove_all(temp_dir);
}

TEST(LammpsDumpReaderTests, ThrowsOnDecreasingTimesteps) {
  const std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "lammps_oor_test";
  std::filesystem::create_directories(temp_dir);

  const std::string file1 = (temp_dir / "part1.dump").string();
  const std::string file2 = (temp_dir / "part2.dump").string();

  writeDumpFrame(file1, 500, 1.0F);
  writeDumpFrame(file2, 100, 2.0F);

  LammpsDumpReader reader;
  EXPECT_THROW(reader.readTrajectory(std::vector<std::string>{file1, file2}), std::runtime_error);

  std::filesystem::remove_all(temp_dir);
}
