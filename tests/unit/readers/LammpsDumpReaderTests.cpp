#include "readers/LammpsDumpReader.hpp"

#include <gtest/gtest.h>

using namespace correlation::readers;

TEST(LammpsDumpReaderTests, Properties) {
  const LammpsDumpReader reader;
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
  const auto traj = reader.readTrajectory(data_dir + "clean.dump");

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
  const auto traj = reader.readTrajectory(data_dir + "clean_tri.dump");

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
  const auto traj = reader.readTrajectory(std::vector<std::string>{file10, file1, file2});

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
  const auto traj = reader.readTrajectory((temp_dir / "chunk_*.dump").string());

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
  const auto traj = reader.readTrajectory(std::vector<std::string>{file1, file2});

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

namespace {

std::vector<std::string> collectPlasticDumpFiles() {
  std::filesystem::path plastic_dir;
  for (const auto &candidate :
       {"examples/plastic", "../examples/plastic", "../../examples/plastic"}) {
    if (std::filesystem::exists(std::filesystem::path(candidate) / "dump.compresion.0")) {
      plastic_dir = candidate;
      break;
    }
  }
  std::vector<std::string> files;
  if (plastic_dir.empty()) {
    return files;
  }
  for (const auto &entry : std::filesystem::directory_iterator(plastic_dir)) {
    if (entry.is_regular_file() &&
        entry.path().filename().string().starts_with("dump.compresion.")) {
      files.push_back(entry.path().string());
    }
  }
  return files;
}

void verifyFrameCellMetrics(const correlation::core::Trajectory &traj, size_t index, double exp_a,
                            double exp_b, double exp_c, double exp_vol) {
  EXPECT_NEAR(traj.getLatticeParameters(index)[0], exp_a, 1e-3);
  EXPECT_NEAR(traj.getLatticeParameters(index)[1], exp_b, 1e-3);
  EXPECT_NEAR(traj.getLatticeParameters(index)[2], exp_c, 1e-3);
  EXPECT_NEAR(traj.getVolume(index), exp_vol, 1.0);
}

void verifyPlasticCellMetrics(const correlation::core::Trajectory &traj) {
  EXPECT_EQ(traj.getFrameLatticeParameters().size(), 81U);
  EXPECT_EQ(traj.getFrameVolumes().size(), 81U);
  verifyFrameCellMetrics(traj, 0, 70.6508, 70.6508, 70.6508, 352656.30);
  verifyFrameCellMetrics(traj, 40, 42.3905, 84.8389, 97.5743, 350912.55);
  verifyFrameCellMetrics(traj, 80, 14.1302, 137.6423, 180.1181, 350313.08);
}

void verifyPlasticElementList(const correlation::core::Cell &first_frame,
                              const correlation::core::Cell &last_frame) {
  EXPECT_EQ(first_frame.elements().size(), 4U);
  EXPECT_EQ(last_frame.elements().size(), 4U);
  EXPECT_EQ(first_frame.elements()[0].symbol, last_frame.elements()[0].symbol);
  EXPECT_EQ(first_frame.elements()[1].symbol, last_frame.elements()[1].symbol);
  EXPECT_EQ(first_frame.elements()[2].symbol, last_frame.elements()[2].symbol);
  EXPECT_EQ(first_frame.elements()[3].symbol, last_frame.elements()[3].symbol);
}

void verifyPlasticAtomConsistency(const correlation::core::Cell &first_frame,
                                  const correlation::core::Cell &last_frame) {
  verifyPlasticElementList(first_frame, last_frame);
  EXPECT_EQ(first_frame.atomCount(), 16000U);
  EXPECT_EQ(last_frame.atomCount(), 16000U);
  EXPECT_EQ(first_frame.atoms()[0].element().id.value, last_frame.atoms()[0].element().id.value);
  EXPECT_EQ(first_frame.atoms()[7999].element().id.value,
            last_frame.atoms()[7999].element().id.value);
  EXPECT_EQ(first_frame.atoms()[15999].element().id.value,
            last_frame.atoms()[15999].element().id.value);
}

} // namespace

TEST(LammpsDumpReaderTests, SortsAtomsByIdAndShiftsByBoxOrigin) {
  const std::filesystem::path temp_dir =
      std::filesystem::temp_directory_path() / "lammps_id_sort_origin_test";
  std::filesystem::create_directories(temp_dir);

  const std::string file1 = (temp_dir / "dump.compresion.0").string();
  const std::string file2 = (temp_dir / "dump.compresion.10000").string();

  // Frame 0: unsorted IDs (3, 1, 2) with types (4, 1, 3) and non-zero box origin (-0.5, 10.0, -5.0)
  {
    std::ofstream out(file1);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n3\nITEM: BOX BOUNDS pp pp pp\n"
        << "-0.5 9.5\n10.0 30.0\n-5.0 15.0\nITEM: ATOMS id type x y z\n"
        << "3 4 4.5 15.0 0.0\n"
        << "1 1 0.5 12.0 -4.0\n"
        << "2 3 2.5 14.0 -2.0\n";
  }

  // Frame 1: different unsorted ID order (2, 3, 1) and compressed box bounds
  {
    std::ofstream out(file2);
    out << "ITEM: TIMESTEP\n10000\nITEM: NUMBER OF ATOMS\n3\nITEM: BOX BOUNDS pp pp pp\n"
        << "0.5 8.5\n5.0 35.0\n-10.0 20.0\nITEM: ATOMS id type x y z\n"
        << "2 3 3.5 9.0 -7.0\n"
        << "3 4 5.5 10.0 -5.0\n"
        << "1 1 1.5 7.0 -9.0\n";
  }

  LammpsDumpReader reader;
  const auto traj = reader.readTrajectory(std::vector<std::string>{file2, file1});

  EXPECT_EQ(traj.getFrameCount(), 2U);
  EXPECT_DOUBLE_EQ(traj.getTimeStep(), 1.0);

  const auto &frame_0 = traj.getFrame(0);
  const auto &frame_1 = traj.getFrame(1);
  EXPECT_EQ(frame_0.atomCount(), 3U);
  EXPECT_EQ(frame_1.atomCount(), 3U);

  // Atom ordering and element symbols match across frames (sorted by id: 1->1, 2->3, 3->4)
  EXPECT_EQ(frame_0.atoms()[0].element().symbol, "1");
  EXPECT_EQ(frame_0.atoms()[1].element().symbol, "3");
  EXPECT_EQ(frame_0.atoms()[2].element().symbol, "4");
  EXPECT_EQ(frame_1.atoms()[0].element().symbol, "1");
  EXPECT_EQ(frame_1.atoms()[1].element().symbol, "3");
  EXPECT_EQ(frame_1.atoms()[2].element().symbol, "4");
  EXPECT_EQ(frame_0.atoms()[0].element().id.value, frame_1.atoms()[0].element().id.value);
  EXPECT_EQ(frame_0.atoms()[1].element().id.value, frame_1.atoms()[1].element().id.value);
  EXPECT_EQ(frame_0.atoms()[2].element().id.value, frame_1.atoms()[2].element().id.value);

  // Shifted coordinates: x - xlo, y - ylo, z - zlo
  EXPECT_NEAR(frame_0.atoms()[0].position().x(), 1.0, 1e-6);
  EXPECT_NEAR(frame_0.atoms()[0].position().y(), 2.0, 1e-6);
  EXPECT_NEAR(frame_0.atoms()[0].position().z(), 1.0, 1e-6);

  // Per-cell lattice parameters and volumes tracked in Trajectory
  EXPECT_NEAR(traj.getLatticeParameters(0)[0], 10.0, 1e-6);
  EXPECT_NEAR(traj.getLatticeParameters(0)[1], 20.0, 1e-6);
  EXPECT_NEAR(traj.getLatticeParameters(0)[2], 20.0, 1e-6);
  EXPECT_NEAR(traj.getVolume(0), 4000.0, 1e-4);

  EXPECT_NEAR(traj.getLatticeParameters(1)[0], 8.0, 1e-6);
  EXPECT_NEAR(traj.getLatticeParameters(1)[1], 30.0, 1e-6);
  EXPECT_NEAR(traj.getLatticeParameters(1)[2], 30.0, 1e-6);
  EXPECT_NEAR(traj.getVolume(1), 7200.0, 1e-4);

  std::filesystem::remove_all(temp_dir);
}

TEST(LammpsDumpReaderTests, ReadsPlasticCompressionExampleSeriesAndTracksCellVolumeAndLattice) {
  const std::vector<std::string> files = collectPlasticDumpFiles();
  EXPECT_EQ(files.size(), 81U);

  LammpsDumpReader reader;
  const auto traj = reader.readTrajectory(files);

  EXPECT_EQ(traj.getFrameCount(), 81U);
  EXPECT_DOUBLE_EQ(traj.getTimeStep(), 1.0);
  verifyPlasticCellMetrics(traj);
  verifyPlasticAtomConsistency(traj.getFrame(0), traj.getFrame(80));
}
