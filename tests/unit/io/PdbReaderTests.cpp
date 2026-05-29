// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/PdbReader.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>

namespace correlation::testing {

class PdbReaderTests : public ::testing::Test {
protected:
  std::string makePdbAtomLine(const std::string &name, double x, double y, double z, const std::string &element) {
    std::string line(80, ' ');
    line.replace(0, 4, "ATOM");
    line.replace(12, std::min<size_t>(4, name.length()), name);
    char buf[10];
    std::snprintf(buf, sizeof(buf), "%8.3f", x);
    line.replace(30, 8, buf);
    std::snprintf(buf, sizeof(buf), "%8.3f", y);
    line.replace(38, 8, buf);
    std::snprintf(buf, sizeof(buf), "%8.3f", z);
    line.replace(46, 8, buf);
    line.replace(76, std::min<size_t>(2, element.length()), element);
    return line;
  }

  std::string makePdbCryst1Line(double a, double b, double c) {
    std::string line(80, ' ');
    line.replace(0, 6, "CRYST1");
    char buf[10];
    std::snprintf(buf, sizeof(buf), "%9.3f", a);
    line.replace(6, 9, buf);
    std::snprintf(buf, sizeof(buf), "%9.3f", b);
    line.replace(15, 9, buf);
    std::snprintf(buf, sizeof(buf), "%9.3f", c);
    line.replace(24, 9, buf);
    return line;
  }

  void SetUp() override {
    // Create standard single structure PDB
    std::ofstream pdb1("test_single.pdb");
    pdb1 << makePdbCryst1Line(20.0, 25.0, 30.0) << "\n";
    pdb1 << makePdbAtomLine("N", 1.0, 2.0, 3.0, "N") << "\n";
    pdb1 << makePdbAtomLine("CA", 2.0, 3.0, 4.0, "C") << "\n";
    pdb1 << "END\n";
    pdb1.close();

    // Create PDB with fallback element parsing (no element at 76-78)
    std::ofstream pdb2("test_fallback.pdb");
    pdb2 << makePdbCryst1Line(10.0, 10.0, 10.0) << "\n";
    pdb2 << makePdbAtomLine("O1", 0.0, 0.0, 0.0, "") << "\n"; // should fallback to "O"
    pdb2 << makePdbAtomLine("1H", 1.0, 1.0, 1.0, "") << "\n"; // should fallback to "H" via PDB name fixes
    pdb2 << "END\n";
    pdb2.close();

    // Create trajectory PDB
    std::ofstream pdb3("test_traj.pdb");
    pdb3 << makePdbCryst1Line(15.0, 15.0, 15.0) << "\n";
    pdb3 << "MODEL        1\n";
    pdb3 << makePdbAtomLine("H", 0.0, 0.0, 0.0, "H") << "\n";
    pdb3 << "ENDMDL\n";
    pdb3 << "MODEL        2\n";
    pdb3 << makePdbAtomLine("H", 1.0, 1.0, 1.0, "H") << "\n";
    pdb3 << "ENDMDL\n";
    pdb3.close();
  }

  void TearDown() override {
    std::filesystem::remove("test_single.pdb");
    std::filesystem::remove("test_fallback.pdb");
    std::filesystem::remove("test_traj.pdb");
  }
};

TEST_F(PdbReaderTests, ReadSingleStructure) {
  correlation::readers::PdbReader reader;
  auto cell = reader.readStructure("test_single.pdb");

  // Check lattice
  auto params = cell.lattice_parameters();
  EXPECT_DOUBLE_EQ(params[0], 20.0);
  EXPECT_DOUBLE_EQ(params[1], 25.0);
  EXPECT_DOUBLE_EQ(params[2], 30.0);

  // Check atoms
  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "N");
  EXPECT_DOUBLE_EQ(cell.atoms()[0].position().x(), 1.0);
  EXPECT_DOUBLE_EQ(cell.atoms()[0].position().y(), 2.0);
  EXPECT_DOUBLE_EQ(cell.atoms()[0].position().z(), 3.0);

  EXPECT_EQ(cell.atoms()[1].element().symbol, "C");
  EXPECT_DOUBLE_EQ(cell.atoms()[1].position().x(), 2.0);
}

TEST_F(PdbReaderTests, ReadFallbackElements) {
  correlation::readers::PdbReader reader;
  auto cell = reader.readStructure("test_fallback.pdb");

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "O");
  EXPECT_EQ(cell.atoms()[1].element().symbol, "H");
}

TEST_F(PdbReaderTests, ReadTrajectoryFrames) {
  correlation::readers::PdbReader reader;
  auto traj = reader.readTrajectory("test_traj.pdb");

  EXPECT_EQ(traj.getFrameCount(), 2);

  // Frame 1
  auto cell1 = traj.getFrame(0);
  EXPECT_EQ(cell1.atomCount(), 1);
  EXPECT_DOUBLE_EQ(cell1.atoms()[0].position().x(), 0.0);

  // Frame 2
  auto cell2 = traj.getFrame(1);
  EXPECT_EQ(cell2.atomCount(), 1);
  EXPECT_DOUBLE_EQ(cell2.atoms()[0].position().x(), 1.0);
}

TEST_F(PdbReaderTests, ReadTrajectorySingleFrameFallback) {
  correlation::readers::PdbReader reader;
  auto traj = reader.readTrajectory("test_single.pdb");

  EXPECT_EQ(traj.getFrameCount(), 1);
  EXPECT_EQ(traj.getFrame(0).atomCount(), 2);
}

TEST_F(PdbReaderTests, ThrowOnNonExistentFile) {
  correlation::readers::PdbReader reader;
  EXPECT_THROW(reader.readStructure("non_existent_file.pdb"), std::runtime_error);
  EXPECT_THROW(reader.readTrajectory("non_existent_file.pdb"), std::runtime_error);
}

} // namespace correlation::testing
