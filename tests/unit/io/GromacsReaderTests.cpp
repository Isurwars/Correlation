// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "readers/GromacsReader.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>

class GromacsReaderTests : public ::testing::Test {
protected:
  void SetUp() override {
    // Create a temporary multi-frame GROMACS file
    std::ofstream gro("test_multi.gro");
    ASSERT_TRUE(gro.is_open());

    // Frame 1
    gro << "Water molecule\n";
    gro << "    3\n";
    gro << "    1SOL     OW    1   0.100   0.200   0.300  0.1000  0.2000  0.3000\n";
    gro << "    1SOL    HW1    2   0.150   0.250   0.350 -0.1000 -0.2000 -0.3000\n";
    gro << "    1SOL    HW2    3   0.050   0.150   0.250  0.0500  0.0500  0.0500\n";
    gro << "   1.00000   2.00000   3.00000\n";

    // Frame 2
    gro << "Water molecule t= 1.0\n";
    gro << "    3\n";
    gro << "    1SOL     OW    1   0.110   0.210   0.310  0.1000  0.2000  0.3000\n";
    gro << "    1SOL    HW1    2   0.160   0.260   0.360 -0.1000 -0.2000 -0.3000\n";
    gro << "    1SOL    HW2    3   0.060   0.160   0.260  0.0500  0.0500  0.0500\n";
    gro << "   1.10000   2.10000   3.10000\n";
    
    gro.close();
  }

  void TearDown() override {
    std::filesystem::remove("test_multi.gro");
  }
};

TEST_F(GromacsReaderTests, ReadMultiFrameTrajectory) {
  correlation::readers::GromacsReader reader;
  EXPECT_TRUE(reader.isTrajectory());
  
  auto traj = reader.readTrajectory("test_multi.gro");
  EXPECT_EQ(traj.getFrameCount(), 2);

  // Frame 1
  auto frame1 = traj.getFrame(0);
  EXPECT_EQ(frame1.atomCount(), 3);
  EXPECT_EQ(frame1.atoms()[0].element().symbol, "O"); // "OW" -> "O"
  EXPECT_EQ(frame1.atoms()[1].element().symbol, "H"); // "HW1" -> "H"
  
  auto pos1 = frame1.atoms()[0].position();
  // GROMACS is in nm, we convert to Angstroms (* 10)
  EXPECT_DOUBLE_EQ(pos1.x(), 1.0);
  EXPECT_DOUBLE_EQ(pos1.y(), 2.0);
  EXPECT_DOUBLE_EQ(pos1.z(), 3.0);

  auto lat1 = frame1.lattice_parameters();
  EXPECT_DOUBLE_EQ(lat1[0], 10.0);
  EXPECT_DOUBLE_EQ(lat1[1], 20.0);
  EXPECT_DOUBLE_EQ(lat1[2], 30.0);

  // Frame 2
  auto frame2 = traj.getFrame(1);
  EXPECT_EQ(frame2.atomCount(), 3);
  auto pos2 = frame2.atoms()[0].position();
  EXPECT_DOUBLE_EQ(pos2.x(), 1.1);
  EXPECT_DOUBLE_EQ(pos2.y(), 2.1);
  EXPECT_DOUBLE_EQ(pos2.z(), 3.1);

  auto lat2 = frame2.lattice_parameters();
  EXPECT_DOUBLE_EQ(lat2[0], 11.0);
}

TEST_F(GromacsReaderTests, ReadStructureReturnsLastFrame) {
  correlation::readers::GromacsReader reader;
  auto cell = reader.readStructure("test_multi.gro");
  
  EXPECT_EQ(cell.atomCount(), 3);
  auto pos = cell.atoms()[0].position();
  // Should be Frame 2
  EXPECT_DOUBLE_EQ(pos.x(), 1.1);
}
