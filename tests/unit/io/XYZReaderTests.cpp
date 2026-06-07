// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/XYZReader.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>

class XYZReaderTests : public ::testing::Test {
protected:
  void SetUp() override {
    // Create standard XYZ
    std::ofstream xyz1("test_standard.xyz");
    xyz1 << "2\n";
    xyz1 << "Standard XYZ Comment\n";
    xyz1 << "C 0.0 0.0 0.0\n";
    xyz1 << "O 1.2 0.0 0.0\n";
    xyz1.close();

    // Create Extended XYZ
    std::ofstream xyz2("test_extended.xyz");
    xyz2 << "2\n";
    xyz2
        << "Properties=species:S:1:pos:R:3:forces:R:3 energy=-15.2 Lattice=\"10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0\"\n";
    xyz2 << "C 0.0 0.0 0.0 0.1 0.2 0.3\n";
    xyz2 << "O 1.2 0.0 0.0 -0.1 -0.2 -0.3\n";
    xyz2.close();

    // Create Extended XYZ with different column layout
    std::ofstream xyz3("test_extended_cols.xyz");
    xyz3 << "2\n";
    xyz3 << "Properties=id:I:1:pos:R:3:type:S:1 energy=-10.0\n";
    xyz3 << "1 0.0 0.0 0.0 C\n";
    xyz3 << "2 1.2 0.0 0.0 O\n";
    xyz3.close();
  }

  void TearDown() override {
    std::filesystem::remove("test_standard.xyz");
    std::filesystem::remove("test_extended.xyz");
    std::filesystem::remove("test_extended_cols.xyz");
  }
};

TEST_F(XYZReaderTests, ReadStandardXYZ) {
  correlation::readers::XYZReader reader;
  auto traj = reader.readTrajectory("test_standard.xyz");
  EXPECT_EQ(traj.getFrameCount(), 1);

  auto cell = traj.getFrame(0);
  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
  EXPECT_EQ(cell.atoms()[1].element().symbol, "O");
}

TEST_F(XYZReaderTests, ReadExtendedXYZ) {
  correlation::readers::XYZReader reader;
  auto traj = reader.readTrajectory("test_extended.xyz");
  auto cell = traj.getFrame(0);

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");

  // Check energy
  EXPECT_DOUBLE_EQ(cell.getEnergy(), -15.2);

  // Check lattice
  auto params = cell.lattice_parameters();
  EXPECT_DOUBLE_EQ(params[0], 10.0);
  EXPECT_DOUBLE_EQ(params[1], 10.0);
  EXPECT_DOUBLE_EQ(params[2], 10.0);
}

TEST_F(XYZReaderTests, ReadExtendedXYZCustomColumns) {
  correlation::readers::XYZReader reader;
  auto traj = reader.readTrajectory("test_extended_cols.xyz");
  auto cell = traj.getFrame(0);

  EXPECT_EQ(cell.atomCount(), 2);
  // Column 4 is species ("type")
  EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
  EXPECT_EQ(cell.atoms()[1].element().symbol, "O");

  // Check energy
  EXPECT_DOUBLE_EQ(cell.getEnergy(), -10.0);

  // Pos should be in cols 1,2,3
  EXPECT_DOUBLE_EQ(cell.atoms()[1].position().x(), 1.2);
}
