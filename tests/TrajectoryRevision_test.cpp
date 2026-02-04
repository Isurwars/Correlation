// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <chrono>

#include "../include/Cell.hpp"
#include "../include/Trajectory.hpp"
#include "../include/FileIO.hpp"

// Test fixture for Trajectory tests
class TrajectoryRevisionTest : public ::testing::Test {
protected:
  // Helper to create a dummy frame with a specific position for an atom
  Cell createFrame(double x, double y, double z) {
      // Create a simple cubic lattice
      linalg::Vector3<double> a = {10.0, 0.0, 0.0};
      linalg::Vector3<double> b = {0.0, 10.0, 0.0};
      linalg::Vector3<double> c = {0.0, 0.0, 10.0};
      
      Cell frame(a, b, c);
      
      // Add one atom
      // addAtom automatically handles element registration
      frame.addAtom("H", {x, y, z});
      
      return frame;
  }

  void SetUp() override {
    // Create a temporary directory for test files
    std::filesystem::create_directory("test_data");
  }

  void TearDown() override {
    // Clean up test files
    std::filesystem::remove_all("test_data");
  }
};

TEST_F(TrajectoryRevisionTest, RemoveDuplicatedFrames) {
  // Create a few cells
  Cell c1;
  c1.addAtom("H", {0.0, 0.0, 0.0});
  c1.addAtom("O", {1.0, 0.0, 0.0});

  Cell c2; // Identical to c1
  c2.addAtom("H", {0.0, 0.0, 0.0});
  c2.addAtom("O", {1.0, 0.0, 0.0});

  Cell c3; // Different positions
  c3.addAtom("H", {0.1, 0.0, 0.0});
  c3.addAtom("O", {1.1, 0.0, 0.0});

  Cell c4; // Identical to c3
  c4.addAtom("H", {0.1, 0.0, 0.0});
  c4.addAtom("O", {1.1, 0.0, 0.0});

  Cell c5; // Different atoms (simulation: just partial difference or different position)
  // To satisfy Trajectory validation, it must have same atom count and symbols in order
  c5.addAtom("H", {0.2, 0.0, 0.0});
  c5.addAtom("O", {1.2, 0.0, 0.0});

  std::vector<Cell> frames = {c1, c2, c3, c4, c5};
  Trajectory traj(frames, 1.0);

  // Duplicates are removed in constructor now
  EXPECT_EQ(traj.getFrameCount(), 3);

  // Redundant call should not change anything
  traj.removeDuplicatedFrames();

  EXPECT_EQ(traj.getFrameCount(), 3);
  // Expected frames: c1, c3, c5
  
  const auto& new_frames = traj.getFrames();
  EXPECT_EQ(new_frames[0].atoms()[0].position().x(), 0.0);
  EXPECT_EQ(new_frames[1].atoms()[0].position().x(), 0.1);
  EXPECT_EQ(new_frames[2].atoms()[0].position().x(), 0.2);
}

TEST_F(TrajectoryRevisionTest, ParseEnergyFromArc) {
  // Create a mock .arc file with energy
  std::string filename = "test_data/energy_test.arc";
  std::ofstream out(filename);
  out << "!BIOSYM archive 3\n";
  out << "PBC=ON\n";
  out << "                                      -12345.6789\n";
  out << "!DATE     Feb 13 17:21:41 2025\n";
  out << "PBC   10.0 10.0 10.0 90.0 90.0 90.0\n";
  out << "He     0.0 0.0 0.0 XXXX 1 xx He 0.0\n";
  out << "end\n";
  out << "end\n";
  out.close();

  auto traj = FileIO::readTrajectory(filename, FileIO::FileType::Arc);
  const auto& frames = traj.getFrames();
  ASSERT_EQ(frames.size(), 1);
  
  EXPECT_DOUBLE_EQ(frames[0].getEnergy(), -12345.6789);
}

TEST_F(TrajectoryRevisionTest, ParseMultipleFramesWithEnergy) {
    std::string filename = "test_data/multi_energy.arc";
    std::ofstream out(filename);
    
    // Frame 1
    out << "!BIOSYM archive 3\n";
    out << "PBC=ON\n";
    out << " -100.0\n"; 
    out << "!DATE\n";
    out << "PBC 10.0 10.0 10.0 90.0 90.0 90.0\n";
    out << "He 0.0 0.0 0.0 XXXX 1 xx He 0.0\n";
    out << "end\n";
    out << "end\n";

    // Frame 2
    out << "!BIOSYM archive 3\n";
    out << "PBC=ON\n";
    out << " -200.5\n"; 
    out << "!DATE\n";
    out << "PBC 10.0 10.0 10.0 90.0 90.0 90.0\n";
    out << "He 1.0 1.0 1.0 XXXX 1 xx He 0.0\n";
    out << "end\n";
    out << "end\n";
    out.close();

    auto traj2 = FileIO::readTrajectory(filename, FileIO::FileType::Arc);
    const auto& frames = traj2.getFrames();
    ASSERT_EQ(frames.size(), 2);
    EXPECT_DOUBLE_EQ(frames[0].getEnergy(), -100.0);
    EXPECT_DOUBLE_EQ(frames[1].getEnergy(), -200.5);
}

TEST_F(TrajectoryRevisionTest, RemovesConsecutiveTriplicates) {
    std::vector<Cell> frames;
    
    // Frame 0: (0,0,0)
    frames.push_back(createFrame(0.0, 0.0, 0.0));
    
    // Frame 1: (0,0,0) - Duplicate of Frame 0
    frames.push_back(createFrame(0.0, 0.0, 0.0));
    
    // Frame 2: (1,1,1) - Distinct
    frames.push_back(createFrame(1.0, 1.0, 1.0));
    
    // Frame 3: (1,1,1) - Duplicate of Frame 2
    frames.push_back(createFrame(1.0, 1.0, 1.0));
    
    // Frame 4: (1,1,1) - Triplicate of Frame 2 (Duplicate of Frame 3)
    frames.push_back(createFrame(1.0, 1.0, 1.0));
    
    // Frame 5: (2,2,2) - Distinct
    frames.push_back(createFrame(2.0, 2.0, 2.0));

    // Initialize Trajectory
    Trajectory traj(frames, 1.0);
    
    // We expect frames 0, 2, and 5 to remain.
    // That is 3 frames total.
    EXPECT_EQ(traj.getFrameCount(), 3);
    EXPECT_EQ(traj.getRemovedFrameCount(), 3); // 3 frames removed (1, 3, 4)
    
    const auto& result_frames = traj.getFrames();
    ASSERT_EQ(result_frames.size(), 3);
    
    // Verify positions to ensure correct frames were kept
    EXPECT_NEAR(result_frames[0].atoms()[0].position()[0], 0.0, 1e-6);
    EXPECT_NEAR(result_frames[1].atoms()[0].position()[0], 1.0, 1e-6);
    EXPECT_NEAR(result_frames[2].atoms()[0].position()[0], 2.0, 1e-6);
}

TEST_F(TrajectoryRevisionTest, HandlesNoDuplicates) {
    std::vector<Cell> frames;
    frames.push_back(createFrame(0.0, 0.0, 0.0));
    frames.push_back(createFrame(1.0, 1.0, 1.0));
    frames.push_back(createFrame(2.0, 2.0, 2.0));

    Trajectory traj(frames, 1.0);
    EXPECT_EQ(traj.getFrameCount(), 3);
}

TEST_F(TrajectoryRevisionTest, HandlesAllDuplicates) {
    std::vector<Cell> frames;
    frames.push_back(createFrame(0.0, 0.0, 0.0));
    frames.push_back(createFrame(0.0, 0.0, 0.0));
    frames.push_back(createFrame(0.0, 0.0, 0.0));

    Trajectory traj(frames, 1.0);
    // Should keep only the first one
    EXPECT_EQ(traj.getFrameCount(), 1);
}


