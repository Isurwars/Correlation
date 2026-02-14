// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <string>
#include <vector>

#include "../include/Cell.hpp"
#include "../include/FileIO.hpp"
#include "../include/Trajectory.hpp"

// Test fixture for Trajectory tests
class Test03_Trajectory : public ::testing::Test {
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

// --- Constructor Tests ---

TEST_F(Test03_Trajectory, DefaultConstructorInitializesCorrectly) {
  Trajectory traj;
  EXPECT_EQ(traj.getFrameCount(), 0);
  EXPECT_DOUBLE_EQ(traj.getTimeStep(), 1.0); // Default timestep
  EXPECT_TRUE(traj.getFrames().empty());
}

TEST_F(Test03_Trajectory, ParameterizedConstructorSetsFramesAndTimestep) {
  std::vector<Cell> frames;
  frames.push_back(createFrame(0.0, 0.0, 0.0));
  frames.push_back(createFrame(1.0, 1.0, 1.0));

  Trajectory traj(frames, 0.5);

  EXPECT_EQ(traj.getFrameCount(), 2);
  EXPECT_DOUBLE_EQ(traj.getTimeStep(), 0.5);
  EXPECT_EQ(traj.getFrames().size(), 2);
}

// --- Accessor Tests ---

TEST_F(Test03_Trajectory, AccessorsWorkCorrectly) {
  Trajectory traj;
  traj.setTimeStep(2.0);
  EXPECT_DOUBLE_EQ(traj.getTimeStep(), 2.0);

  Cell frame = createFrame(0, 0, 0);
  traj.addFrame(frame);
  EXPECT_EQ(traj.getFrameCount(), 1);
  EXPECT_EQ(traj.getFrames().size(), 1);
}

// --- Frame Management Tests ---

TEST_F(Test03_Trajectory, AddFrameAddsFrameToTrajectory) {
  Trajectory traj;
  Cell frame = createFrame(1.0, 2.0, 3.0);
  traj.addFrame(frame);

  ASSERT_EQ(traj.getFrameCount(), 1);
  EXPECT_NEAR(traj.getFrames()[0].atoms()[0].position().x(), 1.0, 1e-9);
}

// --- Velocity Calculation Tests ---

TEST_F(Test03_Trajectory, CalculateVelocitiesComputesCorrectVelocities) {
  std::vector<Cell> frames;
  // Particle moving at constant velocity (1, 0, 0) per frame
  // Time step = 1.0
  frames.push_back(createFrame(0.0, 0.0, 0.0));
  frames.push_back(createFrame(1.0, 0.0, 0.0));
  frames.push_back(createFrame(2.0, 0.0, 0.0));

  Trajectory traj(frames, 1.0);
  traj.calculateVelocities();

  const auto &velocities = traj.getVelocities();

  // Velocities are typically calculated for frames 0 to N-1 (or similar logic
  // depending on implementation) Assuming finite difference forward or central.
  // Let's check the size first.
  ASSERT_FALSE(velocities.empty());

  // Check velocity of the first atom in the first frame
  // v = (r(t+dt) - r(t)) / dt = (1.0 - 0.0) / 1.0 = 1.0
  EXPECT_NEAR(velocities[0][0].x(), 1.0, 1e-6);
  EXPECT_NEAR(velocities[0][0].y(), 0.0, 1e-6);
  EXPECT_NEAR(velocities[0][0].z(), 0.0, 1e-6);
}

TEST_F(Test03_Trajectory, CalculateVelocitiesHandlesPBC) {
  // Create a cubic lattice 10x10x10
  linalg::Vector3<double> a = {10.0, 0.0, 0.0};
  linalg::Vector3<double> b = {0.0, 10.0, 0.0};
  linalg::Vector3<double> c = {0.0, 0.0, 10.0};

  Cell frame1(a, b, c);
  frame1.addAtom("H", {9.0, 5.0, 5.0}); // Near boundary

  Cell frame2(a, b, c);
  frame2.addAtom(
      "H", {1.0, 5.0, 5.0}); // Wrapped around (moved +2.0 units effectively)

  std::vector<Cell> frames = {frame1, frame2};
  Trajectory traj(frames, 1.0);

  traj.calculateVelocities();
  const auto &velocities = traj.getVelocities();

  // Displacement should be 2.0 (9.0 -> 11.0 which wraps to 1.0)
  // Correct velocity vector logic should see the shortest path across PBC
  // 1.0 - 9.0 = -8.0.
  // With PBC (L=10): -8.0 + 10.0 = +2.0.
  EXPECT_NEAR(velocities[0][0].x(), 2.0, 1e-6);
}

// --- Bond Cutoff Tests ---

TEST_F(Test03_Trajectory, PrecomputeBondCutoffsCalculatesCorrectly) {
  // Setup a frame with H and O
  Cell frame = createFrame(0, 0, 0);
  frame.addAtom("O", {0, 0, 0}); // Adding another atom type

  std::vector<Cell> frames = {frame};
  Trajectory traj(frames, 1.0);

  // Manually set cutoffs for test (usually done via PhysicalData or auto-guess)
  // Here we are testing if the Trajectory class storage works.
  // BUT precomputeBondCutoffs uses CovalentRadii from PhysicalData usually.

  traj.precomputeBondCutoffs();

  // Element IDs: H=0, O=1 (approx, depends on registration order/data)
  // Let's verify we can retrieve something non-zero.

  // Retrieve H-H cutoff (id depends on order of insertion in Cell usually or
  // global ID) In `createFrame`, we added H. In this test we added O. The
  // element IDs are assigned sequentially in the Cell if not global. Actually
  // Element::id is global if PhysicalData is used, but here `createFrame` uses
  // string lookup. Let's rely on the fact that H and O have standard IDs if
  // PhysicalData is linked, or dynamic IDs if not. Let's check if the vector is
  // populated.

  const auto &cutoffs = traj.getBondCutoffsSQ();
  ASSERT_FALSE(cutoffs.empty());

  // Check that we have valid cutoffs for H-H, H-O, O-O
  // IDs: H is usually small, O is 8.
  // access via getBondCutoff(id1, id2)

  int id_H = frame.findElement("H")->id.value;
  int id_O = frame.findElement("O")->id.value;

  double cutoff_HH_sq = traj.getBondCutoffSQ(id_H, id_H);
  double cutoff_HO_sq = traj.getBondCutoffSQ(id_H, id_O);
  double cutoff_HH = traj.getBondCutoff(id_H, id_H);

  EXPECT_GT(cutoff_HH_sq, 0.0);
  EXPECT_GT(cutoff_HO_sq, 0.0);
  EXPECT_DOUBLE_EQ(cutoff_HH, std::sqrt(cutoff_HH_sq));
}

TEST_F(Test03_Trajectory, SetBondCutoffsManuallyWorks) {
  Trajectory traj;
  // Set a dummy 2x2 cutoff matrix with squared values
  // We want cutoffs: 1.5, 2.0, 2.5
  // So we pass: 1.5^2=2.25, 2.0^2=4.0, 2.5^2=6.25
  std::vector<std::vector<double>> manual_cutoffs = {{2.25, 4.0}, {4.0, 6.25}};
  traj.setBondCutoffsSQ(manual_cutoffs);

  EXPECT_DOUBLE_EQ(traj.getBondCutoffSQ(0, 0), 2.25);
  EXPECT_DOUBLE_EQ(traj.getBondCutoff(0, 0), 1.5);
  EXPECT_DOUBLE_EQ(traj.getBondCutoff(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(traj.getBondCutoff(1, 1), 2.5);
}

// --- Duplicate Frame Removal Tests (Existing, Renamed) ---

TEST_F(Test03_Trajectory, RemoveDuplicatedFrames) {
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

  Cell c5; // Different atoms (simulation: just partial difference or different
           // position)
  // To satisfy Trajectory validation, it must have same atom count and symbols
  // in order
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

  const auto &new_frames = traj.getFrames();
  EXPECT_EQ(new_frames[0].atoms()[0].position().x(), 0.0);
  EXPECT_EQ(new_frames[1].atoms()[0].position().x(), 0.1);
  EXPECT_EQ(new_frames[2].atoms()[0].position().x(), 0.2);
}

TEST_F(Test03_Trajectory, ParseEnergyFromArc) {
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
  const auto &frames = traj.getFrames();
  ASSERT_EQ(frames.size(), 1);

  EXPECT_DOUBLE_EQ(frames[0].getEnergy(), -12345.6789);
}

TEST_F(Test03_Trajectory, ParseMultipleFramesWithEnergy) {
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
  const auto &frames = traj2.getFrames();
  ASSERT_EQ(frames.size(), 2);
  EXPECT_DOUBLE_EQ(frames[0].getEnergy(), -100.0);
  EXPECT_DOUBLE_EQ(frames[1].getEnergy(), -200.5);
}

TEST_F(Test03_Trajectory, RemovesConsecutiveTriplicates) {
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

  const auto &result_frames = traj.getFrames();
  ASSERT_EQ(result_frames.size(), 3);

  // Verify positions to ensure correct frames were kept
  EXPECT_NEAR(result_frames[0].atoms()[0].position()[0], 0.0, 1e-6);
  EXPECT_NEAR(result_frames[1].atoms()[0].position()[0], 1.0, 1e-6);
  EXPECT_NEAR(result_frames[2].atoms()[0].position()[0], 2.0, 1e-6);
}

TEST_F(Test03_Trajectory, HandlesNoDuplicates) {
  std::vector<Cell> frames;
  frames.push_back(createFrame(0.0, 0.0, 0.0));
  frames.push_back(createFrame(1.0, 1.0, 1.0));
  frames.push_back(createFrame(2.0, 2.0, 2.0));

  Trajectory traj(frames, 1.0);
  EXPECT_EQ(traj.getFrameCount(), 3);
}

TEST_F(Test03_Trajectory, HandlesAllDuplicates) {
  std::vector<Cell> frames;
  frames.push_back(createFrame(0.0, 0.0, 0.0));
  frames.push_back(createFrame(0.0, 0.0, 0.0));
  frames.push_back(createFrame(0.0, 0.0, 0.0));

  Trajectory traj(frames, 1.0);
  // Should keep only the first one
  EXPECT_EQ(traj.getFrameCount(), 1);
}
