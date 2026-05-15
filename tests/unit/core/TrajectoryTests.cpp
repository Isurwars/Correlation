// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/LinearAlgebra.hpp"
#include "readers/FileReader.hpp"

#include <gtest/gtest.h>
#include <vector>
#include <filesystem>
#include <fstream>

namespace correlation::testing {

using namespace correlation::core;

class TrajectoryTests : public ::testing::Test {
protected:
  Cell createSimpleFrame(double x, double y, double z) {
    Cell frame({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
    frame.addAtom("H", {x, y, z});
    return frame;
  }

  void SetUp() override {
    std::filesystem::create_directory("test_data");
  }

  void TearDown() override {
    std::filesystem::remove_all("test_data");
  }
};

// --- Unitary Tests: Constructors & Basic Accessors ---

TEST_F(TrajectoryTests, DefaultConstructorInitializesCorrectly) {
  Trajectory traj;
  EXPECT_EQ(traj.getFrameCount(), 0);
  EXPECT_DOUBLE_EQ(traj.getTimeStep(), 1.0);
  EXPECT_TRUE(traj.getFrames().empty());
}

TEST_F(TrajectoryTests, ParameterizedConstructorSetsFramesAndTimestep) {
  std::vector<Cell> frames;
  frames.push_back(createSimpleFrame(0.0, 0.0, 0.0));
  frames.push_back(createSimpleFrame(1.0, 1.0, 1.0));

  Trajectory traj(frames, 0.5);

  EXPECT_EQ(traj.getFrameCount(), 2);
  EXPECT_DOUBLE_EQ(traj.getTimeStep(), 0.5);
}

TEST_F(TrajectoryTests, AccessorsWorkCorrectly) {
  Trajectory traj;
  traj.setTimeStep(2.0);
  EXPECT_DOUBLE_EQ(traj.getTimeStep(), 2.0);

  Cell frame = createSimpleFrame(0, 0, 0);
  traj.addFrame(frame);
  EXPECT_EQ(traj.getFrameCount(), 1);
}

// --- Unitary Tests: Frame Management ---

TEST_F(TrajectoryTests, AddFrameAddsFrameToTrajectory) {
  Trajectory traj;
  Cell frame = createSimpleFrame(1.0, 2.0, 3.0);
  traj.addFrame(frame);

  ASSERT_EQ(traj.getFrameCount(), 1);
  EXPECT_NEAR(traj.getFrames()[0].atoms()[0].position().x(), 1.0, 1e-9);
}

TEST_F(TrajectoryTests, AddFrameThrowsOnAtomCountMismatch) {
  Trajectory traj;
  traj.addFrame(createSimpleFrame(0, 0, 0));
  Cell bad_frame({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  bad_frame.addAtom("H", {0, 0, 0});
  bad_frame.addAtom("H", {1, 1, 1});
  EXPECT_THROW(traj.addFrame(bad_frame), std::runtime_error);
}

TEST_F(TrajectoryTests, AddFrameThrowsOnElementMismatch) {
  Trajectory traj;
  traj.addFrame(createSimpleFrame(0, 0, 0));
  Cell bad_frame({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  bad_frame.addAtom("O", {0, 0, 0});
  EXPECT_THROW(traj.addFrame(bad_frame), std::runtime_error);
}

TEST_F(TrajectoryTests, AddFrameThrowsOnAtomOrderMismatch) {
  Trajectory traj;
  Cell f1({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  f1.addAtom("H", {0, 0, 0});
  f1.addAtom("O", {1, 1, 1});
  traj.addFrame(f1);

  Cell f2({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  f2.addAtom("O", {1, 1, 1});
  f2.addAtom("H", {0, 0, 0});
  EXPECT_THROW(traj.addFrame(f2), std::runtime_error);
}

TEST_F(TrajectoryTests, AddFrameThrowsOnElementCountMismatch) {
  Trajectory traj;
  Cell f1({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  f1.addAtom("H", {0, 0, 0});
  f1.addAtom("O", {1, 1, 1});
  traj.addFrame(f1);

  Cell f2({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  f2.addAtom("H", {0, 0, 0});
  f2.addAtom("H", {1, 1, 1});
  EXPECT_THROW(traj.addFrame(f2), std::runtime_error);
}

// --- Unitary Tests: Deduplication ---

TEST_F(TrajectoryTests, RemoveDuplicatedFrames) {
  std::vector<Cell> frames = {createSimpleFrame(0,0,0), createSimpleFrame(0,0,0), createSimpleFrame(1,1,1)};
  Trajectory traj;
  for(const auto& f : frames) {
    // AddFrame doesn't deduplicate automatically, but Trajectory(vector) does.
    // Let's test the explicit call.
    traj.getFrames().push_back(f);
  }
  EXPECT_EQ(traj.getFrameCount(), 3);
  traj.removeDuplicatedFrames();
  EXPECT_EQ(traj.getFrameCount(), 2);
}

TEST_F(TrajectoryTests, RemovesConsecutiveTriplicates) {
  std::vector<Cell> frames = {createSimpleFrame(0,0,0), createSimpleFrame(0,0,0), createSimpleFrame(0,0,0)};
  Trajectory traj(frames, 1.0);
  EXPECT_EQ(traj.getFrameCount(), 1);
}

TEST_F(TrajectoryTests, HandlesNoDuplicates) {
  std::vector<Cell> frames = {createSimpleFrame(0,0,0), createSimpleFrame(1,1,1)};
  Trajectory traj(frames, 1.0);
  EXPECT_EQ(traj.getFrameCount(), 2);
}

TEST_F(TrajectoryTests, HandlesAllDuplicates) {
  std::vector<Cell> frames = {createSimpleFrame(0,0,0), createSimpleFrame(0,0,0)};
  Trajectory traj(frames, 1.0);
  EXPECT_EQ(traj.getFrameCount(), 1);
}

// --- Unitary Tests: Physics & State ---

TEST_F(TrajectoryTests, CalculateVelocitiesComputesCorrectVelocities) {
  std::vector<Cell> frames = {createSimpleFrame(0,0,0), createSimpleFrame(1,0,0), createSimpleFrame(2,0,0)};
  Trajectory traj(frames, 1.0);
  traj.calculateVelocities();
  EXPECT_NEAR(traj.getVelocities()[0][0].x(), 1.0, 1e-6);
}

TEST_F(TrajectoryTests, CalculateVelocitiesHandlesPBC) {
  Cell f1({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}}); f1.addAtom("H", {9.0, 5.0, 5.0});
  Cell f2({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}}); f2.addAtom("H", {1.0, 5.0, 5.0});
  Trajectory traj({f1, f2}, 1.0);
  traj.calculateVelocities();
  EXPECT_NEAR(traj.getVelocities()[0][0].x(), 2.0, 1e-6);
}

TEST_F(TrajectoryTests, SetBondCutoffsManuallyWorks) {
  Trajectory traj;
  traj.setBondCutoffsSQ({{2.25, 4.0}, {4.0, 6.25}});
  EXPECT_DOUBLE_EQ(traj.getBondCutoffSQ(0, 0), 2.25);
}

TEST_F(TrajectoryTests, PrecomputeBondCutoffsCalculatesCorrectly) {
  Cell frame({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  frame.addAtom("H", {0, 0, 0});
  Trajectory traj({frame}, 1.0);
  traj.precomputeBondCutoffs();
  EXPECT_GT(traj.getBondCutoffSQ(0, 0), 0.0);
}

TEST_F(TrajectoryTests, CalculateVelocitiesDoesNotCrashOnEmptyTrajectory) {
  Trajectory traj;
  EXPECT_NO_THROW(traj.calculateVelocities());
}

// --- I/O Integration Tests ---

TEST_F(TrajectoryTests, ParseEnergyFromArc) {
  std::string filename = "test_data/energy_test.arc";
  std::ofstream out(filename);
  out << "!BIOSYM archive 3\nPBC=ON\n                                      -12345.6789\n!DATE\nPBC 10.0 10.0 10.0 90.0 90.0 90.0\nHe 0.0 0.0 0.0 XXXX 1 xx He 0.0\nend\nend\n";
  out.close();
  auto traj = correlation::readers::readTrajectory(filename, correlation::readers::FileType::Arc);
  EXPECT_DOUBLE_EQ(traj.getFrames()[0].getEnergy(), -12345.6789);
}

TEST_F(TrajectoryTests, ParseMultipleFramesWithEnergy) {
  std::string filename = "test_data/multi.arc";
  std::ofstream out(filename);
  out << "!BIOSYM archive 3\nPBC=ON\n -100.0\n!DATE\nPBC 10.0 10.0 10.0 90.0 90.0 90.0\nHe 0.0 0.0 0.0 XXXX 1 xx He 0.0\nend\nend\n";
  out << "!BIOSYM archive 3\nPBC=ON\n -200.5\n!DATE\nPBC 10.0 10.0 10.0 90.0 90.0 90.0\nHe 1.0 1.0 1.0 XXXX 1 xx He 0.0\nend\nend\n";
  out.close();
  auto traj = correlation::readers::readTrajectory(filename, correlation::readers::FileType::Arc);
  ASSERT_EQ(traj.getFrameCount(), 2);
  EXPECT_DOUBLE_EQ(traj.getFrames()[0].getEnergy(), -100.0);
  EXPECT_DOUBLE_EQ(traj.getFrames()[1].getEnergy(), -200.5);
}

} // namespace correlation::testing
