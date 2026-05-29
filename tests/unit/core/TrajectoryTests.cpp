// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/LinearAlgebra.hpp"
#include "readers/FileReader.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {

using namespace correlation::core;

class TrajectoryTests : public ::testing::Test {
protected:
  Cell createSimpleFrame(double x, double y, double z) {
    Cell frame({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
    frame.addAtom("H", {x, y, z});
    return frame;
  }

  void SetUp() override { std::filesystem::create_directory("test_data"); }

  void TearDown() override { std::filesystem::remove_all("test_data"); }
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
  std::vector<Cell> frames = {createSimpleFrame(0, 0, 0), createSimpleFrame(0, 0, 0), createSimpleFrame(1, 1, 1)};
  Trajectory traj;
  for (const auto &f : frames) {
    // AddFrame doesn't deduplicate automatically, but Trajectory(vector) does.
    // Let's test the explicit call.
    traj.getFrames().push_back(f);
  }
  EXPECT_EQ(traj.getFrameCount(), 3);
  traj.removeDuplicatedFrames();
  EXPECT_EQ(traj.getFrameCount(), 2);
}

TEST_F(TrajectoryTests, RemovesConsecutiveTriplicates) {
  std::vector<Cell> frames = {createSimpleFrame(0, 0, 0), createSimpleFrame(0, 0, 0), createSimpleFrame(0, 0, 0)};
  Trajectory traj(frames, 1.0);
  EXPECT_EQ(traj.getFrameCount(), 1);
}

TEST_F(TrajectoryTests, HandlesNoDuplicates) {
  std::vector<Cell> frames = {createSimpleFrame(0, 0, 0), createSimpleFrame(1, 1, 1)};
  Trajectory traj(frames, 1.0);
  EXPECT_EQ(traj.getFrameCount(), 2);
}

TEST_F(TrajectoryTests, HandlesAllDuplicates) {
  std::vector<Cell> frames = {createSimpleFrame(0, 0, 0), createSimpleFrame(0, 0, 0)};
  Trajectory traj(frames, 1.0);
  EXPECT_EQ(traj.getFrameCount(), 1);
}

// --- Unitary Tests: Physics & State ---

TEST_F(TrajectoryTests, CalculateVelocitiesComputesCorrectVelocities) {
  std::vector<Cell> frames = {createSimpleFrame(0, 0, 0), createSimpleFrame(1, 0, 0), createSimpleFrame(2, 0, 0)};
  Trajectory traj(frames, 1.0);
  traj.calculateVelocities();
  EXPECT_NEAR(traj.getFrame(0).atoms()[0].velocity().x(), 1.0, 1e-6);
}

TEST_F(TrajectoryTests, CalculateVelocitiesHandlesPBC) {
  Cell f1({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  f1.addAtom("H", {9.0, 5.0, 5.0});
  Cell f2({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  f2.addAtom("H", {1.0, 5.0, 5.0});
  Trajectory traj({f1, f2}, 1.0);
  traj.calculateVelocities();
  EXPECT_NEAR(traj.getFrame(0).atoms()[0].velocity().x(), 2.0, 1e-6);
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

TEST_F(TrajectoryTests, CalculateVelocitiesHandlesZeroOrNegativeTimeStep) {
  std::vector<Cell> frames = {createSimpleFrame(0, 0, 0), createSimpleFrame(1, 0, 0)};

  // Set pre-existing velocities to a known non-zero value to verify they're preserved
  frames[0].atoms()[0].setVelocity({42.0, 42.0, 42.0});
  frames[1].atoms()[0].setVelocity({42.0, 42.0, 42.0});

  Trajectory traj_zero(frames, 0.0);
  traj_zero.calculateVelocities();
  // With zero timestep, calculateVelocities() returns early — pre-existing velocities stay
  EXPECT_DOUBLE_EQ(traj_zero.getFrame(0).atoms()[0].velocity().x(), 42.0);
  EXPECT_DOUBLE_EQ(traj_zero.getFrame(0).atoms()[0].velocity().y(), 42.0);

  Trajectory traj_neg(frames, -0.5);
  traj_neg.calculateVelocities();
  // Same for negative timestep — function returns early
  EXPECT_DOUBLE_EQ(traj_neg.getFrame(0).atoms()[0].velocity().x(), 42.0);
}

TEST_F(TrajectoryTests, ConstructorThrowsOnMismatchedFrames) {
  std::vector<Cell> frames;
  frames.push_back(createSimpleFrame(0.0, 0.0, 0.0));
  Cell bad_frame = createSimpleFrame(1.0, 1.0, 1.0);
  bad_frame.addAtom("O", {2.0, 2.0, 2.0});
  frames.push_back(bad_frame);
  EXPECT_THROW(Trajectory traj(frames, 1.0), std::runtime_error);
}

TEST_F(TrajectoryTests, GetBondCutoffOutOfBoundsReturnsZero) {
  Cell frame({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  frame.addAtom("H", {0, 0, 0});
  Trajectory traj({frame}, 1.0);

  EXPECT_DOUBLE_EQ(traj.getBondCutoffSQ(10, 0), 0.0);
  EXPECT_DOUBLE_EQ(traj.getBondCutoff(0, 10), 0.0);
}

TEST_F(TrajectoryTests, LazyTrajectoryLoadingAndAccess) {
  std::string filename = "test_data/lazy_test.xyz";
  std::ofstream out(filename);
  out << "1\nLattice=\"10 0 0 0 10 0 0 0 10\" energy=-1.0\nSi 1.0 1.0 1.0\n";
  out << "1\nLattice=\"10 0 0 0 10 0 0 0 10\" energy=-2.0\nSi 2.0 2.0 2.0\n";
  out << "1\nLattice=\"10 0 0 0 10 0 0 0 10\" energy=-3.0\nSi 3.0 3.0 3.0\n";
  out.close();

  // Load trajectory using XYZReader (which returns a lazy trajectory)
  auto traj = correlation::readers::readTrajectory(filename, correlation::readers::FileType::Xyz);
  EXPECT_EQ(traj.getFrameCount(), 3);

  // Verify that the first frame can be accessed using firstFrame()
  const auto &first = traj.firstFrame();
  EXPECT_EQ(first.atomCount(), 1);
  EXPECT_EQ(first.atoms()[0].element().symbol, "Si");
  EXPECT_DOUBLE_EQ(first.atoms()[0].position().x(), 1.0);

  // Verify that subsequent frames can be accessed using getFrame()
  Cell f1 = traj.getFrame(1);
  EXPECT_EQ(f1.atomCount(), 1);
  EXPECT_DOUBLE_EQ(f1.atoms()[0].position().x(), 2.0);

  Cell f2 = traj.getFrame(2);
  EXPECT_EQ(f2.atomCount(), 1);
  EXPECT_DOUBLE_EQ(f2.atoms()[0].position().x(), 3.0);

  // Materialize and verify
  EXPECT_NO_THROW(traj.removeDuplicatedFrames());
  EXPECT_EQ(traj.getFrameCount(), 3);
  EXPECT_DOUBLE_EQ(traj.getFrames()[1].atoms()[0].position().x(), 2.0);
}

// --- Extreme / Edge-Case Tests ---

TEST_F(TrajectoryTests, FirstFrameThrowsOnEmptyTrajectory) {
  Trajectory traj;
  EXPECT_THROW(traj.firstFrame(), std::runtime_error);
}

TEST_F(TrajectoryTests, GetFrameThrowsOnOutOfRange) {
  Trajectory traj;
  traj.addFrame(createSimpleFrame(0, 0, 0));

  // Valid index
  EXPECT_NO_THROW(traj.getFrame(0));

  // Out of range
  EXPECT_THROW(traj.getFrame(1), std::out_of_range);
  EXPECT_THROW(traj.getFrame(100), std::out_of_range);
}

TEST_F(TrajectoryTests, CalculateVelocitiesSingleFrame) {
  // A single-frame trajectory cannot compute velocities (need at least 2 frames)
  std::vector<Cell> frames = {createSimpleFrame(5.0, 3.0, 1.0)};
  Trajectory traj(frames, 1.0);
  traj.calculateVelocities();

  // Velocity should remain at default {0,0,0} since there's nothing to diff
  EXPECT_DOUBLE_EQ(traj.getFrame(0).atoms()[0].velocity().x(), 0.0);
  EXPECT_DOUBLE_EQ(traj.getFrame(0).atoms()[0].velocity().y(), 0.0);
  EXPECT_DOUBLE_EQ(traj.getFrame(0).atoms()[0].velocity().z(), 0.0);
}

TEST_F(TrajectoryTests, RemoveDuplicatedFramesCounter) {
  std::vector<Cell> frames = {
      createSimpleFrame(0, 0, 0), createSimpleFrame(0, 0, 0), // dup
      createSimpleFrame(0, 0, 0),                             // dup
      createSimpleFrame(1, 1, 1), createSimpleFrame(1, 1, 1), // dup
  };
  Trajectory traj(frames, 1.0);

  // Constructor deduplicates: 5 frames → 2 unique
  EXPECT_EQ(traj.getFrameCount(), 2);
  EXPECT_EQ(traj.getRemovedFrameCount(), 3);
}

TEST_F(TrajectoryTests, RemoveDuplicatedFramesAlternating) {
  // Alternating pattern: no consecutive duplicates
  std::vector<Cell> frames = {
      createSimpleFrame(0, 0, 0),
      createSimpleFrame(1, 1, 1),
      createSimpleFrame(0, 0, 0),
      createSimpleFrame(1, 1, 1),
  };
  Trajectory traj(frames, 1.0);

  // No consecutive duplicates, so nothing removed
  EXPECT_EQ(traj.getFrameCount(), 4);
  EXPECT_EQ(traj.getRemovedFrameCount(), 0);
}

} // namespace correlation::testing
