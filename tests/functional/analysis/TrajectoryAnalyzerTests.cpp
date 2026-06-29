// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/TrajectoryAnalyzer.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <gtest/gtest.h>

namespace correlation::analysis {

TEST(TrajectoryAnalyzerTests, BasicUsage) {
  // Create a dummy cell
  correlation::core::Cell cell;
  cell.setLatticeParameters({10, 10, 10, 90, 90, 90});
  cell.addAtom("H", {0.5, 0.5, 0.5});
  cell.addAtom("H", {1.5, 0.5, 0.5}); // Distance 1.0

  // Create a trajectory with multiple frames
  correlation::core::Trajectory trajectory;
  trajectory.addFrame(cell);
  trajectory.addFrame(cell);

  double const neighbor_cutoff = 2.0;
  std::vector<std::vector<double>> const bond_cutoffs = {{1.1}}; // Assuming H-H is index 0-0, simplified

  TrajectoryAnalyzer const analyzer(trajectory, neighbor_cutoff, bond_cutoffs);

  EXPECT_EQ(analyzer.getNumFrames(), 2);
  EXPECT_EQ(analyzer.getNeighborCutoff(), neighbor_cutoff);
  EXPECT_EQ(analyzer.getBondCutoffsSQ()[0][0], 1.1);

  // Verify that StructureAnalyzers can be created dynamically
  for (size_t i = 0; i < analyzer.getNumFrames(); ++i) {
    auto analyzer_ptr = analyzer.createAnalyzer(i);
    // Just checking if we can access them and they are not null
    EXPECT_TRUE(analyzer_ptr != nullptr);
    // In a real test we might check neighbor counts if we knew indices
  }
}

TEST(TrajectoryAnalyzerTests, StartAndEndFrameLimits) {
  correlation::core::Cell cell;
  cell.setLatticeParameters({10, 10, 10, 90, 90, 90});
  cell.addAtom("H", {0.5, 0.5, 0.5});

  correlation::core::Trajectory trajectory;
  trajectory.setTimeStep(2.5); // 2.5 fs
  trajectory.addFrame(cell);
  trajectory.addFrame(cell);
  trajectory.addFrame(cell);
  trajectory.addFrame(cell); // 4 frames total

  double const neighbor_cutoff = 2.0;
  std::vector<std::vector<double>> const bond_cutoffs = {{1.1}};

  // Case 1: Analyze frames from index 1 to 3 (exclusive, so frames 1 and 2, getNumFrames() = 2)
  {
    TrajectoryAnalyzer const analyzer(trajectory, neighbor_cutoff, bond_cutoffs, StartFrame{1}, EndFrame{3});
    EXPECT_EQ(analyzer.getStartFrame(), 1);
    EXPECT_EQ(analyzer.getNumFrames(), 2);
    EXPECT_DOUBLE_EQ(analyzer.getTimeStep(), 2.5);
    EXPECT_TRUE(analyzer.getIgnorePeriodicSelfInteractions());
  }

  // Case 2: Analyze all frames using default/out-of-bounds limits
  {
    TrajectoryAnalyzer const analyzer(trajectory, neighbor_cutoff, bond_cutoffs, StartFrame{0}, EndFrame{999}, false);
    EXPECT_EQ(analyzer.getStartFrame(), 0);
    EXPECT_EQ(analyzer.getNumFrames(), 4);
    EXPECT_FALSE(analyzer.getIgnorePeriodicSelfInteractions());
  }

  // Case 3: Start frame out of bounds (should clamp to n_frames and start == end == 4)
  {
    TrajectoryAnalyzer const analyzer(trajectory, neighbor_cutoff, bond_cutoffs, StartFrame{10}, EndFrame{15});
    EXPECT_EQ(analyzer.getStartFrame(), 4);
    EXPECT_EQ(analyzer.getNumFrames(), 0);
  }
}

TEST(TrajectoryAnalyzerTests, ProgressCallbackIsCalled) {
  correlation::core::Cell cell;
  cell.setLatticeParameters({10, 10, 10, 90, 90, 90});
  cell.addAtom("H", {0.5, 0.5, 0.5});

  correlation::core::Trajectory trajectory;
  trajectory.addFrame(cell);

  double const neighbor_cutoff = 2.0;
  std::vector<std::vector<double>> const bond_cutoffs = {{1.1}};

  bool callback_called = false;
  float progress_val = 0.0F;
  std::string progress_msg;

  auto callback = [&](float progress, const std::string &message) {
    callback_called = true;
    progress_val = progress;
    progress_msg = message;
  };

  TrajectoryAnalyzer const analyzer(trajectory, neighbor_cutoff, bond_cutoffs, StartFrame{0}, EndFrame{1}, true,
                                    callback);

  EXPECT_TRUE(callback_called);
  EXPECT_EQ(progress_val, 1.0F);
  EXPECT_EQ(progress_msg, "TrajectoryAnalyzer initialized.");
}

TEST(TrajectoryAnalyzerTests, CreateAnalyzerOutOfBoundsReturnsNullptr) {
  correlation::core::Cell cell;
  cell.setLatticeParameters({10, 10, 10, 90, 90, 90});
  cell.addAtom("H", {0.5, 0.5, 0.5});

  correlation::core::Trajectory trajectory;
  trajectory.addFrame(cell); // 1 frame total

  double const neighbor_cutoff = 2.0;
  std::vector<std::vector<double>> const bond_cutoffs = {{1.1}};

  TrajectoryAnalyzer const analyzer(trajectory, neighbor_cutoff, bond_cutoffs);

  // Index 0 is valid
  auto analyzer_ptr = analyzer.createAnalyzer(0);
  EXPECT_NE(analyzer_ptr, nullptr);

  // Index 1 is out of bounds
  auto analyzer_invalid = analyzer.createAnalyzer(1);
  EXPECT_EQ(analyzer_invalid, nullptr);
}

} // namespace correlation::analysis
