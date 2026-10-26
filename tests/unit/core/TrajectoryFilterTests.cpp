// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "core/TrajectoryFilter.hpp"

#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::core;

TEST(TrajectoryFilterTests, PreservesUniqueFrames) {
  Cell frame0({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  frame0.addAtom("Si", {0.0, 0.0, 0.0});

  Cell frame1({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  frame1.addAtom("Si", {1.0, 0.0, 0.0});

  const std::vector<Cell> frames{frame0, frame1};
  const auto filtered = TrajectoryFilter::filterConsecutiveDuplicates(frames);
  EXPECT_EQ(filtered.size(), 2);
}

TEST(TrajectoryFilterTests, FiltersConsecutiveDuplicates) {
  Cell frame0({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  frame0.addAtom("Si", {0.0, 0.0, 0.0});

  const Cell frame1 = frame0; // duplicate of frame0

  Cell frame2({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  frame2.addAtom("Si", {2.0, 0.0, 0.0});

  const Cell frame3 = frame2; // duplicate of frame2

  const std::vector<Cell> frames{frame0, frame1, frame2, frame3};
  const auto filtered = TrajectoryFilter::filterConsecutiveDuplicates(frames);
  EXPECT_EQ(filtered.size(), 2);
}

TEST(TrajectoryFilterTests, RemovesDuplicatesInTrajectory) {
  Cell frame0({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  frame0.addAtom("Si", {0.0, 0.0, 0.0});

  const Cell frame1 = frame0; // duplicate

  Cell frame2({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  frame2.addAtom("Si", {2.0, 0.0, 0.0});

  // Bypass constructor deduplication by adding frame manually
  Trajectory traj;
  traj.addFrame(frame0);
  traj.addFrame(frame1);
  traj.addFrame(frame2);
  EXPECT_EQ(traj.getFrameCount(), 3);

  const size_t removed = TrajectoryFilter::removeDuplicatedFrames(traj);
  EXPECT_EQ(removed, 1);
  EXPECT_EQ(traj.getFrameCount(), 2);
}

} // namespace correlation::testing
