/**
 * @file TrajectoryStreamerTests.cpp
 * @brief Unit tests for TrajectoryStreamer frame scrubbing and in-place scene updates.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "render3d/Scene3D.hpp"
#include "render3d/TrajectoryStreamer.hpp"
#include <gtest/gtest.h>

namespace correlation::render3d {

TEST(TrajectoryStreamerTests, NullTrajectoryHandling) {
  TrajectoryStreamer streamer;
  EXPECT_FALSE(streamer.hasTrajectory());
  EXPECT_EQ(streamer.totalFrames(), 0);
  EXPECT_EQ(streamer.currentFrameIndex(), 0);

  Scene3D scene;
  EXPECT_FALSE(streamer.seekFrame(0, scene));
  EXPECT_FALSE(streamer.stepNext(scene));
  EXPECT_FALSE(streamer.stepPrev(scene));
}

TEST(TrajectoryStreamerTests, ScrubAndLoopFrames) {
  const math::Matrix3<real_t> lattice(math::Vector3<real_t>(10.0, 0.0, 0.0),
                                      math::Vector3<real_t>(0.0, 10.0, 0.0),
                                      math::Vector3<real_t>(0.0, 0.0, 10.0));

  std::vector<core::Cell> frames;
  for (int frame_idx = 0; frame_idx < 3; ++frame_idx) {
    core::Cell cell;
    cell.updateLattice(lattice);
    const auto coord = static_cast<real_t>(frame_idx + 1);
    cell.addAtom("C", math::Vector3<real_t>(coord, coord, coord));
    frames.push_back(cell);
  }

  const auto traj = std::make_shared<core::Trajectory>(std::move(frames), 1.0);
  TrajectoryStreamer streamer(traj);

  EXPECT_TRUE(streamer.hasTrajectory());
  EXPECT_EQ(streamer.totalFrames(), 3);

  // Initialize scene from first frame
  const auto scene = Scene3D::fromCell(traj->getFrame(0));
  ASSERT_NE(scene, nullptr);
  ASSERT_EQ(scene->atoms().size(), 1);
  EXPECT_FLOAT_EQ(scene->atoms()[0].position.x(), 1.0F);

  // Seek to frame 1
  EXPECT_TRUE(streamer.seekFrame(1, *scene));
  EXPECT_EQ(streamer.currentFrameIndex(), 1);
  EXPECT_FLOAT_EQ(scene->atoms()[0].position.x(), 2.0F);

  // Step next to frame 2
  EXPECT_TRUE(streamer.stepNext(*scene));
  EXPECT_EQ(streamer.currentFrameIndex(), 2);
  EXPECT_FLOAT_EQ(scene->atoms()[0].position.x(), 3.0F);

  // Step next with loop wraps to frame 0
  EXPECT_TRUE(streamer.stepNext(*scene, true));
  EXPECT_EQ(streamer.currentFrameIndex(), 0);
  EXPECT_FLOAT_EQ(scene->atoms()[0].position.x(), 1.0F);

  // Step prev with loop wraps to frame 2
  EXPECT_TRUE(streamer.stepPrev(*scene, true));
  EXPECT_EQ(streamer.currentFrameIndex(), 2);
  EXPECT_FLOAT_EQ(scene->atoms()[0].position.x(), 3.0F);

  // Seek out of range fails
  EXPECT_FALSE(streamer.seekFrame(99, *scene));
  EXPECT_EQ(streamer.currentFrameIndex(), 2); // Unchanged
}

} // namespace correlation::render3d
