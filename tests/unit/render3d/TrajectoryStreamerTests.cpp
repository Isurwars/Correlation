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

namespace {
std::shared_ptr<core::Trajectory> createSampleTrajectory(int frame_count) {
  const math::Matrix3<real_t> lattice(math::Vector3<real_t>(10.0, 0.0, 0.0),
                                      math::Vector3<real_t>(0.0, 10.0, 0.0),
                                      math::Vector3<real_t>(0.0, 0.0, 10.0));

  std::vector<core::Cell> frames;
  frames.reserve(static_cast<std::size_t>(frame_count));
  for (int frame_idx = 0; frame_idx < frame_count; ++frame_idx) {
    core::Cell cell;
    cell.updateLattice(lattice);
    const auto coord = static_cast<real_t>(frame_idx + 1);
    cell.addAtom("C", math::Vector3<real_t>(coord, coord, coord));
    frames.push_back(cell);
  }

  return std::make_shared<core::Trajectory>(std::move(frames), 1.0);
}
} // namespace

TEST(TrajectoryStreamerTests, ScrubFrames) {
  const auto traj = createSampleTrajectory(3);
  TrajectoryStreamer streamer(traj);

  EXPECT_TRUE(streamer.hasTrajectory());
  EXPECT_EQ(streamer.totalFrames(), 3);

  const auto scene = Scene3D::fromCell(traj->getFrame(0));
  ASSERT_NE(scene, nullptr);
  ASSERT_EQ(scene->atoms().size(), 1);
  EXPECT_FLOAT_EQ(scene->atoms()[0].position.x(), 1.0F);

  EXPECT_TRUE(streamer.seekFrame(1, *scene));
  EXPECT_EQ(streamer.currentFrameIndex(), 1);
  EXPECT_FLOAT_EQ(scene->atoms()[0].position.x(), 2.0F);

  EXPECT_TRUE(streamer.stepNext(*scene));
  EXPECT_EQ(streamer.currentFrameIndex(), 2);
  EXPECT_FLOAT_EQ(scene->atoms()[0].position.x(), 3.0F);

  EXPECT_FALSE(streamer.seekFrame(99, *scene));
  EXPECT_EQ(streamer.currentFrameIndex(), 2);
}

TEST(TrajectoryStreamerTests, LoopFrames) {
  const auto traj = createSampleTrajectory(3);
  TrajectoryStreamer streamer(traj);
  const auto scene = Scene3D::fromCell(traj->getFrame(0));
  ASSERT_NE(scene, nullptr);

  // Position at last frame (index 2)
  EXPECT_TRUE(streamer.seekFrame(2, *scene));

  // Step next with loop wraps to frame 0
  EXPECT_TRUE(streamer.stepNext(*scene, true));
  EXPECT_EQ(streamer.currentFrameIndex(), 0);
  EXPECT_FLOAT_EQ(scene->atoms()[0].position.x(), 1.0F);

  // Step prev with loop wraps back to frame 2
  EXPECT_TRUE(streamer.stepPrev(*scene, true));
  EXPECT_EQ(streamer.currentFrameIndex(), 2);
  EXPECT_FLOAT_EQ(scene->atoms()[0].position.x(), 3.0F);
}

TEST(TrajectoryStreamerTests, StreamWithInCellWrapping) {
  const math::Matrix3<real_t> lattice(math::Vector3<real_t>(10.0, 0.0, 0.0),
                                      math::Vector3<real_t>(0.0, 10.0, 0.0),
                                      math::Vector3<real_t>(0.0, 0.0, 10.0));

  std::vector<core::Cell> frames;
  for (int frame_idx = 0; frame_idx < 2; ++frame_idx) {
    core::Cell cell;
    cell.updateLattice(lattice);
    // Frame 0: inside cell (5.0)
    // Frame 1: diffused outside cell (15.5)
    const auto coord = static_cast<real_t>(frame_idx == 0 ? 5.0 : 15.5);
    cell.addAtom("C", math::Vector3<real_t>(coord, coord, coord));
    frames.push_back(cell);
  }

  const auto traj = std::make_shared<core::Trajectory>(std::move(frames), 1.0);
  TrajectoryStreamer streamer(traj);

  // Initialize scene
  const auto scene = Scene3D::fromCell(traj->getFrame(0));
  ASSERT_NE(scene, nullptr);

  // Seek to frame 1 in standard mode (in_cell = false)
  EXPECT_TRUE(streamer.seekFrame(1, *scene, false, 1.15, false));
  EXPECT_FLOAT_EQ(scene->atoms()[0].position.x(), 15.5F);

  // Seek to frame 1 with in_cell = true -> folded to 5.5
  EXPECT_TRUE(streamer.seekFrame(1, *scene, false, 1.15, true));
  EXPECT_NEAR(scene->atoms()[0].position.x(), 5.5F, 1e-4F);
}

} // namespace correlation::render3d
