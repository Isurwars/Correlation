/**
 * @file TrajectoryStreamer.cpp
 * @brief Implementation of TrajectoryStreamer frame navigation and playback.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/TrajectoryStreamer.hpp"

namespace correlation::render3d {

TrajectoryStreamer::TrajectoryStreamer(std::shared_ptr<const core::Trajectory> trajectory) noexcept
    : trajectory_(std::move(trajectory)) {}

void TrajectoryStreamer::setTrajectory(
    std::shared_ptr<const core::Trajectory> trajectory) noexcept {
  trajectory_ = std::move(trajectory);
  current_frame_ = 0;
}

std::size_t TrajectoryStreamer::totalFrames() const noexcept {
  if (trajectory_ == nullptr) {
    return 0;
  }
  return trajectory_->getFrameCount();
}

bool TrajectoryStreamer::seekFrame(std::size_t frame_index, Scene3D &scene, bool recompute_bonds,
                                   double bond_cutoff_scale) {
  if (trajectory_ == nullptr || frame_index >= trajectory_->getFrameCount()) {
    return false;
  }

  const core::Cell frame = trajectory_->getFrame(frame_index);
  scene.updateFromCell(frame, recompute_bonds, bond_cutoff_scale);
  current_frame_ = frame_index;
  return true;
}

bool TrajectoryStreamer::stepNext(Scene3D &scene, bool loop, bool recompute_bonds,
                                  double bond_cutoff_scale) {
  const std::size_t total = totalFrames();
  if (total == 0) {
    return false;
  }

  if (current_frame_ + 1 < total) {
    return seekFrame(current_frame_ + 1, scene, recompute_bonds, bond_cutoff_scale);
  }

  if (loop) {
    return seekFrame(0, scene, recompute_bonds, bond_cutoff_scale);
  }

  return false;
}

bool TrajectoryStreamer::stepPrev(Scene3D &scene, bool loop, bool recompute_bonds,
                                  double bond_cutoff_scale) {
  const std::size_t total = totalFrames();
  if (total == 0) {
    return false;
  }

  if (current_frame_ > 0) {
    return seekFrame(current_frame_ - 1, scene, recompute_bonds, bond_cutoff_scale);
  }

  if (loop) {
    return seekFrame(total - 1, scene, recompute_bonds, bond_cutoff_scale);
  }

  return false;
}

} // namespace correlation::render3d
