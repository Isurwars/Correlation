/**
 * @file TrajectoryFilter.cpp
 * @brief Implementation of trajectory frame filtering and deduplication.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/TrajectoryFilter.hpp"
#include "core/Trajectory.hpp"

namespace correlation::core {

namespace {

[[nodiscard]] bool areFramesEqual(const Cell &current, const Cell &last, real_t tolerance) {
  if (current.atomCount() != last.atomCount()) {
    return false;
  }
  if (std::abs(current.volume() - last.volume()) > tolerance) {
    return false;
  }
  const auto &curr_params = current.latticeParameters();
  const auto &last_params = last.latticeParameters();
  for (size_t k = 0; k < curr_params.size(); ++k) {
    if (std::abs(curr_params.at(k) - last_params.at(k)) > tolerance) {
      return false;
    }
  }
  const auto &current_atoms = current.atoms();
  const auto &last_atoms = last.atoms();

  for (size_t j = 0; j < current_atoms.size(); ++j) {
    if (math::norm(current_atoms[j].position() - last_atoms[j].position()) > tolerance) {
      return false;
    }
  }
  return true;
}

} // namespace

std::vector<Cell> TrajectoryFilter::filterConsecutiveDuplicates(std::span<const Cell> frames,
                                                                real_t tolerance) {
  if (frames.size() < 2) {
    return {frames.begin(), frames.end()};
  }

  std::vector<Cell> unique_frames;
  unique_frames.reserve(frames.size());
  unique_frames.push_back(frames[0]);

  for (size_t i = 1; i < frames.size(); ++i) {
    const auto &current_frame = frames[i];
    const auto &last_unique_frame = unique_frames.back();

    if (!areFramesEqual(current_frame, last_unique_frame, tolerance)) {
      unique_frames.push_back(current_frame);
    }
  }

  return unique_frames;
}

size_t TrajectoryFilter::removeDuplicatedFrames(Trajectory &trajectory, real_t tolerance) {
  const size_t initial_count = trajectory.getFrameCount();
  trajectory.removeDuplicatedFrames(tolerance);
  return initial_count - trajectory.getFrameCount();
}

} // namespace correlation::core
