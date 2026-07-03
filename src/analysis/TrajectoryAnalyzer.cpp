/**
 * @file TrajectoryAnalyzer.cpp
 * @brief Implementation of trajectory-level analysis coordination.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/TrajectoryAnalyzer.hpp"
#include <algorithm>

namespace correlation::analysis {


TrajectoryAnalyzer::TrajectoryAnalyzer(correlation::core::Trajectory &trajectory, double neighbor_cutoff,
                                       const std::vector<std::vector<double>> &bond_cutoffs, StartFrame start_frame,
                                       EndFrame end_frame, bool ignore_periodic_self_interactions,
                                       const std::function<void(float, const std::string &)> &progress_callback)
    : trajectory_(&trajectory), time_step_(trajectory.getTimeStep()), neighbor_cutoff_(neighbor_cutoff),
      bond_cutoffs_(bond_cutoffs), ignore_periodic_self_interactions_(ignore_periodic_self_interactions) {

  size_t n_frames = trajectory.getFrameCount();

  effective_end_ =
      (end_frame.value == static_cast<size_t>(-1) || end_frame.value >= n_frames) ? n_frames : end_frame.value;

  start_frame_ = (start_frame.value >= n_frames) ? n_frames : start_frame.value;

  start_frame_ = std::min(start_frame_, effective_end_);

  if (progress_callback) {
    progress_callback(1.0F, "TrajectoryAnalyzer initialized.");
  }
}

std::unique_ptr<StructureAnalyzer> TrajectoryAnalyzer::createAnalyzer(size_t frame_idx) const {
  if (frame_idx >= trajectory_->getFrameCount()) {
    return nullptr;
  }
  return std::make_unique<StructureAnalyzer>(trajectory_->getFrame(frame_idx), neighbor_cutoff_, bond_cutoffs_,
                                             ignore_periodic_self_interactions_);
}

} // namespace correlation::analysis
