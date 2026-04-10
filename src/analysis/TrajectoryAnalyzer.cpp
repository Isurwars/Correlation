/**
 * @file TrajectoryAnalyzer.cpp
 * @brief Implementation of trajectory-level analysis coordination.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "analysis/TrajectoryAnalyzer.hpp"

namespace correlation::analysis {

//---------------------------------------------------------------------------//
//----------------------------- Constructors --------------------------------//
//---------------------------------------------------------------------------//

TrajectoryAnalyzer::TrajectoryAnalyzer(
    correlation::core::Trajectory &trajectory, double neighbor_cutoff,
    const std::vector<std::vector<double>> &bond_cutoffs, size_t start_frame,
    long long end_frame, bool ignore_periodic_self_interactions,
    std::function<void(float, const std::string &)> progress_callback)
    : trajectory_(trajectory), time_step_(trajectory.getTimeStep()),
      neighbor_cutoff_(neighbor_cutoff), bond_cutoffs_(bond_cutoffs),
      ignore_periodic_self_interactions_(ignore_periodic_self_interactions) {

  auto &frames = trajectory.getFrames();
  size_t n_frames = frames.size();

  effective_end_ =
      (end_frame == -1 || static_cast<size_t>(end_frame) >= n_frames)
          ? n_frames
          : static_cast<size_t>(end_frame);

  start_frame_ = (start_frame >= n_frames) ? n_frames : start_frame;

  if (progress_callback) {
    progress_callback(1.0f, "TrajectoryAnalyzer initialized.");
  }
}

std::unique_ptr<StructureAnalyzer>
TrajectoryAnalyzer::createAnalyzer(size_t frame_idx) const {
  if (frame_idx >= trajectory_.getFrames().size())
    return nullptr;
  return std::make_unique<StructureAnalyzer>(
      trajectory_.getFrames()[frame_idx], neighbor_cutoff_, bond_cutoffs_,
      ignore_periodic_self_interactions_);
}

} // namespace correlation::analysis
