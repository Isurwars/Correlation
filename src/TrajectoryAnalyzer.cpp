// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "TrajectoryAnalyzer.hpp"

#include <iostream>

TrajectoryAnalyzer::TrajectoryAnalyzer(
    Trajectory &trajectory, double neighbor_cutoff,
    const std::vector<std::vector<double>> &bond_cutoffs,
    size_t start_frame, long long end_frame,
    bool ignore_periodic_self_interactions,
    std::function<void(float)> progress_callback)
    : time_step_(trajectory.getTimeStep()), neighbor_cutoff_(neighbor_cutoff),
      bond_cutoffs_(bond_cutoffs),
      ignore_periodic_self_interactions_(ignore_periodic_self_interactions) {

  auto &frames = trajectory.getFrames();
  size_t n_frames = frames.size();
  
  size_t effective_end = (end_frame == -1 || static_cast<size_t>(end_frame) >= n_frames) 
                       ? n_frames 
                       : static_cast<size_t>(end_frame);
  
  if (start_frame >= n_frames) return; // Nothing to do

  analyzers_.reserve(effective_end - start_frame);

  for (size_t i = start_frame; i < effective_end; ++i) {
    if (progress_callback) {
        float p = static_cast<float>(i - start_frame) / static_cast<float>(effective_end - start_frame);
        progress_callback(p);
    }
    analyzers_.push_back(std::make_unique<StructureAnalyzer>(
        frames[i], neighbor_cutoff, bond_cutoffs, ignore_periodic_self_interactions));
  }
  
  if (progress_callback) {
      progress_callback(1.0f);
  }
}
