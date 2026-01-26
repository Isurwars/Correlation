// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/TrajectoryAnalyzer.hpp"

#include <iostream>

TrajectoryAnalyzer::TrajectoryAnalyzer(
    Trajectory &trajectory, double neighbor_cutoff,
    const std::vector<std::vector<double>> &bond_cutoffs,
    bool ignore_periodic_self_interactions)
    : time_step_(trajectory.getTimeStep()), neighbor_cutoff_(neighbor_cutoff),
      bond_cutoffs_(bond_cutoffs),
      ignore_periodic_self_interactions_(ignore_periodic_self_interactions) {

  auto &frames = trajectory.getFrames();
  analyzers_.reserve(frames.size());

  for (auto &cell : frames) {
    if (!bond_cutoffs.empty()) {
      cell.setBondCutoffs(bond_cutoffs);
    }
    analyzers_.push_back(std::make_unique<StructureAnalyzer>(
        cell, neighbor_cutoff, ignore_periodic_self_interactions));
  }
}
