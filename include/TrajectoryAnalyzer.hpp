// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <memory>
#include <vector>

#include "StructureAnalyzer.hpp"
#include "Trajectory.hpp"
#include <functional>

class TrajectoryAnalyzer {
public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//
  TrajectoryAnalyzer(Trajectory &trajectory, double neighbor_cutoff,
                     const std::vector<std::vector<double>> &bond_cutoffs,
                     size_t start_frame = 0, long long end_frame = -1,
                     bool ignore_periodic_self_interactions = true,
                     std::function<void(float)> progress_callback = nullptr);

  //-------------------------------------------------------------------------//
  //------------------------------- Accessors -------------------------------//
  //-------------------------------------------------------------------------//
  const std::vector<std::unique_ptr<StructureAnalyzer>> &getAnalyzers() const {
    return analyzers_;
  }

  double getTimeStep() const { return time_step_; }
  double getNeighborCutoff() const { return neighbor_cutoff_; }
  const std::vector<std::vector<double>> &getBondCutoffsSQ() const {
    return bond_cutoffs_;
  }
  bool getIgnorePeriodicSelfInteractions() const {
    return ignore_periodic_self_interactions_;
  }

private:
  std::vector<std::unique_ptr<StructureAnalyzer>> analyzers_;
  double time_step_;
  double neighbor_cutoff_;
  std::vector<std::vector<double>> bond_cutoffs_;
  bool ignore_periodic_self_interactions_;
};
