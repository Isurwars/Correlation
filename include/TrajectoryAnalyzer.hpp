/**
 * @file TrajectoryAnalyzer.hpp
 * @brief correlation::core::Trajectory-level analysis coordination (averaging,
 * frame iteration).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "StructureAnalyzer.hpp"
#include "core/Trajectory.hpp"

#include <functional>
#include <memory>
#include <vector>

class TrajectoryAnalyzer {
public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//
  TrajectoryAnalyzer(correlation::core::Trajectory &trajectory,
                     double neighbor_cutoff,
                     const std::vector<std::vector<double>> &bond_cutoffs,
                     size_t start_frame = 0, long long end_frame = -1,
                     bool ignore_periodic_self_interactions = true,
                     std::function<void(float, const std::string &)>
                         progress_callback = nullptr);

  //-------------------------------------------------------------------------//
  //------------------------------- Accessors -------------------------------//
  //-------------------------------------------------------------------------//
  std::unique_ptr<StructureAnalyzer> createAnalyzer(size_t frame_idx) const;

  [[nodiscard]] size_t getNumFrames() const {
    return effective_end_ - start_frame_;
  }
  [[nodiscard]] size_t getStartFrame() const { return start_frame_; }

  [[nodiscard]] double getTimeStep() const { return time_step_; }
  [[nodiscard]] double getNeighborCutoff() const { return neighbor_cutoff_; }
  [[nodiscard]] const std::vector<std::vector<double>> &
  getBondCutoffsSQ() const {
    return bond_cutoffs_;
  }
  [[nodiscard]] bool getIgnorePeriodicSelfInteractions() const {
    return ignore_periodic_self_interactions_;
  }

private:
  correlation::core::Trajectory &trajectory_;
  size_t start_frame_;
  size_t effective_end_;
  double time_step_;
  double neighbor_cutoff_;
  std::vector<std::vector<double>> bond_cutoffs_;
  bool ignore_periodic_self_interactions_;
};
