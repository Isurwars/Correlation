/**
 * @file TrajectoryAnalyzer.hpp
 * @brief correlation::core::Trajectory-level analysis coordination (averaging,
 * frame iteration).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "analysis/StructureAnalyzer.hpp"
#include "core/Trajectory.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace correlation::analysis {

/**
 * @class TrajectoryAnalyzer
 * @brief Orchestrates structural analysis across multiple frames of a trajectory.
 *
 * This class handles iteration over trajectory frames, providing a mechanism
 * to create frame-specific StructureAnalyzer instances. It also manages
 * high-level analysis settings such as start/end frames and periodic
 * interaction toggles.
 */
class TrajectoryAnalyzer {
public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//
  /**
   * @brief Constructs a TrajectoryAnalyzer for a given trajectory.
   *
   * @param trajectory The trajectory to analyze.
   * @param neighbor_cutoff Cutoff distance for neighbor detection (Angstrom).
   * @param bond_cutoffs Per-element-pair bond cutoffs (Angstrom).
   * @param start_frame Index of the first frame to analyze.
   * @param end_frame Index of the last frame to analyze (-1 for all).
   * @param ignore_periodic_self_interactions If true, atoms don't interact with their own periodic images.
   * @param progress_callback Callback for tracking progress (0.0 to 1.0).
   */
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
  /**
   * @brief Factory method to create a StructureAnalyzer for a specific frame.
   * @param frame_idx The frame index within the trajectory.
   * @return A unique_ptr to a configured StructureAnalyzer.
   */
  std::unique_ptr<StructureAnalyzer> createAnalyzer(size_t frame_idx) const;

  /** @return Total number of frames after accounting for start/end limits. */
  [[nodiscard]] size_t getNumFrames() const {
    return effective_end_ - start_frame_;
  }
  /** @return The configured starting frame index. */
  [[nodiscard]] size_t getStartFrame() const { return start_frame_; }

  /** @return Time step extracted from the trajectory data. */
  [[nodiscard]] double getTimeStep() const { return time_step_; }
  /** @return The global neighbor cutoff radius. */
  [[nodiscard]] double getNeighborCutoff() const { return neighbor_cutoff_; }
  /** @return The bond cutoffs used for topological analysis. */
  [[nodiscard]] const std::vector<std::vector<double>> &
  getBondCutoffsSQ() const {
    return bond_cutoffs_;
  }
  /** @return True if periodic self-interactions are being ignored. */
  [[nodiscard]] bool getIgnorePeriodicSelfInteractions() const {
    return ignore_periodic_self_interactions_;
  }

private:
  correlation::core::Trajectory &trajectory_; ///< Reference to the source trajectory.
  size_t start_frame_;                        ///< Analysis window start.
  size_t effective_end_;                      ///< Analysis window end (exclusive).
  double time_step_;                          ///< Time step between frames.
  double neighbor_cutoff_;                   ///< Neighbor search radius.
  std::vector<std::vector<double>> bond_cutoffs_; ///< Squared bond cutoffs per pair.
  bool ignore_periodic_self_interactions_;    ///< Interaction guard flag.
};

} // namespace correlation::analysis
