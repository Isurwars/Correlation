/**
 * @file TrajectoryAnalyzer.hpp
 * @brief correlation::core::Trajectory-level analysis coordination (averaging,
 * frame iteration).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
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
  /** @name Constructors */
  ///@{
  TrajectoryAnalyzer(correlation::core::Trajectory &trajectory, double neighbor_cutoff,
                     const std::vector<std::vector<double>> &bond_cutoffs, StartFrame start_frame = {0},
                     EndFrame end_frame = {static_cast<size_t>(-1)}, bool ignore_periodic_self_interactions = true,
                     const std::function<void(float, const std::string &)> &progress_callback = nullptr);

  ///@}

  /** @name Accessors */
  ///@{
  /**
   * @brief Factory method to create a StructureAnalyzer for a specific frame.
   * @param frame_idx The frame index within the trajectory.
   * @return A unique_ptr to a configured StructureAnalyzer.
   */
  std::unique_ptr<StructureAnalyzer> createAnalyzer(size_t frame_idx) const;

  /** @return Total number of frames after accounting for start/end limits. */
  [[nodiscard]] size_t getNumFrames() const { return effective_end_ - start_frame_; }
  /** @return The configured starting frame index. */
  [[nodiscard]] size_t getStartFrame() const { return start_frame_; }

  /** @return Time step extracted from the trajectory data. */
  [[nodiscard]] double getTimeStep() const { return time_step_; }
  /** @return The global neighbor cutoff radius. */
  [[nodiscard]] double getNeighborCutoff() const { return neighbor_cutoff_; }
  /** @return The bond cutoffs used for topological analysis. */
  [[nodiscard]] const std::vector<std::vector<double>> &getBondCutoffsSQ() const { return bond_cutoffs_; }
  /** @return True if periodic self-interactions are being ignored. */
  [[nodiscard]] bool getIgnorePeriodicSelfInteractions() const { return ignore_periodic_self_interactions_; }

  ///@}

private:
  correlation::core::Trajectory *trajectory_;     ///< Pointer to the source trajectory.
  size_t start_frame_;                            ///< Analysis window start.
  size_t effective_end_;                          ///< Analysis window end (exclusive).
  double time_step_;                              ///< Time step between frames.
  double neighbor_cutoff_;                        ///< Neighbor search radius.
  std::vector<std::vector<double>> bond_cutoffs_; ///< Squared bond cutoffs per pair.
  bool ignore_periodic_self_interactions_;        ///< Interaction guard flag.
};

} // namespace correlation::analysis
