/**
 * @file TrajectoryFilter.hpp
 * @brief Domain utility for trajectory frame filtering and deduplication.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Cell.hpp"

#include <span>
#include <vector>

namespace correlation::core {

class Trajectory;

/**
 * @brief Domain filter responsible for removing duplicate or invalid frames from a trajectory sequence.
 */
class TrajectoryFilter {
public:
  /**
   * @brief Filters consecutive duplicated frames from a sequence of Cell snapshots.
   *
   * Two consecutive frames are considered duplicates if all atom counts match and
   * all Cartesian atom positions match within the specified distance tolerance.
   *
   * @param[in] frames Span of immutable cell frames.
   * @param[in] tolerance Euclidean distance threshold in Angstroms for position equivalence.
   * @return Vector containing only unique consecutive frames.
   */
  [[nodiscard]] static std::vector<Cell>
  filterConsecutiveDuplicates(std::span<const Cell> frames,
                              real_t tolerance = static_cast<real_t>(1e-5));

  /**
   * @brief Filters and removes consecutive duplicate frames in-place from a Trajectory.
   * @param[in,out] trajectory Trajectory container to be deduplicated.
   * @param[in] tolerance Euclidean distance threshold in Angstroms.
   * @return Number of duplicate frames removed.
   */
  static size_t removeDuplicatedFrames(Trajectory &trajectory,
                                       real_t tolerance = static_cast<real_t>(1e-5));
};

} // namespace correlation::core
