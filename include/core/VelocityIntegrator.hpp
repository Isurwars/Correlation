/**
 * @file VelocityIntegrator.hpp
 * @brief Domain utility for finite-difference atomic velocity calculation across trajectories.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Cell.hpp"

#include <span>

namespace correlation::core {

class Trajectory;

/**
 * @brief Domain calculator responsible for finite-difference velocity calculation under periodic
 * boundary conditions.
 */
class VelocityIntegrator {
public:
  /**
   * @brief Calculates velocities for all atoms across the provided sequence of frames using finite
   * differences.
   *
   * Utilizes forward difference on the first frame, backward difference on the last frame,
   * and central differences on all interior frames under minimum-image periodic boundary
   * conditions.
   *
   * @param[in,out] frames Span of mutable cell snapshots.
   * @param[in] time_step Constant time step between adjacent frames in femtoseconds.
   */
  static void calculateVelocities(std::span<Cell> frames, real_t time_step);

  /**
   * @brief Calculates and populates velocities on all frames within a Trajectory instance.
   * @param[in,out] trajectory Trajectory container whose frames will be updated.
   */
  static void calculateVelocities(Trajectory &trajectory);
};

} // namespace correlation::core
