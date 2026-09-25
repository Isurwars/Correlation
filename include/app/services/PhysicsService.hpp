/**
 * @file PhysicsService.hpp
 * @brief Stateless physical simulation and bond cutoff calculation service.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
#include "core/Cell.hpp"
#include "math/Precision.hpp"

#include <cstddef>

namespace correlation::app {

/**
 * @class PhysicsService
 * @brief Stateless utility service providing atomic mass, time step, and cutoff calculations.
 */
class PhysicsService {
public:
  /**
   * @brief Computes a recommended simulation time step based on atomic masses in a cell.
   * @param[in] cell Pointer to the cell whose elements determine the time step.
   * @return Recommended time step in femtoseconds (fs).
   */
  [[nodiscard]] static real_t computeRecommendedTimeStep(const correlation::core::Cell *cell);

  /**
   * @brief Scales a bond cutoff matrix by a linear scaling factor.
   * @param[in] cutoffs Source bond cutoff matrix.
   * @param[in] scale_factor Multiplier applied to cutoff distances (squared internally).
   * @return Scaled bond cutoff matrix.
   */
  [[nodiscard]] static correlation::analysis::BondCutoffMatrix
  scaleBondCutoffs(const correlation::analysis::BondCutoffMatrix &cutoffs, real_t scale_factor);

  /**
   * @brief Constructs a uniform symmetric bond cutoff matrix for a given element count.
   * @param[in] num_elements Number of unique elements.
   * @param[in] min_cutoff Minimum cutoff distance in Å.
   * @param[in] max_cutoff Maximum cutoff distance in Å.
   * @return Symmetrized uniform bond cutoff matrix.
   */
  [[nodiscard]] static correlation::analysis::BondCutoffMatrix
  buildUniformBondCutoffs(size_t num_elements, real_t min_cutoff, real_t max_cutoff);
};

} // namespace correlation::app
