/**
 * @file BondCutoffService.hpp
 * @brief Service providing atomic pair cutoff calculations, scaling, and uniform matrices.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

namespace correlation::app {

/**
 * @class BondCutoffService
 * @brief Service calculating and managing covalent bond cutoffs on trajectories and cells.
 */
class BondCutoffService {
public:
  BondCutoffService() = default;
  ~BondCutoffService() = default;

  BondCutoffService(const BondCutoffService &) = default;
  BondCutoffService &operator=(const BondCutoffService &) = default;
  BondCutoffService(BondCutoffService &&) noexcept = default;
  BondCutoffService &operator=(BondCutoffService &&) noexcept = default;

  /**
   * @brief Calculates recommended bond cutoffs (min and max) from trajectory.
   * @param[in,out] trajectory Trajectory instance on which cutoffs are precomputed.
   * @return Computed BondCutoffMatrix or empty matrix if trajectory is null/empty.
   */
  [[nodiscard]] static correlation::analysis::BondCutoffMatrix
  getRecommendedBondCutoffs(correlation::core::Trajectory *trajectory);

  /**
   * @brief Scales recommended covalent cutoffs by a multiplicative factor and assigns them.
   * @param[in,out] trajectory Target trajectory instance.
   * @param[in] scale_factor Scaling factor (> 0).
   * @return Newly scaled BondCutoffMatrix.
   */
  static correlation::analysis::BondCutoffMatrix
  applyScaledBondCutoffs(correlation::core::Trajectory *trajectory, real_t scale_factor);

  /**
   * @brief Constructs uniform min and max cutoff distances across all atom pairs.
   * @param[in] cell Active unit cell defining element count.
   * @param[in,out] trajectory Optional trajectory to update with new cutoffs.
   * @param[in] min_cutoff Minimum cutoff distance in Å.
   * @param[in] max_cutoff Maximum cutoff distance in Å.
   * @return Symmetrized uniform BondCutoffMatrix.
   */
  static correlation::analysis::BondCutoffMatrix
  setUniformBondCutoff(const correlation::core::Cell *cell,
                       correlation::core::Trajectory *trajectory, real_t min_cutoff,
                       real_t max_cutoff);

  /**
   * @brief Gets maximum cutoff for a pair of element type indices.
   */
  [[nodiscard]] static real_t getBondCutoff(const correlation::core::Trajectory *trajectory,
                                            size_t type1, size_t type2) noexcept;

  /**
   * @brief Gets minimum cutoff for a pair of element type indices.
   */
  [[nodiscard]] static real_t getMinBondCutoff(const correlation::core::Trajectory *trajectory,
                                               size_t type1, size_t type2) noexcept;

  /**
   * @brief Sets the bond cutoffs on the target trajectory.
   */
  static void setBondCutoffs(correlation::core::Trajectory *trajectory,
                             const correlation::analysis::BondCutoffMatrix &cutoffs);
};

} // namespace correlation::app
