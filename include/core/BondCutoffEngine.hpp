/**
 * @file BondCutoffEngine.hpp
 * @brief Domain utility for chemical bond cutoff calculation and validation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
#include "core/Cell.hpp"

#include <span>

namespace correlation::core {

/**
 * @brief Domain engine responsible for computing, validating, and querying chemical bond cutoffs.
 */
class BondCutoffEngine {
public:
  /**
   * @brief Pre-calculates the squared bond cutoff distances for every pair of chemical elements.
   * @param[in] elements List or span of unique chemical elements.
   * @return Matrix of calculated BondCutoffRange entries.
   */
  [[nodiscard]] static correlation::analysis::BondCutoffMatrix
  compute(std::span<const Element> elements);

  /**
   * @brief Pre-calculates the bond cutoff matrix from the elements present in a Cell.
   * @param[in] cell The structure containing unique chemical elements.
   * @return Matrix of calculated BondCutoffRange entries.
   */
  [[nodiscard]] static correlation::analysis::BondCutoffMatrix compute(const Cell &cell);

  /**
   * @brief Validates a bond cutoff matrix against domain bounds.
   * @param[in] cutoffs Matrix of bond cutoff ranges to inspect.
   * @throws std::invalid_argument If any cutoff bound is negative or min_sq > max_sq.
   */
  static void validate(const correlation::analysis::BondCutoffMatrix &cutoffs);

  /**
   * @brief Gets the bond cutoff range for two given element type indices.
   * @param[in] cutoffs Matrix of bond cutoff ranges.
   * @param[in] type1 ID of first element type.
   * @param[in] type2 ID of second element type.
   * @return Const reference to the BondCutoffRange, or a zeroed default range if out of bounds.
   */
  [[nodiscard]] static const correlation::analysis::BondCutoffRange &
  getBondCutoffRange(const correlation::analysis::BondCutoffMatrix &cutoffs, size_t type1,
                     size_t type2);

  /**
   * @brief Gets the maximum squared bond cutoff distance for two given element types.
   * @param[in] cutoffs Matrix of bond cutoff ranges.
   * @param[in] type1 ID of first element type.
   * @param[in] type2 ID of second element type.
   * @return The squared maximum bond cutoff distance.
   */
  [[nodiscard]] static real_t
  getBondCutoffSQ(const correlation::analysis::BondCutoffMatrix &cutoffs, size_t type1,
                  size_t type2);

  /**
   * @brief Gets the minimum squared bond cutoff distance for two given element types.
   * @param[in] cutoffs Matrix of bond cutoff ranges.
   * @param[in] type1 ID of first element type.
   * @param[in] type2 ID of second element type.
   * @return The squared minimum bond cutoff distance.
   */
  [[nodiscard]] static real_t
  getMinBondCutoffSQ(const correlation::analysis::BondCutoffMatrix &cutoffs, size_t type1,
                     size_t type2);

  /**
   * @brief Gets the maximum bond cutoff distance for two given element types.
   * @param[in] cutoffs Matrix of bond cutoff ranges.
   * @param[in] type1 ID of first element type.
   * @param[in] type2 ID of second element type.
   * @return The maximum bond cutoff distance.
   */
  [[nodiscard]] static real_t getBondCutoff(const correlation::analysis::BondCutoffMatrix &cutoffs,
                                            size_t type1, size_t type2);

  /**
   * @brief Gets the minimum bond cutoff distance for two given element types.
   * @param[in] cutoffs Matrix of bond cutoff ranges.
   * @param[in] type1 ID of first element type.
   * @param[in] type2 ID of second element type.
   * @return The minimum bond cutoff distance.
   */
  [[nodiscard]] static real_t
  getMinBondCutoff(const correlation::analysis::BondCutoffMatrix &cutoffs, size_t type1,
                   size_t type2);
};

} // namespace correlation::core
