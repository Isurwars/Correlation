/**
 * @file BondCutoffEngine.cpp
 * @brief Implementation of domain chemical bond cutoff calculation and validation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/BondCutoffEngine.hpp"
#include "physics/PhysicalData.hpp"

#include <cmath>
#include <stdexcept>

namespace correlation::core {

namespace {

[[nodiscard]] real_t safeGetRadius(std::string_view symbol) {
  try {
    return physics::getCovalentRadius(symbol);
  } catch (const std::out_of_range &) {
    return static_cast<real_t>(1.5); // Default covalent radius for unknown elements
  }
}

} // namespace

correlation::analysis::BondCutoffMatrix
BondCutoffEngine::compute(std::span<const Element> elements) {
  const size_t num_elements = elements.size();
  correlation::analysis::BondCutoffMatrix cutoffs(
      num_elements, std::vector<correlation::analysis::BondCutoffRange>(num_elements));

  for (size_t i = 0; i < num_elements; ++i) {
    const real_t radius_a = safeGetRadius(elements[i].symbol);
    for (size_t j = i; j < num_elements; ++j) {
      const real_t radius_b = safeGetRadius(elements[j].symbol);
      const real_t sum_radii = radius_a + radius_b;
      const real_t min_bond_dist = sum_radii * static_cast<real_t>(0.6);
      const real_t min_bond_dist_sq = min_bond_dist * min_bond_dist;
      const real_t max_bond_dist = sum_radii * static_cast<real_t>(1.3);
      const real_t max_bond_dist_sq = max_bond_dist * max_bond_dist;

      cutoffs[i][j] = correlation::analysis::BondCutoffRange{
          .min_sq = min_bond_dist_sq,
          .max_sq = max_bond_dist_sq,
      };
      cutoffs[j][i] = cutoffs[i][j];
    }
  }

  return cutoffs;
}

correlation::analysis::BondCutoffMatrix BondCutoffEngine::compute(const Cell &cell) {
  return compute(cell.elementSpan());
}

void BondCutoffEngine::validate(const correlation::analysis::BondCutoffMatrix &cutoffs) {
  for (const auto &row : cutoffs) {
    for (const auto &range : row) {
      if (range.min_sq < 0.0 || range.max_sq < 0.0) {
        throw std::invalid_argument("Bond cutoff bounds cannot be negative.");
      }
      if (range.min_sq > range.max_sq) {
        throw std::invalid_argument(
            "Minimum bond cutoff cannot be greater than maximum bond cutoff.");
      }
    }
  }
}

const correlation::analysis::BondCutoffRange &
BondCutoffEngine::getBondCutoffRange(const correlation::analysis::BondCutoffMatrix &cutoffs,
                                     size_t type1, size_t type2) {
  static const correlation::analysis::BondCutoffRange DEFAULT_RANGE{
      .min_sq = static_cast<real_t>(0.0),
      .max_sq = static_cast<real_t>(0.0),
  };
  if (type1 >= cutoffs.size() || type2 >= cutoffs[type1].size()) {
    return DEFAULT_RANGE;
  }
  return cutoffs[type1][type2];
}

real_t BondCutoffEngine::getBondCutoffSQ(const correlation::analysis::BondCutoffMatrix &cutoffs,
                                        size_t type1, size_t type2) {
  return getBondCutoffRange(cutoffs, type1, type2).max_sq;
}

real_t BondCutoffEngine::getMinBondCutoffSQ(const correlation::analysis::BondCutoffMatrix &cutoffs,
                                           size_t type1, size_t type2) {
  return getBondCutoffRange(cutoffs, type1, type2).min_sq;
}

real_t BondCutoffEngine::getBondCutoff(const correlation::analysis::BondCutoffMatrix &cutoffs,
                                      size_t type1, size_t type2) {
  return std::sqrt(getBondCutoffSQ(cutoffs, type1, type2));
}

real_t BondCutoffEngine::getMinBondCutoff(const correlation::analysis::BondCutoffMatrix &cutoffs,
                                         size_t type1, size_t type2) {
  return std::sqrt(getMinBondCutoffSQ(cutoffs, type1, type2));
}

} // namespace correlation::core
