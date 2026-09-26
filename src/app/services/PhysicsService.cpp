/**
 * @file PhysicsService.cpp
 * @brief Implementation of stateless physical simulation and bond cutoff calculation service.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/services/PhysicsService.hpp"
#include "app/core/AppOptions.hpp"
#include "physics/PhysicalData.hpp"

#include <cmath>
#include <iostream>
#include <limits>

namespace correlation::app {

real_t PhysicsService::computeRecommendedTimeStep(const correlation::core::Cell *cell) {
  if ((cell == nullptr) || cell->elements().empty()) {
    return AppDefaults::TIME_STEP;
  }

  real_t min_mass = std::numeric_limits<real_t>::max();
  bool found = false;

  for (const auto &element : cell->elements()) {
    try {
      const real_t mass = correlation::physics::getAtomicMass(element.symbol);
      if (mass < min_mass) {
        min_mass = mass;
        found = true;
      }
    } catch (const std::out_of_range &err) {
      std::cerr << "Warning: Unknown element symbol '" << element.symbol
                << "' ignored in mass calculation: " << err.what() << '\n';
    }
  }

  if (found && min_mass > static_cast<real_t>(0.0)) {
    return static_cast<real_t>(std::sqrt(9.0 * min_mass / 5.0));
  }

  return AppDefaults::TIME_STEP;
}

correlation::analysis::BondCutoffMatrix
PhysicsService::scaleBondCutoffs(const correlation::analysis::BondCutoffMatrix &cutoffs,
                                 real_t scale_factor) {
  if (cutoffs.empty() || scale_factor <= static_cast<real_t>(0.0)) {
    return cutoffs;
  }

  auto scaled = cutoffs;
  const real_t factor_sq = scale_factor * scale_factor;
  for (auto &row : scaled) {
    for (auto &range : row) {
      range.min_sq *= factor_sq;
      range.max_sq *= factor_sq;
    }
  }
  return scaled;
}

correlation::analysis::BondCutoffMatrix
PhysicsService::buildUniformBondCutoffs(size_t num_elements, real_t min_cutoff, real_t max_cutoff) {
  if (num_elements == 0) {
    return {};
  }

  const real_t min_sq = (min_cutoff > static_cast<real_t>(0.0)) ? (min_cutoff * min_cutoff)
                                                                : static_cast<real_t>(0.0);
  const real_t max_sq = (max_cutoff > min_cutoff) ? (max_cutoff * max_cutoff) : min_sq;

  return correlation::analysis::BondCutoffMatrix(
      num_elements, std::vector<correlation::analysis::BondCutoffRange>(
                        num_elements, {.min_sq = min_sq, .max_sq = max_sq}));
}

} // namespace correlation::app
