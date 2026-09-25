/**
 * @file BondCutoffService.cpp
 * @brief Implementation of the BondCutoffService.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/services/BondCutoffService.hpp"
#include "app/services/PhysicsService.hpp"

namespace correlation::app {

correlation::analysis::BondCutoffMatrix
BondCutoffService::getRecommendedBondCutoffs(correlation::core::Trajectory *trajectory) {
  if (trajectory == nullptr || trajectory->getFrameCount() == 0) {
    return {};
  }
  trajectory->precomputeBondCutoffs();
  return trajectory->getBondCutoffs();
}

correlation::analysis::BondCutoffMatrix
BondCutoffService::applyScaledBondCutoffs(correlation::core::Trajectory *trajectory,
                                          real_t scale_factor) {
  auto cutoffs =
      PhysicsService::scaleBondCutoffs(getRecommendedBondCutoffs(trajectory), scale_factor);
  setBondCutoffs(trajectory, cutoffs);
  return cutoffs;
}

correlation::analysis::BondCutoffMatrix
BondCutoffService::setUniformBondCutoff(const correlation::core::Cell *cell,
                                        correlation::core::Trajectory *trajectory,
                                        real_t min_cutoff, real_t max_cutoff) {
  if (cell == nullptr) {
    return {};
  }
  auto cutoffs =
      PhysicsService::buildUniformBondCutoffs(cell->elements().size(), min_cutoff, max_cutoff);
  setBondCutoffs(trajectory, cutoffs);
  return cutoffs;
}

real_t BondCutoffService::getBondCutoff(const correlation::core::Trajectory *trajectory,
                                        size_t type1, size_t type2) noexcept {
  if (trajectory == nullptr) {
    return 0.0;
  }
  return trajectory->getBondCutoff(type1, type2);
}

real_t BondCutoffService::getMinBondCutoff(const correlation::core::Trajectory *trajectory,
                                           size_t type1, size_t type2) noexcept {
  if (trajectory == nullptr) {
    return 0.0;
  }
  return trajectory->getMinBondCutoff(type1, type2);
}

void BondCutoffService::setBondCutoffs(correlation::core::Trajectory *trajectory,
                                       const correlation::analysis::BondCutoffMatrix &cutoffs) {
  if (trajectory != nullptr) {
    trajectory->setBondCutoffs(cutoffs);
  }
}

} // namespace correlation::app
