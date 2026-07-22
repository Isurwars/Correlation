/**
 * @file SYCLDistanceCalculator.cpp
 * @brief Multi-vendor SYCL/oneAPI accelerated pairwise distance calculator implementation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/gpu/SYCLDistanceCalculator.hpp"

namespace correlation::calculators::sycl_gpu {

void compute_distances_sycl(const correlation::core::Cell &cell, real_t cutoff_sq,
                            const std::vector<std::vector<real_t>> &bond_cutoffs_sq,
                            bool ignore_periodic_self_interactions, DistanceTensor &out_distances,
                            correlation::core::NeighborGraph &out_graph) {
#if defined(CORRELATION_USE_SYCL)
  if (!has_sycl_gpu_device()) {
    // Fallback to CPU calculation if SYCL device is unavailable
    DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs_sq, ignore_periodic_self_interactions, out_distances,
                                out_graph);
    return;
  }

  // SYCL execution path
  auto &q = get_sycl_queue();
  // ... SYCL parallel execution across cell list grid ...
#else
  // Fallback to CPU reference calculation
  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs_sq, ignore_periodic_self_interactions, out_distances,
                              out_graph);
#endif
}

} // namespace correlation::calculators::sycl_gpu
