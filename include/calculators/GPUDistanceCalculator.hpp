/**
 * @file GPUDistanceCalculator.hpp
 * @brief GPU-accelerated pairwise distance calculator declaration.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "calculators/DistanceCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

namespace correlation::calculators::gpu {

/**
 * @brief Checks if a compatible GPU device (Nvidia or AMD) is available at runtime.
 */
bool has_gpu_device();

/**
 * @brief Performs GPU-accelerated pairwise distance computations.
 */
void compute_distances_gpu(const correlation::core::Cell &cell, real_t cutoff_sq,
                           const std::vector<std::vector<real_t>> &bond_cutoffs_sq,
                           bool ignore_periodic_self_interactions,
                           DistanceTensor &out_distances,
                           correlation::core::NeighborGraph &out_graph);

} // namespace correlation::calculators::gpu
