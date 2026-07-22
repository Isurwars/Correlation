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
#include "math/Precision.hpp"

namespace correlation::calculators::gpu {

/**
 * @brief Checks if a compatible GPU device (Nvidia or AMD) is available at runtime.
 */
bool has_gpu_device();

/**
 * @brief Performs GPU-accelerated pairwise distance computations for floating point precision T (float or double).
 */
template <typename T = real_t>
void compute_distances_gpu(const correlation::core::Cell &cell, T cutoff_sq,
                           const std::vector<std::vector<T>> &bond_cutoffs_sq,
                           bool ignore_periodic_self_interactions,
                           DistanceTensor &out_distances,
                           correlation::core::NeighborGraph &out_graph);

extern template void compute_distances_gpu<float>(const correlation::core::Cell &cell, float cutoff_sq,
                                                  const std::vector<std::vector<float>> &bond_cutoffs_sq,
                                                  bool ignore_periodic_self_interactions,
                                                  DistanceTensor &out_distances,
                                                  correlation::core::NeighborGraph &out_graph);

extern template void compute_distances_gpu<double>(const correlation::core::Cell &cell, double cutoff_sq,
                                                   const std::vector<std::vector<double>> &bond_cutoffs_sq,
                                                   bool ignore_periodic_self_interactions,
                                                   DistanceTensor &out_distances,
                                                   correlation::core::NeighborGraph &out_graph);

} // namespace correlation::calculators::gpu
