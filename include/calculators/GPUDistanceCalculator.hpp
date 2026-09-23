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
[[nodiscard]] bool hasGpuDevice();

/**
 * @brief Backward-compatible inline wrapper for hasGpuDevice.
 */
[[nodiscard]] inline bool has_gpu_device() { return hasGpuDevice(); }

/**
 * @brief Performs GPU-accelerated pairwise distance computations for floating point precision T
 * (float or double).
 * @tparam T Floating-point precision type (float or double).
 * @param[in] cell Simulation unit cell.
 * @param[in] cutoff_sq Squared distance cutoff threshold.
 * @param[in] bond_cutoffs_sq Matrix of species-pair squared bond distance cutoffs.
 * @param[in] ignore_periodic_self_interactions If true, suppresses self periodic images.
 * @param[out] out_histograms Optional raw histogram tensor pointer.
 * @param[in] hist_config Histogram binning configuration.
 * @param[out] out_graph Generated NeighborGraph to populate.
 */
template <typename T = real_t>
void computeDistancesGpu(const correlation::core::Cell &cell, T cutoff_sq,
                         const std::vector<std::vector<T>> &bond_cutoffs_sq,
                         bool ignore_periodic_self_interactions, RawHistogramTensor *out_histograms,
                         DistanceCalculationConfig hist_config,
                         correlation::core::NeighborGraph &out_graph);

/**
 * @brief Backward-compatible inline wrapper for computeDistancesGpu.
 */
template <typename T = real_t>
inline void compute_distances_gpu(const correlation::core::Cell &cell, T cutoff_sq,
                                  const std::vector<std::vector<T>> &bond_cutoffs_sq,
                                  bool ignore_periodic_self_interactions,
                                  RawHistogramTensor *out_histograms,
                                  DistanceCalculationConfig hist_config,
                                  correlation::core::NeighborGraph &out_graph) {
  computeDistancesGpu<T>(cell, cutoff_sq, bond_cutoffs_sq, ignore_periodic_self_interactions,
                         out_histograms, hist_config, out_graph);
}

extern template void computeDistancesGpu<float>(
    const correlation::core::Cell &cell, float cutoff_sq,
    const std::vector<std::vector<float>> &bond_cutoffs_sq, bool ignore_periodic_self_interactions,
    RawHistogramTensor *out_histograms, DistanceCalculationConfig hist_config,
    correlation::core::NeighborGraph &out_graph);

extern template void computeDistancesGpu<double>(
    const correlation::core::Cell &cell, double cutoff_sq,
    const std::vector<std::vector<double>> &bond_cutoffs_sq, bool ignore_periodic_self_interactions,
    RawHistogramTensor *out_histograms, DistanceCalculationConfig hist_config,
    correlation::core::NeighborGraph &out_graph);

} // namespace correlation::calculators::gpu
