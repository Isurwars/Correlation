/**
 * @file SYCLDistanceCalculator.hpp
 * @brief Multi-vendor SYCL/oneAPI accelerated pairwise distance calculator header.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "calculators/spatial/DistanceCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"
#include "math/Precision.hpp"

namespace correlation::calculators::sycl_gpu {

/**
 * @brief Performs multi-vendor SYCL GPU-accelerated pairwise distance computations.
 * @param[in] cell Atomic simulation unit cell.
 * @param[in] cutoff_sq Squared distance cutoff threshold.
 * @param[in] bond_cutoffs Matrix of species-pair bond distance thresholds.
 * @param[in] ignore_periodic_self_interactions If true, suppresses self periodic images.
 * @param[out] out_graph Generated NeighborGraph to populate.
 * @param[out] out_histograms Optional raw histogram tensor pointer.
 * @param[in] hist_config Histogram binning configuration.
 */
void compute_distances_sycl(const correlation::core::Cell &cell, real_t cutoff_sq,
                            const correlation::analysis::BondCutoffMatrix &bond_cutoffs,
                            bool ignore_periodic_self_interactions,
                            correlation::core::NeighborGraph &out_graph,
                            RawHistogramTensor *out_histograms = nullptr,
                            DistanceCalculationConfig hist_config = {});

} // namespace correlation::calculators::sycl_gpu
