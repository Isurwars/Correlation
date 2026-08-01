/**
 * @file SYCLAngleCalculator.hpp
 * @brief Multi-vendor SYCL/oneAPI accelerated 3-body bond angle calculator header.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "calculators/AngleCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"
#include "math/Precision.hpp"

namespace correlation::calculators::sycl_gpu {

/**
 * @brief Parameters for SYCL bond angle distribution calculation.
 */
struct SYCLAngleParams {
  real_t min_angle_deg{0.0};   /**< Minimum angle in degrees. */
  real_t max_angle_deg{180.0}; /**< Maximum angle in degrees. */
  real_t bin_width_deg{1.0};   /**< Angular bin width in degrees. */
};

/**
 * @brief Performs multi-vendor SYCL GPU-accelerated bond angle tensor computations.
 * @param[in] cell Simulation box cell containing atomic structure.
 * @param[in] graph Pre-computed neighbor connectivity graph.
 * @param[out] out_angles Output 4D angle tensor populated with angles in radians.
 */
void compute_angle_tensor_sycl(const correlation::core::Cell &cell, const correlation::core::NeighborGraph &graph,
                              AngleTensor &out_angles);

/**
 * @brief Performs multi-vendor SYCL GPU-accelerated bond angle distribution computation.
 * @param[in] cell Simulation box cell containing atomic structure.
 * @param[in] graph Pre-computed neighbor connectivity graph.
 * @param[in] params SYCL angle calculation parameters.
 * @return Histogram profile containing computed bond angle distribution P(theta).
 */
correlation::analysis::Histogram compute_angles_sycl(const correlation::core::Cell &cell,
                                                     const correlation::core::NeighborGraph &graph,
                                                     const SYCLAngleParams &params = {});

} // namespace correlation::calculators::sycl_gpu
