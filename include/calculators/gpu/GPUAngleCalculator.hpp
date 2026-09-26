/**
 * @file GPUAngleCalculator.hpp
 * @brief GPU-accelerated 3-body bond angle calculator declaration.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "calculators/spatial/AngleCalculator.hpp"
#include "calculators/gpu/SYCLAngleCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

namespace correlation::calculators::gpu {

/**
 * @brief Performs GPU-accelerated bond angle tensor computations with automatic fallback to CPU.
 * @param[in] cell Simulation cell containing atomic coordinates and box dimensions.
 * @param[in] graph Precomputed neighbor connectivity graph.
 * @param[out] out_angles Output 4D angle tensor populated with computed angles in radians.
 */
inline void compute_angles_gpu(const correlation::core::Cell &cell,
                               const correlation::core::NeighborGraph &graph,
                               AngleTensor &out_angles) {
  correlation::calculators::sycl_gpu::compute_angle_tensor_sycl(cell, graph, out_angles);
}

} // namespace correlation::calculators::gpu
