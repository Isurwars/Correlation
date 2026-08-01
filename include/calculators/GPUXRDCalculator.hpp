/**
 * @file GPUXRDCalculator.hpp
 * @brief GPU-accelerated X-Ray Diffraction (XRD) calculator declaration.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "calculators/gpu/SYCLXRDCalculator.hpp"
#include "core/Cell.hpp"

namespace correlation::calculators::gpu {

/**
 * @brief Performs GPU-accelerated X-Ray Diffraction (XRD) calculations with automatic fallback to CPU.
 */
inline correlation::analysis::Histogram compute_xrd_gpu(const correlation::core::Cell &cell,
                                                        const sycl_gpu::SYCLXRDParams &params = {}) {
  return sycl_gpu::compute_xrd_sycl(cell, params);
}

} // namespace correlation::calculators::gpu
