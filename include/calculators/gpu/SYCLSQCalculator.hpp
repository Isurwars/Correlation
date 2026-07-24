/**
 * @file SYCLSQCalculator.hpp
 * @brief Multi-vendor SYCL/oneAPI accelerated S(Q) Debye scattering calculator header.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "core/Cell.hpp"
#include "math/Precision.hpp"

namespace correlation::calculators::sycl_gpu {

/**
 * @brief Parameters for SYCL S(Q) structure factor calculation.
 */
struct SYCLSQParams {
  real_t q_min{0.0};       /**< Minimum Q value for structure factor computation. */
  real_t q_max{10.0};      /**< Maximum Q value for structure factor computation. */
  real_t q_bin_width{0.1}; /**< Bin resolution step for Q values. */
};

/**
 * @brief Performs multi-vendor SYCL GPU-accelerated Debye scattering S(Q) calculations.
 * @param[in] cell Simulation box cell containing atomic structure.
 * @param[in] params SYCL S(Q) calculation parameters.
 * @return Histogram profile containing computed S(Q) data.
 */
correlation::analysis::Histogram compute_sq_sycl(const correlation::core::Cell &cell, const SYCLSQParams &params = {});

} // namespace correlation::calculators::sycl_gpu
