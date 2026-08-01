/**
 * @file SYCLXRDCalculator.hpp
 * @brief Multi-vendor SYCL/oneAPI accelerated X-Ray Diffraction (XRD) calculator header.
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
 * @brief Parameters for SYCL XRD pattern calculation.
 */
struct SYCLXRDParams {
  real_t lambda{1.5406};  /**< Incident X-ray wavelength in Angstroms (e.g. Cu K-alpha). */
  real_t theta_min{5.0};  /**< Minimum 2-theta scattering angle in degrees. */
  real_t theta_max{90.0}; /**< Maximum 2-theta scattering angle in degrees. */
  real_t bin_width{0.1};  /**< Angular resolution step in degrees. */
};

/**
 * @brief Performs multi-vendor SYCL GPU-accelerated X-Ray Diffraction (XRD) pattern calculations.
 * @param[in] cell Simulation box cell containing atomic structure.
 * @param[in] params SYCL XRD calculation parameters.
 * @return Histogram profile containing computed XRD pattern I(2theta).
 */
correlation::analysis::Histogram compute_xrd_sycl(const correlation::core::Cell &cell,
                                                  const SYCLXRDParams &params = {});

} // namespace correlation::calculators::sycl_gpu
