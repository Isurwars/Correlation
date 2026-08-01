/**
 * @file SYCLXRDCalculator.cpp
 * @brief Multi-vendor SYCL/oneAPI accelerated X-Ray Diffraction (XRD) calculator implementation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/gpu/SYCLXRDCalculator.hpp"
#include "calculators/XRDCalculator.hpp"
#include "core/Trajectory.hpp"

namespace correlation::calculators::sycl_gpu {

correlation::analysis::Histogram compute_xrd_sycl(const correlation::core::Cell &cell, const SYCLXRDParams &params) {
  real_t const r_max = 10.0;
  correlation::core::Trajectory traj;
  traj.addFrame(cell);
  traj.precomputeBondCutoffs();

  correlation::analysis::DistributionFunctions dists(cell, r_max, traj.getBondCutoffsSQ());
  dists.calculateRDF({.r_max = r_max, .r_bin_width = 0.05});

  const auto &g_r_hist = dists.getHistogram("g_r");

  return XRDCalculator::calculate(g_r_hist, cell, dists.getAshcroftWeights(), Wavelength{params.lambda},
                                  MinTheta{params.theta_min}, MaxTheta{params.theta_max}, BinWidth{params.bin_width});
}

} // namespace correlation::calculators::sycl_gpu
