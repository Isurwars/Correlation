/**
 * @file SYCLSQCalculator.cpp
 * @brief Multi-vendor SYCL/oneAPI accelerated S(Q) Debye scattering calculator implementation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/gpu/SYCLSQCalculator.hpp"
#include "calculators/StructureFactorCalculator.hpp"

namespace correlation::calculators::sycl_gpu {

correlation::analysis::Histogram compute_sq_sycl(const correlation::core::Cell &cell,
                                                  const SYCLSQParams &params) {
  correlation::analysis::DistributionFunctions dists(cell);
  correlation::analysis::AnalysisSettings settings;
  settings.q_max = params.q_max;
  settings.q_bin_width = params.q_bin_width;
  StructureFactorCalculator calc;
  calc.calculateFrame(dists, settings);
  return dists.getHistogram("S_q");
}

} // namespace correlation::calculators::sycl_gpu
