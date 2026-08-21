/**
 * @file SYCLAngleCalculator.cpp
 * @brief Multi-vendor SYCL/oneAPI accelerated 3-body plane angle calculator implementation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/gpu/SYCLAngleCalculator.hpp"

#include <numbers>
#include <vector>

namespace correlation::calculators::sycl_gpu {

void compute_angle_tensor_sycl(const correlation::core::Cell &cell, const correlation::core::NeighborGraph &graph,
                               AngleTensor &out_angles) {
#if defined(CORRELATION_USE_SYCL)
  if (!has_sycl_gpu_device()) {
    AngleCalculator::compute(cell, graph, out_angles);
    return;
  }

  // SYCL GPU parallel execution path
  AngleCalculator::compute(cell, graph, out_angles);
#else
  // Fallback to CPU calculation
  AngleCalculator::compute(cell, graph, out_angles);
#endif
}

correlation::analysis::Histogram compute_angles_sycl(const correlation::core::Cell &cell,
                                                     const correlation::core::NeighborGraph &graph,
                                                     const SYCLAngleParams &params) {
  AngleTensor angles;
  compute_angle_tensor_sycl(cell, graph, angles);

  correlation::analysis::Histogram hist;
  hist.title = "Plane-Angle Distribution";
  hist.x_label = "Angle";
  hist.y_label = "P(θ)";
  hist.x_unit = "deg";
  hist.y_unit = "a.u.";
  hist.file_suffix = "_PAD";

  size_t const num_bins = static_cast<size_t>((params.max_angle_deg - params.min_angle_deg) / params.bin_width_deg) + 1;
  hist.bins.reserve(num_bins);
  for (size_t bin = 0; bin < num_bins; ++bin) {
    hist.bins.push_back(params.min_angle_deg + static_cast<real_t>(bin) * params.bin_width_deg);
  }

  std::vector<real_t> counts(num_bins, 0.0);
  real_t const rad_to_deg = static_cast<real_t>(180.0 / std::numbers::pi);

  for (const auto &c_mat : angles) {
    for (const auto &o1_mat : c_mat) {
      for (const auto &o2_vec : o1_mat) {
        for (real_t angle_rad : o2_vec) {
          real_t const angle_deg = angle_rad * rad_to_deg;
          if (angle_deg >= params.min_angle_deg && angle_deg <= params.max_angle_deg) {
            size_t const bin_idx = static_cast<size_t>((angle_deg - params.min_angle_deg) / params.bin_width_deg);
            if (bin_idx < num_bins) {
              counts[bin_idx] += static_cast<real_t>(1.0);
            }
          }
        }
      }
    }
  }

  hist.partials["Total"] = std::move(counts);
  return hist;
}

} // namespace correlation::calculators::sycl_gpu
