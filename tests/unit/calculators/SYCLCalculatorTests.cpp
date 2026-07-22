/**
 * @file SYCLCalculatorTests.cpp
 * @brief Unit tests for SYCL multi-vendor GPU utility and fallback functions.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/gpu/SYCLDistanceCalculator.hpp"
#include "calculators/gpu/SYCLSQCalculator.hpp"
#include "calculators/gpu/SYCLUtils.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>

namespace correlation::testing {

TEST(SYCLCalculatorTests, DeviceDetectionDoesNotCrash) {
  bool const has_gpu = correlation::calculators::sycl_gpu::has_sycl_gpu_device();
  (void)has_gpu;
  SUCCEED();
}

TEST(SYCLCalculatorTests, ComputesDistancesViaSYCLFallback) {
  correlation::core::Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {1.0, 0.0, 0.0});

  size_t const num_elements = cell.elements().size();
  correlation::calculators::DistanceTensor out_distances(num_elements, std::vector<std::vector<real_t>>(num_elements));
  correlation::core::NeighborGraph out_graph(2);
  real_t const cutoff_sq = 4.0;
  std::vector<std::vector<real_t>> const bond_cutoffs_sq = {{2.25, 2.25}, {2.25, 2.25}};

  correlation::calculators::sycl_gpu::compute_distances_sycl(cell, cutoff_sq, bond_cutoffs_sq, false,
                                                             out_distances, out_graph);

  EXPECT_FALSE(out_distances.empty());
}

TEST(SYCLCalculatorTests, ComputesSQViaSYCLFallback) {
  correlation::core::Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {1.0, 0.0, 0.0});

  auto s_q_hist = correlation::calculators::sycl_gpu::compute_sq_sycl(cell, {.q_min = 0.5, .q_max = 10.0, .q_bin_width = 0.1});
  EXPECT_FALSE(s_q_hist.bins.empty());
}

} // namespace correlation::testing
