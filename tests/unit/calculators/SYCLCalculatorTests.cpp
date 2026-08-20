/**
 * @file SYCLCalculatorTests.cpp
 * @brief Unit tests for SYCL multi-vendor GPU utility, calculator, and fallback functions.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/GPUAngleCalculator.hpp"
#include "calculators/GPUDistanceCalculator.hpp"
#include "calculators/GPUXRDCalculator.hpp"
#include "calculators/gpu/SYCLAngleCalculator.hpp"
#include "calculators/gpu/SYCLDistanceCalculator.hpp"
#include "calculators/gpu/SYCLSQCalculator.hpp"
#include "calculators/gpu/SYCLUtils.hpp"
#include "calculators/gpu/SYCLXRDCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include <gtest/gtest.h>
#include <numbers>

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
  correlation::calculators::RawHistogramTensor out_histograms;
  correlation::core::NeighborGraph out_graph(2);
  real_t const cutoff_sq = 4.0;
  correlation::analysis::BondCutoffMatrix const bond_cutoffs = {{{0.36, 2.25}, {0.36, 2.25}},
                                                                {{0.36, 2.25}, {0.36, 2.25}}};
  correlation::calculators::DistanceCalculationConfig const hist_config{
      .r_max = 2.0,
      .r_bin_width = 0.02,
      .num_bins = 100,
  };

  correlation::calculators::sycl_gpu::compute_distances_sycl(cell, cutoff_sq, bond_cutoffs, false, out_graph,
                                                             &out_histograms, hist_config);

  EXPECT_FALSE(out_histograms.empty());
  EXPECT_EQ(out_histograms.size(), num_elements);
}

TEST(SYCLCalculatorTests, ComputesDistancesGPUWrapper) {
  correlation::core::Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {1.0, 0.0, 0.0});

  size_t const num_elements = cell.elements().size();
  correlation::calculators::RawHistogramTensor out_histograms(
      num_elements, std::vector<std::vector<real_t>>(num_elements, std::vector<real_t>(100, 0.0)));
  correlation::core::NeighborGraph out_graph(2);
  real_t const cutoff_sq = 4.0;
  std::vector<std::vector<real_t>> const bond_cutoffs_sq = {{2.25, 2.25}, {2.25, 2.25}};
  correlation::calculators::DistanceCalculationConfig const hist_config{
      .r_max = 2.0,
      .r_bin_width = 0.02,
      .num_bins = 100,
  };

  correlation::calculators::gpu::compute_distances_gpu(cell, cutoff_sq, bond_cutoffs_sq, false, &out_histograms,
                                                       hist_config, out_graph);

  EXPECT_FALSE(out_histograms.empty());
}

TEST(SYCLCalculatorTests, ComputesSQViaSYCLFallback) {
  correlation::core::Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {1.0, 0.0, 0.0});

  auto s_q_hist = correlation::calculators::sycl_gpu::compute_sq_sycl(cell, {
                                                                                .q_min = 0.5,
                                                                                .q_max = 10.0,
                                                                                .q_bin_width = 0.1,
                                                                            });
  EXPECT_FALSE(s_q_hist.bins.empty());
  EXPECT_GT(s_q_hist.bins.size(), 0U);
}

TEST(SYCLCalculatorTests, ComputesAngleTensorViaSYCLFallback) {
  correlation::core::Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {1.0, 0.0, 0.0});
  cell.addAtom("Si", {0.0, 1.0, 0.0});

  correlation::core::NeighborGraph graph(3);
  graph.addDirectedEdge(0, 1, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(0, 2, 1.0, {0.0, 1.0, 0.0});

  correlation::calculators::AngleTensor angles;
  correlation::calculators::sycl_gpu::compute_angle_tensor_sycl(cell, graph, angles);

  EXPECT_EQ(angles.size(), cell.elements().size());
  EXPECT_FALSE(angles[0][0][0].empty());

  real_t const expected_angle = static_cast<real_t>(std::numbers::pi / 2.0);
  EXPECT_NEAR(angles[0][0][0][0], expected_angle, 1e-4);
}

TEST(SYCLCalculatorTests, ComputesAnglesViaSYCLFallback) {
  correlation::core::Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {1.0, 0.0, 0.0});
  cell.addAtom("Si", {0.0, 1.0, 0.0});

  correlation::core::NeighborGraph graph(3);
  graph.addDirectedEdge(0, 1, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(0, 2, 1.0, {0.0, 1.0, 0.0});

  auto angle_hist = correlation::calculators::sycl_gpu::compute_angles_sycl(cell, graph,
                                                                            {
                                                                                .min_angle_deg = 0.0,
                                                                                .max_angle_deg = 180.0,
                                                                                .bin_width_deg = 1.0,
                                                                            });
  EXPECT_FALSE(angle_hist.bins.empty());
  EXPECT_TRUE(angle_hist.partials.contains("Total"));
}

TEST(SYCLCalculatorTests, ComputesAnglesGPUWrapper) {
  correlation::core::Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {1.0, 0.0, 0.0});
  cell.addAtom("Si", {0.0, 1.0, 0.0});

  correlation::core::NeighborGraph graph(3);
  graph.addDirectedEdge(0, 1, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(0, 2, 1.0, {0.0, 1.0, 0.0});

  correlation::calculators::AngleTensor angles;
  correlation::calculators::gpu::compute_angles_gpu(cell, graph, angles);

  EXPECT_EQ(angles.size(), cell.elements().size());
}

TEST(SYCLCalculatorTests, ComputesXRDViaSYCLFallback) {
  correlation::core::Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {1.0, 0.0, 0.0});

  auto xrd_hist = correlation::calculators::sycl_gpu::compute_xrd_sycl(cell, {
                                                                                 .lambda = 1.5406,
                                                                                 .theta_min = 10.0,
                                                                                 .theta_max = 80.0,
                                                                                 .bin_width = 0.5,
                                                                             });
  EXPECT_FALSE(xrd_hist.bins.empty());
  EXPECT_EQ(xrd_hist.x_label, "2θ");
}

TEST(SYCLCalculatorTests, ComputesXRDGPUWrapper) {
  correlation::core::Cell cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {1.0, 0.0, 0.0});

  auto xrd_hist = correlation::calculators::gpu::compute_xrd_gpu(cell, {
                                                                           .lambda = 1.5406,
                                                                           .theta_min = 10.0,
                                                                           .theta_max = 80.0,
                                                                           .bin_width = 0.5,
                                                                       });
  EXPECT_FALSE(xrd_hist.bins.empty());
}

TEST(SYCLCalculatorTests, SYCLAngleParamsDefaultValues) {
  correlation::calculators::sycl_gpu::SYCLAngleParams params;
  EXPECT_NEAR(params.min_angle_deg, static_cast<real_t>(0.0), 1e-4);
  EXPECT_NEAR(params.max_angle_deg, static_cast<real_t>(180.0), 1e-4);
  EXPECT_NEAR(params.bin_width_deg, static_cast<real_t>(1.0), 1e-4);
}

TEST(SYCLCalculatorTests, SYCLXRDParamsDefaultValues) {
  correlation::calculators::sycl_gpu::SYCLXRDParams params;
  EXPECT_NEAR(params.lambda, static_cast<real_t>(1.5406), 1e-4);
  EXPECT_NEAR(params.theta_min, static_cast<real_t>(5.0), 1e-4);
  EXPECT_NEAR(params.theta_max, static_cast<real_t>(90.0), 1e-4);
  EXPECT_NEAR(params.bin_width, static_cast<real_t>(0.1), 1e-4);
}

} // namespace correlation::testing
