/**
 * @file GPUXRDCalculatorTests.cpp
 * @brief Unit tests for GPUXRDCalculator direct Debye diffraction and CPU fallback.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/DistributionFunctions.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "calculators/GPUXRDCalculator.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>

namespace correlation::calculators {

TEST(GPUXRDCalculatorTests, DiscoveryInCalculatorFactory) {
  const auto &factory = CalculatorFactory::instance();
  const auto *calc = factory.getCalculator("XRD — GPU Accelerated");

  ASSERT_NE(calc, nullptr);
  EXPECT_EQ(calc->getShortName(), "XRD_GPU");
  EXPECT_EQ(calc->getGroup(), "Diffraction");
  EXPECT_TRUE(calc->isFrameCalculator());
  EXPECT_FALSE(calc->isTrajectoryCalculator());
}

TEST(GPUXRDCalculatorTests, FallbackOrGPUExecution) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {2.35, 2.35, 2.35});

  correlation::analysis::DistributionFunctions dists(cell);
  correlation::analysis::AnalysisSettings const settings;

  GPUXRDCalculator const gpu_calc;
  EXPECT_NO_THROW(gpu_calc.calculateFrame(dists, settings));
  EXPECT_TRUE(dists.getAllHistograms().contains("XRD_gpu") || dists.getAllHistograms().contains("XRD"));
}

TEST(GPUXRDCalculatorTests, ComputeXRDDirectFunction) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {2.35, 2.35, 2.35});

  GPUXRDParams const params{
      .lambda = static_cast<real_t>(1.5406),
      .theta_min = static_cast<real_t>(10.0),
      .theta_max = static_cast<real_t>(80.0),
      .bin_width = static_cast<real_t>(0.5),
  };

  correlation::analysis::Histogram const hist = gpu::compute_xrd_gpu(cell, params);
  EXPECT_FALSE(hist.bins.empty());
  EXPECT_EQ(hist.x_label, "2θ");
  EXPECT_TRUE(hist.partials.contains("Total"));
  EXPECT_EQ(hist.partials.at("Total").size(), hist.bins.size());
}

TEST(GPUXRDCalculatorTests, HandlesEmptyCellGracefully) {
  correlation::core::Cell const empty_cell;
  GPUXRDParams const params{
      .lambda = static_cast<real_t>(1.5406),
      .theta_min = static_cast<real_t>(10.0),
      .theta_max = static_cast<real_t>(80.0),
      .bin_width = static_cast<real_t>(1.0),
  };

  correlation::analysis::Histogram const hist = gpu::compute_xrd_gpu(empty_cell, params);
  EXPECT_FALSE(hist.bins.empty());
  EXPECT_TRUE(hist.partials.contains("Total"));
}

TEST(GPUXRDCalculatorTests, ThrowsOnInvalidParams) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  // Invalid lambda
  EXPECT_THROW(gpu::compute_xrd_gpu(cell, {.lambda = -1.0}), std::invalid_argument);

  // Invalid theta bounds
  EXPECT_THROW(gpu::compute_xrd_gpu(cell, {.theta_min = 80.0, .theta_max = 10.0}), std::invalid_argument);

  // Invalid bin width
  EXPECT_THROW(gpu::compute_xrd_gpu(cell, {.bin_width = 0.0}), std::invalid_argument);
}

} // namespace correlation::calculators
