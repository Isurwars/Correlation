/**
 * @file GPUSQCalculatorTests.cpp
 * @brief Unit tests for GPUSQCalculator in float and double precision.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/DistributionFunctions.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "calculators/GPUSQCalculator.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>

namespace correlation::calculators {

TEST(GPUSQCalculatorTests, DiscoveryInCalculatorFactory) {
  const auto &factory = CalculatorFactory::instance();
  const auto *calc = factory.getCalculator("S(Q) — GPU Accelerated");

  ASSERT_NE(calc, nullptr);
  EXPECT_EQ(calc->getShortName(), "S_Q_GPU");
  EXPECT_EQ(calc->getGroup(), "Scattering");
  EXPECT_TRUE(calc->isFrameCalculator());
  EXPECT_FALSE(calc->isTrajectoryCalculator());
}

TEST(GPUSQCalculatorTests, FallbackOrGPUExecution) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {1.5, 1.5, 1.5});

  correlation::analysis::DistributionFunctions dists(cell, 3.0, {{9.0}});
  correlation::analysis::AnalysisSettings settings;
  settings.q_max = 5.0;
  settings.q_bin_width = 0.5;

  GPUSQCalculator gpu_calc;
  EXPECT_NO_THROW(gpu_calc.calculateFrame(dists, settings));
  EXPECT_TRUE(dists.getAllHistograms().contains("S_q") || dists.getAllHistograms().contains("S_Q") ||
              dists.getAllHistograms().contains("S_Q_gpu"));
}

} // namespace correlation::calculators
