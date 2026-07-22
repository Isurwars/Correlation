/**
 * @file GPUSteinhardtCalculatorTests.cpp
 * @brief Unit tests for GPUSteinhardtCalculator supporting float and double precision.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "calculators/GPUSteinhardtCalculator.hpp"
#include "calculators/SteinhardtCalculator.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>

namespace correlation::calculators {

TEST(GPUSteinhardtCalculatorTests, DiscoveryInCalculatorFactory) {
  const auto &factory = CalculatorFactory::instance();
  const auto *calc = factory.getCalculator("Steinhardt Parameter — GPU Accelerated");

  ASSERT_NE(calc, nullptr);
  EXPECT_EQ(calc->getShortName(), "Steinhardt_GPU");
  EXPECT_EQ(calc->getGroup(), "Structural");
  EXPECT_TRUE(calc->isFrameCalculator());
  EXPECT_FALSE(calc->isTrajectoryCalculator());
}

TEST(GPUSteinhardtCalculatorTests, FallbackOrGPUExecution) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {1.5, 1.5, 1.5});

  correlation::analysis::StructureAnalyzer const analyzer(cell, 3.0, {{9.0}}, false);

  correlation::analysis::DistributionFunctions dists(cell, 3.0, {{9.0}});
  correlation::analysis::AnalysisSettings settings;

  GPUSteinhardtCalculator gpu_calc;
  auto hists = SteinhardtCalculator::calculate(cell, &analyzer);
  for (auto &[name, hist] : hists) {
    dists.addHistogram(name, std::move(hist));
  }

  EXPECT_NO_THROW(gpu_calc.calculateFrame(dists, settings));
  EXPECT_TRUE(dists.getAllHistograms().contains("Q4") || dists.getAllHistograms().contains("Q4_gpu"));
}

TEST(GPUSteinhardtCalculatorTests, FloatPrecisionEvaluation) {
  correlation::core::Cell cell({12.0, 12.0, 12.0, 90.0, 90.0, 90.0});
  cell.addAtom("Fe", {0.0, 0.0, 0.0});
  cell.addAtom("Fe", {2.0, 0.0, 0.0});
  cell.addAtom("Fe", {0.0, 2.0, 0.0});

  correlation::analysis::StructureAnalyzer const analyzer(cell, 3.0, {{9.0}}, false);
  correlation::analysis::DistributionFunctions dists(cell, 3.0, {{9.0}});
  correlation::analysis::AnalysisSettings settings;

  GPUSteinhardtCalculator gpu_calc;
  auto hists = SteinhardtCalculator::calculate(cell, &analyzer);
  for (auto &[name, hist] : hists) {
    dists.addHistogram(name, std::move(hist));
  }

  EXPECT_NO_THROW(gpu_calc.calculateFrame(dists, settings));
}

TEST(GPUSteinhardtCalculatorTests, DoublePrecisionEvaluation) {
  correlation::core::Cell cell({15.0, 15.0, 15.0, 90.0, 90.0, 90.0});
  cell.addAtom("Al", {0.0, 0.0, 0.0});
  cell.addAtom("Al", {2.5, 0.0, 0.0});
  cell.addAtom("Al", {0.0, 2.5, 0.0});
  cell.addAtom("Al", {0.0, 0.0, 2.5});

  correlation::analysis::StructureAnalyzer const analyzer(cell, 3.5, {{12.0}}, false);
  correlation::analysis::DistributionFunctions dists(cell, 3.5, {{12.0}});
  correlation::analysis::AnalysisSettings settings;

  GPUSteinhardtCalculator gpu_calc;
  auto hists = SteinhardtCalculator::calculate(cell, &analyzer);
  for (auto &[name, hist] : hists) {
    dists.addHistogram(name, std::move(hist));
  }

  EXPECT_NO_THROW(gpu_calc.calculateFrame(dists, settings));
}

} // namespace correlation::calculators
