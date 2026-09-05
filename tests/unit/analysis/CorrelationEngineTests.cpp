// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/CorrelationEngine.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <gtest/gtest.h>

namespace {

using correlation::analysis::CorrelationEngine;
using correlation::analysis::CorrelationEngineConfig;

TEST(CorrelationEngineTests, ValidateConfigAcceptsValidDefaults) {
  CorrelationEngineConfig const config;
  std::string const err = CorrelationEngine::validateConfig(config);
  EXPECT_TRUE(err.empty()) << "Unexpected error: " << err;
}

TEST(CorrelationEngineTests, ValidateConfigRejectsInvalidRMax) {
  CorrelationEngineConfig config;
  config.settings.r_max = -1.0;
  std::string const err = CorrelationEngine::validateConfig(config);
  EXPECT_FALSE(err.empty());
  EXPECT_NE(err.find("r_max"), std::string::npos);
}

TEST(CorrelationEngineTests, ValidateConfigRejectsInvalidRBinWidth) {
  CorrelationEngineConfig config;
  config.settings.r_max = 10.0;
  config.settings.r_bin_width = 10.0;
  std::string const err = CorrelationEngine::validateConfig(config);
  EXPECT_FALSE(err.empty());
  EXPECT_NE(err.find("r_bin_width"), std::string::npos);
}

TEST(CorrelationEngineTests, ValidateConfigRejectsInvalidRingSize) {
  CorrelationEngineConfig config;
  config.settings.max_ring_size = 2;
  std::string const err = CorrelationEngine::validateConfig(config);
  EXPECT_FALSE(err.empty());
  EXPECT_NE(err.find("max_ring_size"), std::string::npos);
}

TEST(CorrelationEngineTests, ValidateConfigRejectsInvertedFrames) {
  CorrelationEngineConfig config;
  config.min_frame = 10;
  config.max_frame = 5;
  std::string const err = CorrelationEngine::validateConfig(config);
  EXPECT_FALSE(err.empty());
  EXPECT_NE(err.find("min_frame"), std::string::npos);
}

TEST(CorrelationEngineTests, RunAnalysisAbortsOnEmptyTrajectory) {
  correlation::core::Trajectory empty_traj;
  CorrelationEngineConfig const config;

  auto const result = CorrelationEngine::runAnalysis(empty_traj, config);
  ASSERT_FALSE(result.has_value());
  EXPECT_EQ(result.error(), "Analysis aborted: No trajectory loaded.");
}

TEST(CorrelationEngineTests, RunAnalysisFailsOnInvalidConfig) {
  correlation::core::Trajectory traj;
  correlation::core::Cell cell(std::array<real_t, 6>{10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("C", correlation::math::Vector3<real_t>{0.0, 0.0, 0.0});
  traj.addFrame(cell);

  CorrelationEngineConfig config;
  config.settings.r_max = -5.0;

  auto const result = CorrelationEngine::runAnalysis(traj, config);
  ASSERT_FALSE(result.has_value());
  EXPECT_NE(result.error().find("r_max"), std::string::npos);
}

TEST(CorrelationEngineTests, RunAnalysisSucceedsOnSingleFrameTrajectory) {
  correlation::core::Trajectory traj;
  correlation::core::Cell cell(std::array<real_t, 6>{10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("C", correlation::math::Vector3<real_t>{0.0, 0.0, 0.0});
  cell.addAtom("O", correlation::math::Vector3<real_t>{1.2, 0.0, 0.0});
  traj.addFrame(cell);

  CorrelationEngineConfig config;
  config.settings.r_max = 5.0;
  config.settings.r_bin_width = 0.5;
  config.settings.active_calculators["RDF"] = true;

  bool progress_called = false;
  auto const result = CorrelationEngine::runAnalysis(
      traj, config, [&progress_called](float /*p*/, const std::string & /*msg*/) { progress_called = true; });

  ASSERT_TRUE(result.has_value()) << "runAnalysis failed: " << (result ? "" : result.error());
  ASSERT_NE(result.value(), nullptr);
  EXPECT_TRUE(progress_called);
  EXPECT_TRUE(result.value()->getAllHistograms().contains("g_r"));
}

} // namespace
