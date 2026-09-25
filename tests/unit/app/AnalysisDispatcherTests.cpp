// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only

#include "app/services/AnalysisDispatcher.hpp"
#include "app/services/TrajectoryLoader.hpp"

#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace {

std::string getTestDataDir() {
  std::vector<std::string> const candidates = {
      "../../tests/data/",
      "../tests/data/",
      "tests/data/",
      "data/",
  };
  for (const auto &dir : candidates) {
    if (std::filesystem::exists(dir + "xyz/clean.xyz")) {
      return dir;
    }
  }
  return "../../tests/data/";
}

TEST(AnalysisDispatcherTests, InitialStateIsEmpty) {
  correlation::app::AnalysisDispatcher dispatcher;

  EXPECT_FALSE(dispatcher.isCancelled());
  EXPECT_TRUE(dispatcher.getAvailableHistogramNames().empty());
  EXPECT_EQ(dispatcher.getHistogram("g(r)"), nullptr);
  EXPECT_EQ(dispatcher.getDistributionFunctions(), nullptr);
  EXPECT_TRUE(dispatcher.getHistograms().empty());
}

TEST(AnalysisDispatcherTests, ValidateOptionsRejectsInvalidMaxRingSize) {
  correlation::app::ProgramOptions opts;
  opts.max_ring_size = 2;

  std::string err = correlation::app::AnalysisDispatcher::validateOptions(opts);
  EXPECT_FALSE(err.empty());
  EXPECT_EQ(err, "Error: max_ring_size must be an integer >= 3.");
}

TEST(AnalysisDispatcherTests, RunsAnalysisOnLoadedTrajectory) {
  correlation::app::TrajectoryLoader loader;
  std::string const filepath = getTestDataDir() + "car/clean.car";

  auto load_res = loader.loadFile(filepath);
  ASSERT_TRUE(load_res.has_value());

  correlation::app::AnalysisDispatcher dispatcher;
  correlation::app::ProgramOptions opts;
  opts.r_max = 5.0;
  opts.r_bin_width = 0.1;
  opts.active_calculators["RDF"] = true;

  auto run_res = dispatcher.runAnalysis(*loader.trajectoryMut(), opts);
  ASSERT_TRUE(run_res.has_value()) << "runAnalysis failed: " << (run_res ? "" : run_res.error());

  auto names = dispatcher.getAvailableHistogramNames();
  EXPECT_FALSE(names.empty());
  EXPECT_NE(dispatcher.getHistogram("g_r"), nullptr);
  EXPECT_NE(dispatcher.getDistributionFunctions(), nullptr);
}

TEST(AnalysisDispatcherTests, CancellationSetsFlag) {
  correlation::app::AnalysisDispatcher dispatcher;
  EXPECT_FALSE(dispatcher.isCancelled());

  dispatcher.cancelAnalysis();
  EXPECT_TRUE(dispatcher.isCancelled());
}

} // namespace
