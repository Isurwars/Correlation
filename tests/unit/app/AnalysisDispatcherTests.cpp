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

  const auto res = correlation::app::AnalysisDispatcher::validateOptions(opts);
  EXPECT_FALSE(res.has_value());
  EXPECT_EQ(res.error(), "Error: max_ring_size must be an integer >= 3.");
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

TEST(AnalysisDispatcherTests, WriteFilesZipBundle) {
  correlation::app::TrajectoryLoader loader;
  std::string const filepath = getTestDataDir() + "car/clean.car";

  ASSERT_TRUE(loader.loadFile(filepath).has_value());

  correlation::app::AnalysisDispatcher dispatcher;
  correlation::app::ProgramOptions opts;
  opts.r_max = 5.0;
  opts.r_bin_width = 0.1;
  opts.active_calculators["RDF"] = true;

  ASSERT_TRUE(dispatcher.runAnalysis(*loader.trajectoryMut(), opts).has_value());

  opts.output_file_base = "test_dispatcher_bundle";
  opts.use_zip = true;
  opts.use_csv = true;

  const std::string zip_file = "test_dispatcher_bundle.zip";
  std::error_code error_code;
  std::filesystem::remove(zip_file, error_code);

  auto write_res = dispatcher.writeFiles(opts);
  EXPECT_TRUE(write_res.has_value()) << (write_res ? "" : write_res.error());
  EXPECT_TRUE(std::filesystem::exists(zip_file));
  std::filesystem::remove(zip_file, error_code);
}

TEST(AnalysisDispatcherTests, BondDependentCalculatorClassification) {
  // Positive classifications (bond-dependent)
  EXPECT_TRUE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("PAD"));
  EXPECT_TRUE(
      correlation::app::AnalysisDispatcher::isBondDependentCalculator("Plane Angle Distribution"));
  EXPECT_TRUE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("BAD"));
  EXPECT_TRUE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("DAD"));
  EXPECT_TRUE(correlation::app::AnalysisDispatcher::isBondDependentCalculator(
      "Dihedral Angle Distribution"));
  EXPECT_TRUE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("CN"));
  EXPECT_TRUE(
      correlation::app::AnalysisDispatcher::isBondDependentCalculator("Coordination Number"));
  EXPECT_TRUE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("RD"));
  EXPECT_TRUE(
      correlation::app::AnalysisDispatcher::isBondDependentCalculator("Ring Size Distribution"));
  EXPECT_TRUE(
      correlation::app::AnalysisDispatcher::isBondDependentCalculator("Common Neighbor Analysis"));
  EXPECT_TRUE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("CNA"));
  EXPECT_TRUE(
      correlation::app::AnalysisDispatcher::isBondDependentCalculator("Cluster Size Distribution"));
  EXPECT_TRUE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("Cluster Analysis"));
  EXPECT_TRUE(correlation::app::AnalysisDispatcher::isBondDependentCalculator(
      "Hydrogen Bond Distribution"));
  EXPECT_TRUE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("Hydrogen Bond"));

  // Negative classifications (distance-based or dynamic)
  EXPECT_FALSE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("RDF"));
  EXPECT_FALSE(correlation::app::AnalysisDispatcher::isBondDependentCalculator(
      "Radial Distribution Function"));
  EXPECT_FALSE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("g_r"));
  EXPECT_FALSE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("SQ"));
  EXPECT_FALSE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("S(Q)"));
  EXPECT_FALSE(
      correlation::app::AnalysisDispatcher::isBondDependentCalculator("Static Structure Factor"));
  EXPECT_FALSE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("XRD"));
  EXPECT_FALSE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("MSD"));
  EXPECT_FALSE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("VACF"));
  EXPECT_FALSE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("VDOS"));
  EXPECT_FALSE(
      correlation::app::AnalysisDispatcher::isBondDependentCalculator("Local Entropy Fingerprint"));
  EXPECT_FALSE(correlation::app::AnalysisDispatcher::isBondDependentCalculator("Hyperuniformity"));

  const auto names = correlation::app::AnalysisDispatcher::getBondDependentCalculatorNames();
  EXPECT_GE(names.size(), 7);
}

TEST(AnalysisDispatcherTests, DistributionFunctionsRemoveHistogram) {
  correlation::app::TrajectoryLoader loader;
  std::string const filepath = getTestDataDir() + "car/clean.car";
  ASSERT_TRUE(loader.loadFile(filepath).has_value());

  correlation::core::Cell const cell =
      loader.cell() != nullptr ? *loader.cell() : correlation::core::Cell{};
  correlation::analysis::DistributionFunctions dist_funcs(cell);

  correlation::analysis::Histogram hist;
  hist.name = "TestHist";
  dist_funcs.addHistogram("TestHist", std::move(hist));
  EXPECT_NE(dist_funcs.getHistogram("TestHist").name.empty(), true);

  EXPECT_TRUE(dist_funcs.removeHistogram("TestHist"));
  EXPECT_FALSE(dist_funcs.removeHistogram("TestHist"));
  EXPECT_THROW(static_cast<void>(dist_funcs.getHistogram("TestHist")), std::out_of_range);
}

TEST(AnalysisDispatcherTests, SelectiveRecalculationFallbackWhenNoPriorRun) {
  correlation::app::TrajectoryLoader loader;
  std::string const filepath = getTestDataDir() + "car/clean.car";
  ASSERT_TRUE(loader.loadFile(filepath).has_value());

  correlation::app::AnalysisDispatcher dispatcher;
  correlation::app::ProgramOptions opts;
  opts.r_max = 5.0;
  opts.r_bin_width = 0.1;
  opts.active_calculators["RDF"] = true;

  EXPECT_EQ(dispatcher.getDistributionFunctions(), nullptr);

  // Calling runSelectiveRecalculation with no cached df_ should fall back to full runAnalysis
  auto res = dispatcher.runSelectiveRecalculation(*loader.trajectoryMut(), opts);
  ASSERT_TRUE(res.has_value()) << (res ? "" : res.error());
  EXPECT_NE(dispatcher.getDistributionFunctions(), nullptr);
  EXPECT_NE(dispatcher.getHistogram("g_r"), nullptr);
}

TEST(AnalysisDispatcherTests, SelectiveRecalculationPreservesInvariantHistograms) {
  correlation::app::TrajectoryLoader loader;
  std::string const filepath = getTestDataDir() + "car/si_crystal.car";
  ASSERT_TRUE(loader.loadFile(filepath).has_value());

  correlation::app::AnalysisDispatcher dispatcher;
  correlation::app::ProgramOptions opts;
  opts.r_max = 5.0;
  opts.r_bin_width = 0.05;
  opts.angle_bin_width = 1.0;
  opts.active_calculators["RDF"] = true;
  opts.active_calculators["PAD"] = true;

  // First full run with bonding cutoff permitting Si-Si bonds (e.g. 0 to 3.0 Å)
  opts.bond_cutoffs[0][0] = {0.0, 3.0};
  auto run_res = dispatcher.runAnalysis(*loader.trajectoryMut(), opts);
  ASSERT_TRUE(run_res.has_value()) << (run_res ? "" : run_res.error());

  const auto *initial_gr = dispatcher.getHistogram("g_r");
  ASSERT_NE(initial_gr, nullptr);
  const auto initial_gr_partials = initial_gr->partials;

  const auto *initial_pad = dispatcher.getHistogram("PAD");
  ASSERT_NE(initial_pad, nullptr);
  const auto initial_pad_partials = initial_pad->partials;

  // Now change cutoffs to cut off all bonds (0 to 1.0 Å, shorter than Si-Si 2.35 Å)
  opts.bond_cutoffs[0][0] = {0.0, 1.0};

  // Perform selective recalculation
  auto recalc_res = dispatcher.runSelectiveRecalculation(*loader.trajectoryMut(), opts);
  ASSERT_TRUE(recalc_res.has_value()) << (recalc_res ? "" : recalc_res.error());

  // Verify distance-based g_r is intact and unchanged
  const auto *cached_gr = dispatcher.getHistogram("g_r");
  ASSERT_NE(cached_gr, nullptr);
  EXPECT_EQ(cached_gr->partials, initial_gr_partials);

  // Verify bond-dependent PAD was recalculated and updated with the new cutoff
  const auto *recalculated_pad = dispatcher.getHistogram("PAD");
  ASSERT_NE(recalculated_pad, nullptr);
  EXPECT_NE(recalculated_pad->partials, initial_pad_partials);
}

TEST(AnalysisDispatcherTests, SelectiveRecalculationNoOpWhenNoBondDependentCalculators) {
  correlation::app::TrajectoryLoader loader;
  std::string const filepath = getTestDataDir() + "car/si_crystal.car";
  ASSERT_TRUE(loader.loadFile(filepath).has_value());

  correlation::app::AnalysisDispatcher dispatcher;
  correlation::app::ProgramOptions opts;
  opts.r_max = 5.0;
  opts.r_bin_width = 0.05;
  opts.active_calculators["RDF"] = true;

  ASSERT_TRUE(dispatcher.runAnalysis(*loader.trajectoryMut(), opts).has_value());

  // Selective recalculation with only distance-based active calculator should be a no-op
  opts.active_calculators.clear();
  opts.active_calculators["RDF"] = true;
  auto recalc_res = dispatcher.runSelectiveRecalculation(*loader.trajectoryMut(), opts);
  EXPECT_TRUE(recalc_res.has_value());
}

} // namespace
