/**
 * @file PlotExportManagerTests.cpp
 * @brief Unit tests for PlotExportManager.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/DistributionFunctions.hpp"
#include "app/viewmodel/PlotExportManager.hpp"
#include "plotters/PlotTypes.hpp"
#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

namespace correlation::app {
namespace {

correlation::analysis::Histogram createMockHistogram(const std::string &title) {
  correlation::analysis::Histogram hist;
  hist.title = title + " Test Title";
  hist.x_label = "r";
  hist.x_unit = "Å";
  hist.y_label = "g(r)";
  hist.bins = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
  hist.partials["Total"] = {0.0, 0.5, 1.2, 2.5, 1.8, 0.9, 0.4, 0.1};
  hist.partials["A-A"] = {0.0, 0.2, 0.6, 1.2, 0.9, 0.4, 0.2, 0.0};
  hist.partials["A-B"] = {0.0, 0.3, 0.6, 1.3, 0.9, 0.5, 0.2, 0.1};
  return hist;
}

TEST(PlotExportManagerTests, InitialStateNotActive) {
  const PlotExportManager manager;
  EXPECT_FALSE(manager.isDialogActive());
}

TEST(PlotExportManagerTests, SavePlotFailsWithNullHistogram) {
  correlation::plotters::PlotConfig config;
  const auto result = PlotExportManager::savePlot("test.svg", "g_r", nullptr, config, {});
  EXPECT_FALSE(result.has_value());
  EXPECT_EQ(result.error(), "Histogram is null");
}

TEST(PlotExportManagerTests, SavePlotSingleHistogramSvg) {
  const auto hist = createMockHistogram("g_r");
  correlation::plotters::PlotConfig config;
  const std::filesystem::path out_path =
      std::filesystem::temp_directory_path() / "test_export_single.svg";

  if (std::filesystem::exists(out_path)) {
    std::filesystem::remove(out_path);
  }

  const auto result = PlotExportManager::savePlot(out_path.string(), "g_r", &hist, config, {});
  ASSERT_TRUE(result.has_value()) << result.error();
  EXPECT_TRUE(std::filesystem::exists(out_path));
  EXPECT_GT(std::filesystem::file_size(out_path), 0);

  std::filesystem::remove(out_path);
}

TEST(PlotExportManagerTests, SavePlotComparisonOverlaySvg) {
  const auto current_hist = createMockHistogram("g_r");
  auto pinned_hist = createMockHistogram("g_r");
  pinned_hist.bins = {0.0, 0.4, 1.0, 2.1, 1.5, 0.7, 0.3, 0.0};

  PinnedRun run1;
  run1.label = "Run 1";
  run1.histograms["g_r"] = pinned_hist;

  const std::vector<PinnedRun> pinned_runs = {run1};

  correlation::plotters::PlotConfig config;
  const std::filesystem::path out_path =
      std::filesystem::temp_directory_path() / "test_export_comparison.svg";

  if (std::filesystem::exists(out_path)) {
    std::filesystem::remove(out_path);
  }

  const auto result =
      PlotExportManager::savePlot(out_path.string(), "g_r", &current_hist, config, pinned_runs);
  ASSERT_TRUE(result.has_value()) << result.error();
  EXPECT_TRUE(std::filesystem::exists(out_path));
  EXPECT_GT(std::filesystem::file_size(out_path), 0);

  std::filesystem::remove(out_path);
}

} // namespace
} // namespace correlation::app
