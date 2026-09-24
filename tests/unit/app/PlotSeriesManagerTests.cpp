/**
 * @file PlotSeriesManagerTests.cpp
 * @brief Unit tests for PlotSeriesManager.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/PlotSeriesManager.hpp"
#include <gtest/gtest.h>
#include <optional>
#include <stdexcept>

namespace correlation::app {
namespace {

template <typename T> [[nodiscard]] T requireValue(const std::optional<T> &opt) {
  if (!opt.has_value()) {
    throw std::runtime_error("Optional value is unexpectedly empty");
  }
  return *opt;
}

TEST(PlotSeriesManagerTests, InitialStateIsPristine) {
  const PlotSeriesManager manager;
  EXPECT_EQ(manager.getPinnedRunsCount(), 0);
  EXPECT_TRUE(manager.getPinnedRuns().empty());
  EXPECT_FALSE(manager.shouldShowDifference());
  EXPECT_TRUE(manager.getCurveVisibilityMap().empty());
  EXPECT_TRUE(manager.getCustomColors().empty());
  EXPECT_TRUE(manager.getCurrentToggleKeys().empty());
}

TEST(PlotSeriesManagerTests, GenerateToggleItemsHandlesNullAndEmpty) {
  PlotSeriesManager manager;
  const correlation::plotters::PlotConfig config;

  const auto items = manager.generateToggleItems(nullptr, config);
  ASSERT_NE(items, nullptr);
  EXPECT_EQ(items->row_count(), 0);
}

TEST(PlotSeriesManagerTests, ToggleCurveVisibilityModifiesState) {
  PlotSeriesManager manager;
  const correlation::plotters::PlotConfig config;

  correlation::analysis::Histogram hist;
  hist.partials["Total"] = {1.0, 2.0};
  hist.partials["Si-O"] = {0.5, 1.0};

  const auto items = manager.generateToggleItems(&hist, config);
  ASSERT_NE(items, nullptr);
  ASSERT_EQ(items->row_count(), 2);
  EXPECT_EQ(manager.getCurrentToggleKeys().size(), 2);

  // Default is visible
  EXPECT_TRUE(manager.isCurveVisible("Total"));
  EXPECT_TRUE(manager.isCurveVisible("Si-O"));

  // Toggle Total (id 0) to false
  manager.setCurveVisible(0, false);
  EXPECT_FALSE(manager.isCurveVisible("Total"));
  EXPECT_TRUE(manager.isCurveVisible("Si-O"));

  // Toggle all curves to false
  manager.setAllCurvesVisible(false);
  EXPECT_FALSE(manager.isCurveVisible("Total"));
  EXPECT_FALSE(manager.isCurveVisible("Si-O"));

  // Toggle all curves to true
  manager.setAllCurvesVisible(true);
  EXPECT_TRUE(manager.isCurveVisible("Total"));
  EXPECT_TRUE(manager.isCurveVisible("Si-O"));
}

TEST(PlotSeriesManagerTests, CustomColorAssignment) {
  PlotSeriesManager manager;
  const correlation::plotters::PlotConfig config;

  correlation::analysis::Histogram hist;
  hist.partials["Total"] = {1.0};

  ASSERT_NE(manager.generateToggleItems(&hist, config), nullptr);
  manager.setCustomColor(0, "#FF5500");

  const auto &colors = manager.getCustomColors();
  ASSERT_TRUE(colors.contains("Total"));
  EXPECT_EQ(colors.at("Total"), "#FF5500");

  // Regenerate toggle items should reflect custom color
  const auto items = manager.generateToggleItems(&hist, config);
  ASSERT_EQ(items->row_count(), 1);
  const auto item = requireValue(items->row_data(0));
  EXPECT_EQ(item.color_hex, slint::Color::from_rgb_uint8(0xFF, 0x55, 0x00));
}

TEST(PlotSeriesManagerTests, PinnedRunsLifecycle) {
  PlotSeriesManager manager;

  std::map<std::string, correlation::analysis::Histogram> run1;
  correlation::analysis::Histogram hist1;
  hist1.partials["Total"] = {1.0};
  run1["rdf"] = hist1;

  manager.pinCurrentRun(run1);
  EXPECT_EQ(manager.getPinnedRunsCount(), 1);
  ASSERT_EQ(manager.getPinnedRuns().size(), 1);
  EXPECT_EQ(manager.getPinnedRuns()[0].label, "Run 1");

  std::map<std::string, correlation::analysis::Histogram> run2;
  run2["rdf"] = hist1;
  manager.pinCurrentRun(run2);
  EXPECT_EQ(manager.getPinnedRunsCount(), 2);
  EXPECT_EQ(manager.getPinnedRuns()[1].label, "Run 2");

  manager.clearPinnedRuns();
  EXPECT_EQ(manager.getPinnedRunsCount(), 0);
  EXPECT_TRUE(manager.getPinnedRuns().empty());
}

TEST(PlotSeriesManagerTests, DifferenceCurveFlag) {
  PlotSeriesManager manager;
  EXPECT_FALSE(manager.shouldShowDifference());

  manager.setShowDifference(true);
  EXPECT_TRUE(manager.shouldShowDifference());

  manager.setShowDifference(false);
  EXPECT_FALSE(manager.shouldShowDifference());
}

TEST(PlotSeriesManagerTests, ResetRestoresDefaults) {
  PlotSeriesManager manager;
  const correlation::plotters::PlotConfig config;

  correlation::analysis::Histogram hist;
  hist.partials["Total"] = {1.0};
  ASSERT_NE(manager.generateToggleItems(&hist, config), nullptr);

  manager.setCustomColor(0, "#123456");
  manager.setCurveVisible(0, false);
  manager.setShowDifference(true);

  std::map<std::string, correlation::analysis::Histogram> run;
  run["test"] = hist;
  manager.pinCurrentRun(run);

  EXPECT_EQ(manager.getPinnedRunsCount(), 1);
  EXPECT_TRUE(manager.shouldShowDifference());
  EXPECT_FALSE(manager.getCustomColors().empty());
  EXPECT_FALSE(manager.getCurveVisibilityMap().empty());

  manager.reset();

  EXPECT_EQ(manager.getPinnedRunsCount(), 0);
  EXPECT_FALSE(manager.shouldShowDifference());
  EXPECT_TRUE(manager.getCustomColors().empty());
  EXPECT_TRUE(manager.getCurveVisibilityMap().empty());
  EXPECT_TRUE(manager.getCurrentToggleKeys().empty());
}

TEST(PlotSeriesManagerTests, PartialVisibilityRankThreshold) {
  PlotSeriesManager manager;
  const correlation::plotters::PlotConfig config;

  correlation::analysis::Histogram hist;
  for (int i = 0; i < 9; ++i) {
    hist.partials["P" + std::to_string(i)] = {1.0};
  }

  const auto items = manager.generateToggleItems(&hist, config);
  ASSERT_NE(items, nullptr);
  ASSERT_EQ(items->row_count(), 9);

  // First 7 sorted partials (P0..P6) are visible by default
  for (int i = 0; i < 7; ++i) {
    EXPECT_TRUE(requireValue(items->row_data(i)).visible);
  }
  // 8th and 9th (P7, P8) default to hidden (rank >= 7)
  EXPECT_FALSE(requireValue(items->row_data(7)).visible);
  EXPECT_FALSE(requireValue(items->row_data(8)).visible);
}

} // namespace
} // namespace correlation::app
