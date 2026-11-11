/**
 * @file PlotSeriesManagerTests.cpp
 * @brief Unit tests for PlotSeriesManager.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/formatters/PlotSeriesManager.hpp"
#include <gtest/gtest.h>
#include <optional>
#include <stdexcept>

using correlation::real_t;

namespace correlation::app {
namespace {

template <typename T> [[nodiscard]] T requireValue(const std::optional<T> &opt) {
  if (!opt.has_value()) {
    throw std::runtime_error("Optional value is unexpectedly empty");
  }
  return *opt;
}

void verifyRankThresholdVisibilities(const std::shared_ptr<slint::Model<CurveToggleItem>> &items) {
  ASSERT_EQ(items->row_count(), 12);
  for (int i = 0; i < 10; ++i) {
    EXPECT_TRUE(requireValue(items->row_data(i)).visible);
  }
  EXPECT_FALSE(requireValue(items->row_data(10)).visible);
  EXPECT_FALSE(requireValue(items->row_data(11)).visible);
}

void verifyFeFilterMatches(const std::shared_ptr<slint::Model<CurveToggleItem>> &items) {
  ASSERT_EQ(items->row_count(), 4);
  EXPECT_FALSE(requireValue(items->row_data(0)).matches_filter);
  EXPECT_TRUE(requireValue(items->row_data(1)).matches_filter);
  EXPECT_TRUE(requireValue(items->row_data(2)).matches_filter);
  EXPECT_FALSE(requireValue(items->row_data(3)).matches_filter);
}

void verifyAllMatchFilter(const std::shared_ptr<slint::Model<CurveToggleItem>> &items) {
  for (std::size_t i = 0; i < items->row_count(); ++i) {
    EXPECT_TRUE(requireValue(items->row_data(i)).matches_filter);
  }
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
  hist.partials["O-O"] = {0.2, 0.4};

  const auto items = manager.generateToggleItems(&hist, config);
  ASSERT_NE(items, nullptr);
  ASSERT_EQ(items->row_count(), 3);
  EXPECT_EQ(manager.getCurrentToggleKeys().size(), 3);

  // Default is visible
  EXPECT_TRUE(manager.isCurveVisible("Total"));
  EXPECT_TRUE(manager.isCurveVisible("Si-O"));
  EXPECT_TRUE(manager.isCurveVisible("O-O"));

  // Toggle Total (id 0) to false
  manager.setCurveVisible(0, false);
  EXPECT_FALSE(manager.isCurveVisible("Total"));
  EXPECT_TRUE(manager.isCurveVisible("Si-O"));

  // Toggle all curves to false
  manager.setAllCurvesVisible(false);
  EXPECT_FALSE(manager.isCurveVisible("Total"));
  EXPECT_FALSE(manager.isCurveVisible("Si-O"));
  EXPECT_FALSE(manager.isCurveVisible("O-O"));

  // Toggle all curves to true
  manager.setAllCurvesVisible(true);
  EXPECT_TRUE(manager.isCurveVisible("Total"));
  EXPECT_TRUE(manager.isCurveVisible("Si-O"));
  EXPECT_TRUE(manager.isCurveVisible("O-O"));
}

TEST(PlotSeriesManagerTests, SinglePartialOmitsRedundantCurveAndReportsOnlyTotal) {
  PlotSeriesManager manager;
  const correlation::plotters::PlotConfig config;

  correlation::analysis::Histogram hist;
  hist.partials["Total"] = {1.0, 2.0};
  hist.partials["Si-Si"] = {1.0, 2.0};

  const auto items = manager.generateToggleItems(&hist, config);
  ASSERT_NE(items, nullptr);
  ASSERT_EQ(items->row_count(), 1);
  EXPECT_EQ(manager.getCurrentToggleKeys().size(), 1);
  EXPECT_EQ(std::string(requireValue(items->row_data(0)).label.data()), "Total");
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
  for (int i = 0; i < 12; ++i) {
    hist.partials["P" + (i < 10 ? std::string("0") : std::string("")) + std::to_string(i)] = {1.0};
  }

  const auto items = manager.generateToggleItems(&hist, config);
  ASSERT_NE(items, nullptr);
  verifyRankThresholdVisibilities(items);
}

TEST(PlotSeriesManagerTests, ToggleItemsMatchSvgOrderAndColors) {
  PlotSeriesManager manager;
  const correlation::plotters::PlotConfig config;

  correlation::analysis::Histogram hist;
  hist.bins = {1.0, 2.0, 3.0, 4.0, 5.0};
  hist.partials["Total"] = {0.1, 0.5, 1.2, 0.8, 0.2};
  hist.partials["Si-O"] = {0.0, 0.3, 0.9, 0.4, 0.1};
  hist.partials["O-O"] = {0.0, 0.2, 0.3, 0.4, 0.1};

  const auto items = manager.generateToggleItems(&hist, config);
  ASSERT_NE(items, nullptr);
  ASSERT_EQ(items->row_count(), 3);

  const auto item0 = requireValue(items->row_data(0));
  EXPECT_EQ(std::string(item0.label.data()), "Total");
  EXPECT_EQ(item0.color_hex, slint::Color::from_rgb_uint8(0x00, 0x00, 0x00));

  const auto item1 = requireValue(items->row_data(1));
  EXPECT_EQ(std::string(item1.label.data()), "Si-O");
  EXPECT_EQ(item1.color_hex, slint::Color::from_rgb_uint8(0xE6, 0x9F, 0x00));
}

TEST(PlotSeriesManagerTests, AshcroftWeightsRankPartialsDescending) {
  PlotSeriesManager manager;
  const correlation::plotters::PlotConfig config;

  correlation::analysis::Histogram hist;
  hist.partials["Total"] = {1.0};
  hist.partials["A-A"] = {10.0};
  hist.partials["B-B"] = {1.0};
  hist.partials["C-C"] = {5.0};

  const std::map<std::string, real_t> weights = {
      {"A-A", 0.1},
      {"B-B", 0.8},
      {"C-C", 0.5},
  };

  const auto items = manager.generateToggleItems(&hist, config, weights);
  ASSERT_NE(items, nullptr);
  ASSERT_EQ(items->row_count(), 4);

  // Total anchored on top
  EXPECT_EQ(std::string(requireValue(items->row_data(0)).label.data()), "Total");
  // B-B has highest Ashcroft weight (0.8)
  EXPECT_EQ(std::string(requireValue(items->row_data(1)).label.data()), "B-B");
  // C-C has second highest Ashcroft weight (0.5)
  EXPECT_EQ(std::string(requireValue(items->row_data(2)).label.data()), "C-C");
  // A-A has lowest Ashcroft weight (0.1)
  EXPECT_EQ(std::string(requireValue(items->row_data(3)).label.data()), "A-A");
}

TEST(PlotSeriesManagerTests, ColorsRemainStableAcrossVisibilityToggles) {
  PlotSeriesManager manager;
  const correlation::plotters::PlotConfig config;

  correlation::analysis::Histogram hist;
  hist.partials["Total"] = {1.0};
  hist.partials["Si-Si"] = {2.0};
  hist.partials["Si-O"] = {1.5};
  hist.partials["O-O"] = {0.8};

  [[maybe_unused]] const auto items_before = manager.generateToggleItems(&hist, config);
  const auto assigned_colors_before = manager.getAssignedCurveColors();

  // Hide Si-Si (toggle off)
  manager.setCurveVisible(1, false);
  [[maybe_unused]] const auto items_after = manager.generateToggleItems(&hist, config);
  const auto assigned_colors_after = manager.getAssignedCurveColors();

  // Color assignments must be completely identical across toggles
  EXPECT_EQ(assigned_colors_before, assigned_colors_after);
  EXPECT_EQ(assigned_colors_after.at("Total"), "#000000");
  EXPECT_FALSE(assigned_colors_after.at("Si-Si").empty());
  EXPECT_FALSE(assigned_colors_after.at("Si-O").empty());
  EXPECT_FALSE(assigned_colors_after.at("O-O").empty());
}

TEST(PlotSeriesManagerTests, FilterCurvesSubstringCaseInsensitive) {
  PlotSeriesManager manager;
  const correlation::plotters::PlotConfig config;

  correlation::analysis::Histogram hist;
  hist.partials["Total"] = {1.0};
  hist.partials["Fe-Fe"] = {2.0};
  hist.partials["Fe-Si"] = {1.5};
  hist.partials["Si-Si"] = {0.8};

  manager.setFilterText("fe");
  EXPECT_EQ(manager.getFilterText(), "fe");

  const auto items = manager.generateToggleItems(&hist, config);
  ASSERT_NE(items, nullptr);
  verifyFeFilterMatches(items);

  // Clear filter
  manager.setFilterText("");
  const auto items_cleared = manager.generateToggleItems(&hist, config);
  verifyAllMatchFilter(items_cleared);
}

TEST(PlotSeriesManagerTests, InvertCurveVisibilityGlobalAndFilterScoped) {
  PlotSeriesManager manager;
  const correlation::plotters::PlotConfig config;

  correlation::analysis::Histogram hist;
  hist.partials["Total"] = {1.0};
  hist.partials["Fe-Fe"] = {2.0};
  hist.partials["Si-Si"] = {1.5};

  [[maybe_unused]] const auto items = manager.generateToggleItems(&hist, config);
  EXPECT_TRUE(manager.isCurveVisible("Total"));
  EXPECT_TRUE(manager.isCurveVisible("Fe-Fe"));
  EXPECT_TRUE(manager.isCurveVisible("Si-Si"));

  // Global Invert: all true become false
  manager.invertCurveVisibility();
  EXPECT_FALSE(manager.isCurveVisible("Total"));
  EXPECT_FALSE(manager.isCurveVisible("Fe-Fe"));
  EXPECT_FALSE(manager.isCurveVisible("Si-Si"));

  // Turn Fe-Fe on
  manager.setCurveVisible(1, true);
  EXPECT_TRUE(manager.isCurveVisible("Fe-Fe"));

  // Filter-scoped Invert on "Fe": only Fe-Fe inverts (to false)
  manager.setFilterText("Fe");
  manager.invertCurveVisibility();
  EXPECT_FALSE(manager.isCurveVisible("Fe-Fe"));
  EXPECT_FALSE(manager.isCurveVisible("Si-Si")); // remains false
  EXPECT_FALSE(manager.isCurveVisible("Total")); // remains false

  // Filter-scoped setAllCurvesVisible(true) on "Fe": only Fe-Fe turns true
  manager.setAllCurvesVisible(true);
  EXPECT_TRUE(manager.isCurveVisible("Fe-Fe"));
  EXPECT_FALSE(manager.isCurveVisible("Si-Si")); // untouched
}

TEST(PlotSeriesManagerTests, DualViewPreservesStableColorsAcrossFrames) {
  PlotSeriesManager manager;
  const correlation::plotters::PlotConfig config;

  correlation::analysis::Histogram frame1_hist;
  frame1_hist.partials["Ensemble Mean"] = {1.0, 1.5, 2.0};
  frame1_hist.partials["Frame 1"] = {100.0, 200.0, 300.0};

  const auto items1 = manager.generateToggleItems(&frame1_hist, config);
  ASSERT_NE(items1, nullptr);
  ASSERT_EQ(items1->row_count(), 2);
  const auto assigned1 = manager.getAssignedCurveColors();
  EXPECT_EQ(assigned1.at("Ensemble Mean"),
            correlation::plotters::detail::color(0, 2, config.palette));
  EXPECT_EQ(assigned1.at("Frame 1"), correlation::plotters::detail::color(1, 2, config.palette));

  // Advance to Frame 2 where frame values are lower than mean
  correlation::analysis::Histogram frame2_hist;
  frame2_hist.partials["Ensemble Mean"] = {1.0, 1.5, 2.0};
  frame2_hist.partials["Frame 2"] = {0.01, 0.02, 0.01};

  const auto items2 = manager.generateToggleItems(&frame2_hist, config);
  ASSERT_NE(items2, nullptr);
  ASSERT_EQ(items2->row_count(), 2);
  const auto assigned2 = manager.getAssignedCurveColors();
  EXPECT_EQ(assigned2.at("Ensemble Mean"),
            correlation::plotters::detail::color(0, 2, config.palette));
  EXPECT_EQ(assigned2.at("Frame 2"), correlation::plotters::detail::color(1, 2, config.palette));
}

TEST(PlotSeriesManagerTests, DualViewCustomColorPersistsAcrossFrames) {
  PlotSeriesManager manager;
  const correlation::plotters::PlotConfig config;

  correlation::analysis::Histogram frame1_hist;
  frame1_hist.partials["Ensemble Mean"] = {1.0, 1.5, 2.0};
  frame1_hist.partials["Frame 1"] = {1.2, 1.4, 2.1};

  [[maybe_unused]] const auto items1 = manager.generateToggleItems(&frame1_hist, config);
  // Set custom color on Frame 1 (id 1)
  manager.setCustomColor(1, "#FF5733");

  const auto assigned1 = manager.getAssignedCurveColors();
  EXPECT_EQ(assigned1.at("Frame 1"), "#FF5733");

  // Advance to Frame 2
  correlation::analysis::Histogram frame2_hist;
  frame2_hist.partials["Ensemble Mean"] = {1.0, 1.5, 2.0};
  frame2_hist.partials["Frame 2"] = {0.9, 1.3, 1.8};

  [[maybe_unused]] const auto items2 = manager.generateToggleItems(&frame2_hist, config);
  const auto assigned2 = manager.getAssignedCurveColors();
  EXPECT_EQ(assigned2.at("Frame 2"), "#FF5733");
  EXPECT_EQ(assigned2.at("Ensemble Mean"),
            correlation::plotters::detail::color(0, 2, config.palette));
}

} // namespace
} // namespace correlation::app
