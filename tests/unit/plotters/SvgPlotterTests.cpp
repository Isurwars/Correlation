// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "plotters/SvgPlotter.hpp"

#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace correlation::testing {

using namespace correlation::plotters;
using namespace correlation::analysis;

TEST(SvgPlotterTests, RendersEmptyHistogramGracefully) {
  Histogram empty_hist;
  empty_hist.title = "Empty Plot";
  empty_hist.x_label = "Distance";
  empty_hist.y_label = "g(r)";

  PlotConfig config;
  config.use_native_text = true;
  std::string svg = renderHistogramAsSvg(empty_hist, config);
  EXPECT_FALSE(svg.empty());
  EXPECT_NE(svg.find("<svg"), std::string::npos);
  EXPECT_NE(svg.find("<text"), std::string::npos);
  EXPECT_NE(svg.find("No data available"), std::string::npos);
  EXPECT_EQ(svg.find("<polyline"), std::string::npos);
}

TEST(SvgPlotterTests, RendersValidHistogramCorrectly) {
  Histogram hist;
  hist.title = "Radial Distribution Function g(r)";
  hist.x_label = "r";
  hist.y_label = "g";
  hist.x_unit = "A";
  hist.y_unit = "A^-1";

  // 5 bins
  hist.bins = {1.0, 2.0, 3.0, 4.0, 5.0};
  hist.partials["Total"] = {0.1, 0.5, 1.2, 0.8, 0.2};
  hist.partials["Si-O"] = {0.0, 0.3, 0.9, 0.4, 0.1};

  // Act
  PlotConfig config;
  config.theme = PlotConfig::Theme::Light;
  config.show_grid = true;
  config.use_native_text = true;
  std::string svg = renderHistogramAsSvg(hist, config);

  // Assert
  EXPECT_FALSE(svg.empty());
  EXPECT_NE(svg.find("<svg"), std::string::npos);
  EXPECT_NE(svg.find("</svg>"), std::string::npos);

  // Verify labels are in the SVG output
  EXPECT_NE(svg.find("<text"), std::string::npos);
  EXPECT_NE(svg.find("r (A)"), std::string::npos);
  EXPECT_NE(svg.find("g (A^-1)"), std::string::npos);

  // Verify that the structure draws polylines/lines
  EXPECT_NE(svg.find("<polyline"), std::string::npos);
  EXPECT_NE(svg.find("stroke=\"#E69F00\""), std::string::npos); // First color Orange
  EXPECT_NE(svg.find("stroke=\"#56B4E9\""), std::string::npos); // Second color Sky Blue
}

TEST(SvgPlotterTests, RendersDarkThemeCorrectly) {
  Histogram hist;
  hist.title = "Dark Theme Plot";
  hist.bins = {1.0, 2.0, 3.0};
  hist.partials["Total"] = {1.0, 2.0, 3.0};

  PlotConfig config;
  config.theme = PlotConfig::Theme::Dark;
  std::string svg = renderHistogramAsSvg(hist, config);

  // Dark theme bg is #1e1e2e
  EXPECT_NE(svg.find("fill=\"#1e1e2e\""), std::string::npos);
}

TEST(SvgPlotterTests, RendersComparisonOverlayCorrectly) {
  Histogram hist1;
  hist1.title = "Comparison";
  hist1.bins = {1.0, 2.0, 3.0};
  hist1.partials["Total"] = {0.5, 1.0, 1.5};

  Histogram hist2;
  hist2.title = "Comparison";
  hist2.bins = {1.0, 2.0, 3.0};
  hist2.partials["Total"] = {0.6, 1.1, 1.6};

  std::vector<LabeledHistogram> datasets = {{"Run A", &hist1}, {"Run B", &hist2}};

  // Act
  std::string svg = renderComparisonSvg(datasets, "Total");

  // Assert
  EXPECT_FALSE(svg.empty());
  EXPECT_NE(svg.find("<svg"), std::string::npos);
  EXPECT_NE(svg.find("<polyline"), std::string::npos);
}

TEST(SvgPlotterTests, RendersWithHoverActive) {
  Histogram hist;
  hist.title = "Hover Plot";
  hist.x_label = "r";
  hist.y_label = "g";
  hist.bins = {1.0, 2.0, 3.0};
  hist.partials["Total"] = {0.5, 1.0, 1.5};

  PlotConfig config;
  config.use_native_text = true;
  HoverInfo hover;
  hover.active = true;
  hover.mouse_x = 200.0;
  hover.mouse_y = 150.0;
  hover.widget_width = 800.0;
  hover.widget_height = 600.0;

  std::string svg = renderHistogramAsSvg(hist, config, hover);
  EXPECT_FALSE(svg.empty());
  // Should draw the dashed guide line
  EXPECT_NE(svg.find("stroke-dasharray=\"4,4\""), std::string::npos);
  // Should draw a bullet marker circle
  EXPECT_NE(svg.find("<circle"), std::string::npos);
  // Should render tooltip text
  EXPECT_NE(svg.find("<text"), std::string::npos);
}

} // namespace correlation::testing
