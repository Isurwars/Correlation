// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
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

TEST(SvgPlotterTests, RendersWithHover2DNearestSnapping) {
  Histogram hist;
  hist.title = "Hover Plot 2D";
  hist.x_label = "r";
  hist.y_label = "g";
  hist.bins = {1.0, 2.0, 3.0};
  hist.partials["Total"] = {10.0, 10.0, 10.0};
  hist.partials["Si-O"] = {1.0, 1.0, 1.0};

  PlotConfig config;
  config.use_native_text = true;

  // xScale maps [1.0, 3.0] to [100.0, 1160.0] (x=2.0 -> 630.0)
  // yScale maps [0.0, 10.5] to [810.0, 50.0] (y=1.0 -> 737.6)
  HoverInfo hover;
  hover.active = true;
  hover.mouse_x = 630.0;
  hover.mouse_y = 737.6;
  hover.widget_width = 1200.0;
  hover.widget_height = 900.0;

  std::string svg = renderHistogramAsSvg(hist, config, hover);
  EXPECT_FALSE(svg.empty());

  // Verify that "Si-O" is identified as the nearest curve
  EXPECT_NE(svg.find("Si-O: 1.0000 (nearest)"), std::string::npos);
}

TEST(SvgPlotterTests, RendersShadedCurveCorrectly) {
  Histogram hist;
  hist.title = "Shaded Plot";
  hist.bins = {1.0, 2.0, 3.0};
  hist.partials["Total"] = {1.0, 2.0, 1.5};

  PlotConfig config;
  config.fill_area = true;
  config.use_native_text = true;

  std::string svg = renderHistogramAsSvg(hist, config);
  EXPECT_FALSE(svg.empty());

  // Verify that linearGradient and polygon elements are in the SVG output
  EXPECT_NE(svg.find("<linearGradient id=\"area-grad-0\""), std::string::npos);
  EXPECT_NE(svg.find("<polygon fill=\"url(#area-grad-0)\""), std::string::npos);
}

TEST(SvgPlotterTests, RendersContinuousColormapsAndCustomColors) {
  Histogram hist;
  hist.title = "Colormap Test";
  hist.bins = {1.0, 2.0, 3.0};
  hist.partials["Total"] = {1.0, 2.0, 1.5};
  hist.partials["O-O"] = {0.5, 0.8, 1.0};

  // Test Magma continuous palette
  PlotConfig magma_config;
  magma_config.palette = PlotConfig::Palette::Magma;
  std::string magma_svg = renderHistogramAsSvg(hist, magma_config);
  EXPECT_FALSE(magma_svg.empty());
  EXPECT_NE(magma_svg.find("stroke=\"#"), std::string::npos);

  // Test custom curve color override
  PlotConfig custom_config;
  std::map<std::string, std::string> custom_colors = {{"Total", "#FF00FF"}};
  std::string custom_svg = renderHistogramAsSvg(hist, custom_config, {}, {}, {}, custom_colors);
  EXPECT_NE(custom_svg.find("stroke=\"#FF00FF\""), std::string::npos);
}

TEST(SvgPlotterTests, RendersCustomMarkerSizeCorrectly) {
  Histogram hist;
  hist.title = "Marker Test";
  hist.bins = {1.0, 2.0, 3.0};
  hist.partials["Total"] = {1.0, 2.0, 1.5};

  PlotConfig config;
  config.show_markers = true;
  config.marker_size = 5.2;

  std::string svg = renderHistogramAsSvg(hist, config);
  EXPECT_FALSE(svg.empty());
  EXPECT_NE(svg.find("r=\"5.2\""), std::string::npos);
}

TEST(SvgPlotterTests, AppliesManualZoomBoundsAndClipping) {
  Histogram hist;
  hist.title = "Zoom Test";
  hist.bins = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  hist.partials["Total"] = {0.0, 1.0, 2.5, 3.0, 1.5, 0.2};

  PlotConfig config;
  config.manual_x_min = 1.5;
  config.manual_x_max = 3.5;
  config.manual_y_min = 0.5;
  config.manual_y_max = 2.8;

  std::string svg = renderHistogramAsSvg(hist, config);
  EXPECT_FALSE(svg.empty());
  EXPECT_NE(svg.find("id=\"plot-area-clip\""), std::string::npos);
  EXPECT_NE(svg.find("clip-path=\"url(#plot-area-clip)\""), std::string::npos);
}

TEST(SvgPlotterTests, RendersReferenceLinesWithLabels) {
  Histogram hist;
  hist.title = "Marker Lines Test";
  hist.bins = {1.0, 2.0, 3.0, 4.0};
  hist.partials["Total"] = {0.5, 1.0, 1.5, 2.0};

  PlotConfig config;
  config.use_native_text = true;
  config.reference_lines.push_back(ReferenceLine{
      .value = 2.5,
      .is_vertical = true,
      .label = "X: 2.500",
      .color_hex = "#0284C7",
  });
  config.reference_lines.push_back(ReferenceLine{
      .value = 1.2,
      .is_vertical = false,
      .label = "Y: 1.200",
      .color_hex = "#E11D48",
  });

  std::string svg = renderHistogramAsSvg(hist, config);
  EXPECT_FALSE(svg.empty());
  EXPECT_NE(svg.find("stroke=\"#0284C7\""), std::string::npos);
  EXPECT_NE(svg.find("stroke=\"#E11D48\""), std::string::npos);
  EXPECT_NE(svg.find("stroke-dasharray=\"4,4\""), std::string::npos);
  EXPECT_NE(svg.find("X: 2.500"), std::string::npos);
  EXPECT_NE(svg.find("Y: 1.200"), std::string::npos);
}

TEST(SvgPlotterTests, ScreenToDataGeometryTransformsAccurately) {
  PlotConfig config;
  config.width = 1000.0;
  config.height = 600.0;

  auto geom = detail::getViewportGeometry(config);
  EXPECT_DOUBLE_EQ(geom.kw, 1000.0);
  EXPECT_DOUBLE_EQ(geom.kh, 600.0);
  EXPECT_DOUBLE_EQ(geom.px0, 100.0);
  EXPECT_DOUBLE_EQ(geom.px1, 960.0);
  EXPECT_DOUBLE_EQ(geom.py0, 50.0);
  EXPECT_DOUBLE_EQ(geom.py1, 530.0);

  detail::NiceScale xs(detail::DataRange{.min = 0.0, .max = 10.0}, 10, true);
  detail::NiceScale ys(detail::DataRange{.min = 0.0, .max = 5.0}, 5, true);

  // Screen matches SVG 1:1 when widget is 1000x600
  auto [data_x, data_y] = detail::screenToData(100.0, 530.0, 1000.0, 600.0, config, xs, ys);
  EXPECT_NEAR(data_x, 0.0, 1e-4);
  EXPECT_NEAR(data_y, 0.0, 1e-4);

  auto [top_x, top_y] = detail::screenToData(960.0, 50.0, 1000.0, 600.0, config, xs, ys);
  EXPECT_NEAR(top_x, 10.0, 1e-4);
  EXPECT_NEAR(top_y, 5.0, 1e-4);
}

} // namespace correlation::testing
