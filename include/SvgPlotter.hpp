/**
 * @file SvgPlotter.hpp
 * @brief Scientifically precise SVG generator for distribution function histograms.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 *
 * @details
 * Produces a self-contained SVG string from a `Histogram` object with
 * "nice" tick intervals, scientific notation, and colorblind-safe palettes.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <format>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "math/PathFont.hpp"
#include "DistributionFunctions.hpp"

namespace SvgPlotter {

/**
 * @brief Theme and layout configuration for the plot.
 */
struct PlotConfig {
  enum class Theme { Light, Dark };

  Theme theme = Theme::Light;
  double width = 1200.0;
  double height = 900.0;
  bool show_grid = true;
  bool show_markers = false;

  // Colors based on theme
  std::string bg_color() const {
    return (theme == Theme::Light) ? "#FFFFFF" : "#1e1e2e";
  }
  std::string axis_color() const {
    return (theme == Theme::Light) ? "#000000" : "#cdd6f4";
  }
  std::string grid_color() const {
    return (theme == Theme::Light) ? "#e0e0e0" : "#45475a";
  }
  std::string text_color() const {
    return (theme == Theme::Light) ? "#333333" : "#a6adc8";
  }
};

// ============================================================================
// Internal helpers
// ============================================================================
namespace detail {

/// Okabe-Ito colorblind-safe palette.
static const std::vector<std::string> kColors = {
    "#E69F00", // Orange
    "#56B4E9", // Sky Blue
    "#009E73", // Bluish Green
    "#F0E442", // Yellow
    "#0072B2", // Blue
    "#D55E00", // Vermillion
    "#CC79A7", // Reddish Purple
};

inline std::string color(std::size_t i) { return kColors[i % kColors.size()]; }

inline double mapValue(double v, double data_min, double data_max,
                       double svg_min, double svg_max) {
  if (std::abs(data_max - data_min) < 1e-15)
    return (svg_min + svg_max) / 2.0;
  return svg_min + (v - data_min) / (data_max - data_min) * (svg_max - svg_min);
}

/**
 * @brief Logic for generating "nice" human-readable tick intervals.
 */
struct NiceScale {
  double min, max, spacing;
  std::vector<double> ticks;

  NiceScale(double actual_min, double actual_max, int max_ticks = 6) {
    if (std::abs(actual_max - actual_min) < 1e-12) {
      min = actual_min - 0.5;
      max = actual_min + 0.5;
      spacing = 0.1;
    } else {
      double range = niceNum(actual_max - actual_min, false);
      spacing = niceNum(range / (max_ticks - 1), true);
      min = std::floor(actual_min / spacing) * spacing;
      max = std::ceil(actual_max / spacing) * spacing;
    }

    for (double v = min; v <= max + 0.5 * spacing; v += spacing) {
      ticks.push_back(v);
    }
  }

private:
  double niceNum(double range, bool round) {
    double exponent = std::floor(std::log10(range));
    double fraction = range / std::pow(10.0, exponent);
    double nice_fraction;

    if (round) {
      if (fraction < 1.5)
        nice_fraction = 1.0;
      else if (fraction < 3.0)
        nice_fraction = 2.0;
      else if (fraction < 7.0)
        nice_fraction = 5.0;
      else
        nice_fraction = 10.0;
    } else {
      if (fraction <= 1.0)
        nice_fraction = 1.0;
      else if (fraction <= 2.0)
        nice_fraction = 2.0;
      else if (fraction <= 5.0)
        nice_fraction = 5.0;
      else
        nice_fraction = 10.0;
    }
    return nice_fraction * std::pow(10.0, exponent);
  }
};

/**
 * @brief Formats a number for SVG display, using scientific notation if needed.
 * Returns a string that might contain <tspan> for superscripts.
 */
inline std::string fmtScientific(double v) {
  double abs_v = std::abs(v);
  if (abs_v < 1e-12)
    return "0";

  // Use scientific notation for very small/large values
  if (abs_v < 0.001 || abs_v >= 10000.0) {
    int exponent = static_cast<int>(std::floor(std::log10(abs_v)));
    double fraction = v / std::pow(10.0, exponent);
    std::string res = std::format("{:.1f}×10", fraction);
    std::string exp_s = std::to_string(exponent);
    for (char c : exp_s) {
      if (c == '-') res += "⁻";
      else if (c == '0') res += "⁰";
      else if (c == '1') res += "¹";
      else if (c == '2') res += "²";
      else if (c == '3') res += "³";
      else if (c == '4') res += "⁴";
      else if (c == '5') res += "⁵";
      else if (c == '6') res += "⁶";
      else if (c == '7') res += "⁷";
      else if (c == '8') res += "⁸";
      else if (c == '9') res += "⁹";
    }
    return res;
  }

  // Otherwise use simple formatting
  std::string s = std::format("{:.2f}", v);
  auto dot = s.find('.');
  if (dot != std::string::npos) {
    std::size_t last = s.find_last_not_of('0');
    if (last != std::string::npos && last > dot)
      s = s.substr(0, last + 1);
    else if (last == dot)
      s = s.substr(0, dot);
  }
  return s;
}

} // namespace detail

// ============================================================================
// Public API
// ============================================================================

/**
 * @brief Renders a `Histogram` as a self-contained SVG string.
 *
 * @param hist     The Histogram to render.
 * @param config   Optional plot configuration (theme, size, etc.).
 * @returns        A complete SVG document as `std::string`.
 */
inline std::string renderHistogramAsSvg(const Histogram &hist, const PlotConfig &config = {}) {
  std::string title = hist.title.empty() ? "Histogram" : hist.title;
  std::string x_label = hist.x_label.empty() ? "x" : hist.x_label;
  std::string y_label = hist.y_label.empty() ? "y" : hist.y_label;

  // Add units if available
  if (!hist.x_unit.empty()) x_label += std::format(" ({})", hist.x_unit);
  if (!hist.y_unit.empty()) {
    std::string y_unit = hist.y_unit;
    // Standardize Å^-1 to Å⁻¹ for path rendering
    if (y_unit == "Å^-1") y_unit = "Å⁻¹";
    y_label += std::format(" ({})", y_unit);
  }

  // ---- Layout configuration --------------------------------------------
  const double kW = config.width;
  const double kH = config.height;
  // Margins (enough room for labels)
  const double kLeft = 100.0;
  const double kRight = 40.0;
  const double kTop = 80.0;
  const double kBot = 90.0;

  const double px0 = kLeft;
  const double px1 = kW - kRight;
  const double py0 = kTop;
  const double py1 = kH - kBot;

  // ---- Data source -----------------------------------------------------
  const auto &partials = hist.smoothed_partials.empty() ? hist.partials : hist.smoothed_partials;
  const auto &xs = hist.bins;

  if (xs.empty() || partials.empty()) {
    return std::format(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {0} {1}\">"
        "<rect width=\"100%\" height=\"100%\" fill=\"{2}\"/>"
        "<path d=\"{3}\" fill=\"{4}\" stroke=\"none\"/></svg>",
        kW, kH, config.bg_color(), 
        PathFont::Roboto::instance().render("No data available", kW/2.0, kH/2.0 + 8.0, 24, "middle"),
        config.text_color());
  }

  // ---- Compute ranges and nice ticks -----------------------------------
  double raw_x_min = xs.front(), raw_x_max = xs.back();
  double raw_y_min = 0.0, raw_y_max = 0.0;
  for (const auto &[key, ys] : partials) {
    for (double v : ys) {
      raw_y_max = std::max(raw_y_max, v);
      raw_y_min = std::min(raw_y_min, v);
    }
  }
  // Add 5% padding to Y
  double y_padding = (raw_y_max - raw_y_min) * 0.05;
  raw_y_max += y_padding;

  detail::NiceScale xScale(raw_x_min, raw_x_max, 7);
  detail::NiceScale yScale(raw_y_min, raw_y_max, 6);

  // ---- Build SVG -------------------------------------------------------
  std::ostringstream svg;
  svg << std::format("<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\">\n", kW, kH);
  svg << std::format("  <rect width=\"100%\" height=\"100%\" fill=\"{}\" rx=\"6\"/>\n", config.bg_color());

  // Grid and Axes
  svg << "  <!-- Grid & Ticks -->\n";
  for (double yv : yScale.ticks) {
    double spy = detail::mapValue(yv, yScale.min, yScale.max, py1, py0);
    // Grid line
    if (config.show_grid) {
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"0.8\" stroke-dasharray=\"3,3\"/>\n",
                         px0, spy, px1, spy, config.grid_color());
    }
    // Tick mark
    svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                       px0 - 8.0, spy, px0, spy, config.axis_color());
    // Label
    svg << std::format("  <path d=\"{}\" fill=\"{}\" stroke=\"none\"/>\n",
                       PathFont::Roboto::instance().render(detail::fmtScientific(yv), px0 - 15.0, spy + 7.0, 20, "end"),
                       config.text_color());
  }

  for (double xv : xScale.ticks) {
    double spx = detail::mapValue(xv, xScale.min, xScale.max, px0, px1);
    if (config.show_grid) {
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"0.8\" stroke-dasharray=\"3,3\"/>\n",
                         spx, py0, spx, py1, config.grid_color());
    }
    svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                       spx, py1, spx, py1 + 8.0, config.axis_color());
    // Label
    svg << std::format("  <path d=\"{}\" fill=\"{}\" stroke=\"none\"/>\n",
                       PathFont::Roboto::instance().render(detail::fmtScientific(xv), spx, py1 + 25.0, 20, "middle"),
                       config.text_color());
  }

  // Draw axis border
  svg << std::format("  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" "
                     "fill=\"none\" stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                     px0, py0, px1 - px0, py1 - py0, config.axis_color());

  // Emphasis line for y=0 or y=1 (if histogram is S(Q))
  double focus_y = (title.find("S(Q)") != std::string::npos || title.find("S(q)") != std::string::npos) ? 1.0 : 0.0;
  if (focus_y >= yScale.min && focus_y <= yScale.max) {
    double spy = detail::mapValue(focus_y, yScale.min, yScale.max, py1, py0);
    svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"1.8\"/>\n",
                       px0, spy, px1, spy, config.axis_color());
  }

  // Polylines
  std::vector<std::pair<std::string, std::string>> legend;
  std::size_t ci = 0;
  for (const auto &[key, ys] : partials) {
    const std::string col = detail::color(ci++);
    svg << "  <polyline fill=\"none\" stroke=\"" << col
        << "\" stroke-width=\"3.0\" stroke-linejoin=\"round\" points=\"";
    std::size_t n = std::min(xs.size(), ys.size());
    for (std::size_t i = 0; i < n; ++i) {
      double sx = detail::mapValue(xs[i], xScale.min, xScale.max, px0, px1);
      double sy = detail::mapValue(ys[i], yScale.min, yScale.max, py1, py0);
      svg << std::format("{:.2f},{:.2f} ", sx, sy);
    }
    svg << "\" />\n";
    legend.push_back({key, col});
  }

  // Legend (Top Right)
  double lx = px1 - 15.0;
  double ly = py0 + 25.0;
  for (auto it = legend.rbegin(); it != legend.rend(); ++it) {
    svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"4.0\"/>\n",
                       lx - 40.0, ly, lx - 10.0, ly, it->second);
    svg << std::format("  <path d=\"{}\" fill=\"{}\" stroke=\"none\"/>\n",
                       PathFont::Roboto::instance().render(it->first, lx - 45.0, ly + 6.0, 18, "end"),
                       config.text_color());
    ly += 28.0;
  }

  // Titles/Labels
  svg << std::format("  <path d=\"{}\" fill=\"{}\" stroke=\"none\"/>\n",
                     PathFont::Roboto::instance().render(x_label, (px0 + px1) / 2.0, py1 + 75.0, 28, "middle"),
                     config.text_color());

  // Y label rotated
  svg << std::format("  <g transform=\"translate({:.1f}, {:.1f}) rotate(-90)\">\n", 40.0, (py0 + py1) / 2.0);
  svg << std::format("  <path d=\"{}\" fill=\"{}\" stroke=\"none\"/>\n",
                     PathFont::Roboto::instance().render(y_label, 0, 0, 28, "middle"),
                     config.text_color());
  svg << "  </g>\n";

  svg << std::format("  <path d=\"{}\" fill=\"{}\" stroke=\"none\"/>\n",
                     PathFont::Roboto::instance().render(title, (px0 + px1) / 2.0, py0 - 45.0, 36, "middle"),
                     config.text_color());

  svg << "</svg>\n";
  return svg.str();
}

} // namespace SvgPlotter
