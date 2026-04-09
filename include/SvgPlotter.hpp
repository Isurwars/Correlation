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

#include "DistributionFunctions.hpp"

namespace SvgPlotter {

/**
 * @brief Theme and layout configuration for the plot.
 */
struct PlotConfig {
  enum class Theme { Light, Dark };

  Theme theme = Theme::Light;
  double width = 800.0;
  double height = 450.0;
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
    return std::format("{:.1f}×10<tspan dy=\"-0.6em\" font-size=\"75%\">{}</tspan><tspan dy=\"0.6em\"> </tspan>",
                       fraction, exponent);
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
  if (!hist.y_unit.empty()) y_label += std::format(" ({})", hist.y_unit);

  // ---- Layout configuration --------------------------------------------
  const double kW = config.width;
  const double kH = config.height;
  // Margins (enough room for labels)
  const double kLeft = 75.0;
  const double kRight = 25.0;
  const double kTop = 50.0;
  const double kBot = 60.0;

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
        "<text x=\"{3}\" y=\"{4}\" fill=\"{5}\" font-size=\"18\" "
        "text-anchor=\"middle\">No data available</text></svg>",
        kW, kH, config.bg_color(), kW / 2.0, kH / 2.0, config.text_color());
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
  svg << std::format("<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\" font-family=\"sans-serif\">\n", kW, kH);
  svg << std::format("  <rect width=\"100%\" height=\"100%\" fill=\"{}\" rx=\"4\"/>\n", config.bg_color());

  // Grid and Axes
  svg << "  <!-- Grid & Ticks -->\n";
  for (double yv : yScale.ticks) {
    double spy = detail::mapValue(yv, yScale.min, yScale.max, py1, py0);
    // Grid line
    if (config.show_grid) {
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"0.5\" stroke-dasharray=\"2,2\"/>\n",
                         px0, spy, px1, spy, config.grid_color());
    }
    // Tick mark
    svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"1\"/>\n",
                       px0 - 5.0, spy, px0, spy, config.axis_color());
    // Label
    svg << std::format("  <text x=\"{:.1f}\" y=\"{:.1f}\" fill=\"{}\" font-size=\"11\" "
                       "text-anchor=\"end\" dominant-baseline=\"middle\">{}</text>\n",
                       px0 - 8.0, spy, config.text_color(), detail::fmtScientific(yv));
  }

  for (double xv : xScale.ticks) {
    double spx = detail::mapValue(xv, xScale.min, xScale.max, px0, px1);
    if (config.show_grid) {
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"0.5\" stroke-dasharray=\"2,2\"/>\n",
                         spx, py0, spx, py1, config.grid_color());
    }
    svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"1\"/>\n",
                       spx, py1, spx, py1 + 5.0, config.axis_color());
    svg << std::format("  <text x=\"{:.1f}\" y=\"{:.1f}\" fill=\"{}\" font-size=\"11\" "
                       "text-anchor=\"middle\">{}</text>\n",
                       spx, py1 + 18.0, config.text_color(), detail::fmtScientific(xv));
  }

  // Draw axis border
  svg << std::format("  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" "
                     "fill=\"none\" stroke=\"{}\" stroke-width=\"1\"/>\n",
                     px0, py0, px1 - px0, py1 - py0, config.axis_color());

  // Emphasis line for y=0 or y=1 (if histogram is S(Q))
  double focus_y = (title.find("S(Q)") != std::string::npos || title.find("S(q)") != std::string::npos) ? 1.0 : 0.0;
  if (focus_y >= yScale.min && focus_y <= yScale.max) {
    double spy = detail::mapValue(focus_y, yScale.min, yScale.max, py1, py0);
    svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"1.2\"/>\n",
                       px0, spy, px1, spy, config.axis_color());
  }

  // Polylines
  std::vector<std::pair<std::string, std::string>> legend;
  std::size_t ci = 0;
  for (const auto &[key, ys] : partials) {
    const std::string col = detail::color(ci++);
    svg << "  <polyline fill=\"none\" stroke=\"" << col
        << "\" stroke-width=\"2.0\" stroke-linejoin=\"round\" points=\"";
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
  double lx = px1 - 10.0;
  double ly = py0 + 15.0;
  for (auto it = legend.rbegin(); it != legend.rend(); ++it) {
    svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"2.5\"/>\n",
                       lx - 25.0, ly, lx - 5.0, ly, it->second);
    svg << std::format("  <text x=\"{:.1f}\" y=\"{:.1f}\" fill=\"{}\" font-size=\"12\" "
                       "text-anchor=\"end\" dominant-baseline=\"middle\">{}</text>\n",
                       lx - 30.0, ly, config.text_color(), it->first);
    ly += 18.0;
  }

  // Titles/Labels
  svg << std::format("  <text x=\"{:.1f}\" y=\"{:.1f}\" fill=\"{}\" font-size=\"14\" "
                     "text-anchor=\"middle\" font-weight=\"bold\">{}</text>\n",
                     (px0 + px1) / 2.0, py1 + 45.0, config.axis_color(), x_label);

  svg << std::format("  <text transform=\"rotate(-90)\" x=\"{:.1f}\" y=\"{:.1f}\" "
                     "fill=\"{}\" font-size=\"14\" text-anchor=\"middle\" font-weight=\"bold\">{}</text>\n",
                     -((py0 + py1) / 2.0), 16.0, config.axis_color(), y_label);

  svg << std::format("  <text x=\"{:.1f}\" y=\"{:.1f}\" fill=\"{}\" font-size=\"18\" "
                     "text-anchor=\"middle\" font-weight=\"bold\">{}</text>\n",
                     (px0 + px1) / 2.0, py0 - 15.0, config.axis_color(), title);

  svg << "</svg>\n";
  return svg.str();
}

} // namespace SvgPlotter
