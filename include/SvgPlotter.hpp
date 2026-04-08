/**
 * @file SvgPlotter.hpp
 * @brief Lightweight SVG generator for distribution function histograms.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 *
 * @details
 * Produces a self-contained SVG string from a `Histogram` object that can be
 * loaded directly by Slint via `slint::Image::load_from_svg_data()`.
 *
 * **Design goals:**
 * - Zero external dependencies beyond the C++ standard library.
 * - Renders all partial keys present in `smoothed_partials` (falling back to
 *   `partials` if smoothing was not applied).
 * - Uses a fixed 800×450 coordinate system with 70 px left / 20 px right /
 *   50 px top / 50 px bottom margins to accommodate axis labels.
 * - Colour-cycles through a fixed Material Design-inspired palette for up to
 *   8 partial keys.
 *
 * **Matplotlib fallback:**
 * If richer visualisation is ever needed (e.g. publication figures), replace
 * the call to `renderHistogramAsSvg` with a call to a function that invokes
 * Python/Matplotlib via `popen()`, writes a PNG to a temp path, and returns
 * the path to be loaded via `slint::Image::load_from_path()`.
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

// ============================================================================
// Internal helpers
// ============================================================================
namespace detail {

/// Material Design-inspired colour palette (hex strings, inline SVG
/// fill/stroke).
static const std::vector<std::string> kColors = {
    "#4FC3F7", // light-blue-300
    "#EF9A9A", // red-200
    "#A5D6A7", // green-200
    "#FFE082", // amber-200
    "#CE93D8", // purple-200
    "#80DEEA", // cyan-200
    "#FFAB91", // deep-orange-200
    "#B0BEC5", // blue-grey-200
};

/// Returns the colour string for the i-th partial (wraps around).
inline std::string color(std::size_t i) { return kColors[i % kColors.size()]; }

/// Map a value in [data_min, data_max] to SVG coordinate space [svg_min,
/// svg_max].
inline double mapValue(double v, double data_min, double data_max,
                       double svg_min, double svg_max) {
  if (std::abs(data_max - data_min) < 1e-15)
    return (svg_min + svg_max) / 2.0;
  return svg_min + (v - data_min) / (data_max - data_min) * (svg_max - svg_min);
}

/// Append a @c \<polyline\> for the given x/y vectors.
inline void appendPolyline(std::ostringstream &svg,
                           const std::vector<double> &xs,
                           const std::vector<double> &ys, double x_min,
                           double x_max, double y_min, double y_max,
                           double plot_x0, double plot_x1, double plot_y0,
                           double plot_y1, const std::string &stroke,
                           const std::string &label) {
  svg << "  <polyline fill=\"none\" stroke=\"" << stroke
      << "\" stroke-width=\"1.8\" stroke-linejoin=\"round\" points=\"";

  std::size_t n = std::min(xs.size(), ys.size());
  for (std::size_t i = 0; i < n; ++i) {
    double px = mapValue(xs[i], x_min, x_max, plot_x0, plot_x1);
    // SVG y-axis is inverted: larger data values map to smaller y coords.
    double py = mapValue(ys[i], y_min, y_max, plot_y1, plot_y0);
    svg << std::format("{:.2f},{:.2f} ", px, py);
  }
  svg << "\" />\n";

  // Legend entry (drawn later by the caller)
  (void)label; // stored separately
}

/// Format a double tick label, stripping unnecessary trailing zeros.
inline std::string fmtTick(double v) {
  if (std::abs(v) < 1e-9)
    return "0";
  std::string s = std::format("{:.2f}", v);
  // Strip trailing zeros after decimal point.
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
 * All partial keys found in `hist.smoothed_partials` are rendered as
 * individual coloured polylines with a legend.  If `smoothed_partials` is
 * empty, the function falls back to `hist.partials`.
 *
 * @param hist     The Histogram to render.
 * @param title    Title text displayed above the plot.
 * @param x_label  Label for the x-axis (e.g., "r (Å)").
 * @param y_label  Label for the y-axis (e.g., "g(r)").
 * @returns        A complete SVG document as `std::string`.
 */
inline std::string renderHistogramAsSvg(const Histogram &hist) {
  // Use metadata if available, otherwise fallback
  std::string title = hist.title.empty() ? "Histogram" : hist.title;
  std::string x_label = hist.x_label.empty() ? "x" : hist.x_label;
  std::string y_label = hist.y_label.empty() ? "y" : hist.y_label;
  // ---- Layout constants ------------------------------------------------
  constexpr double kW = 800.0;
  constexpr double kH = 450.0;
  constexpr double kLeft = 72.0;
  constexpr double kRight = 20.0;
  constexpr double kTop = 40.0;
  constexpr double kBot = 55.0;

  const double px0 = kLeft;
  const double px1 = kW - kRight;
  const double py0 = kTop;
  const double py1 = kH - kBot;

  // ---- Choose data source (smoothed preferred) -------------------------
  const auto &partials =
      hist.smoothed_partials.empty() ? hist.partials : hist.smoothed_partials;
  const auto &xs = hist.bins;

  if (xs.empty() || partials.empty()) {
    // Return a minimal "no data" SVG.
    return std::format(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        "viewBox=\"0 0 {} {}\">"
        "<rect width=\"100%\" height=\"100%\" fill=\"#1e1e2e\"/>"
        "<text x=\"{}\" y=\"{}\" fill=\"#cdd6f4\" font-size=\"18\" "
        "text-anchor=\"middle\">No data available</text></svg>",
        (int)kW, (int)kH, (int)(kW / 2), (int)(kH / 2));
  }

  // ---- Compute data ranges --------------------------------------------
  double x_min = xs.front(), x_max = xs.back();
  double y_min = 0.0, y_max = 0.0;
  for (const auto &[key, ys] : partials) {
    for (double v : ys) {
      y_max = std::max(y_max, v);
      y_min = std::min(y_min, v);
    }
  }
  // Add 5% padding on y_max.
  y_max += (y_max - y_min) * 0.05;
  if (std::abs(y_max - y_min) < 1e-12) {
    y_min = 0.0;
    y_max = 1.0;
  }

  // ---- Build SVG -------------------------------------------------------
  std::ostringstream svg;

  // Header + background
  svg << std::format(
      "<svg xmlns=\"http://www.w3.org/2000/svg\" "
      "viewBox=\"0 0 {} {}\">\n"
      "  <rect width=\"100%\" height=\"100%\" fill=\"#1e1e2e\" rx=\"6\"/>\n",
      (int)kW, (int)kH);

  // Grid lines (y-axis, 5 divisions)
  constexpr int kYDivs = 5;
  svg << "  <!-- Y grid -->\n";
  for (int i = 0; i <= kYDivs; ++i) {
    double yv = y_min + i * (y_max - y_min) / kYDivs;
    double spy = detail::mapValue(yv, y_min, y_max, py1, py0);
    svg << std::format(
        "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
        "stroke=\"#45475a\" stroke-width=\"0.8\" stroke-dasharray=\"4,4\"/>\n",
        px0, spy, px1, spy);
    // Y tick label
    svg << std::format(
        "  <text x=\"{:.1f}\" y=\"{:.1f}\" fill=\"#a6adc8\" font-size=\"11\" "
        "text-anchor=\"end\" dominant-baseline=\"middle\">{}</text>\n",
        px0 - 5.0, spy, detail::fmtTick(yv));
  }

  // Grid lines (x-axis, 6 divisions)
  constexpr int kXDivs = 6;
  svg << "  <!-- X grid -->\n";
  for (int i = 0; i <= kXDivs; ++i) {
    double xv = x_min + i * (x_max - x_min) / kXDivs;
    double spx = detail::mapValue(xv, x_min, x_max, px0, px1);
    svg << std::format(
        "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
        "stroke=\"#45475a\" stroke-width=\"0.8\" stroke-dasharray=\"4,4\"/>\n",
        spx, py0, spx, py1);
    // X tick label
    svg << std::format(
        "  <text x=\"{:.1f}\" y=\"{:.1f}\" fill=\"#a6adc8\" font-size=\"11\" "
        "text-anchor=\"middle\">{}</text>\n",
        spx, py1 + 16.0, detail::fmtTick(xv));
  }

  // Axis borders
  svg << std::format(
      "  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" "
      "fill=\"none\" stroke=\"#6c7086\" stroke-width=\"1\"/>\n",
      px0, py0, px1 - px0, py1 - py0);

  // Data polylines + build legend entries
  std::vector<std::pair<std::string, std::string>> legend; // {label, color}
  std::size_t ci = 0;
  for (const auto &[key, ys] : partials) {
    const std::string col = detail::color(ci++);
    detail::appendPolyline(svg, xs, ys, x_min, x_max, y_min, y_max, px0, px1,
                           py0, py1, col, key);
    legend.push_back({key, col});
  }

  // Legend (top-right inside plot area)
  double leg_x = px1 - 8.0;
  double leg_y = py0 + 10.0;
  for (auto it = legend.rbegin(); it != legend.rend(); ++it) {
    svg << std::format(
        "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
        "stroke=\"{}\" stroke-width=\"2\"/>\n",
        leg_x - 22.0, leg_y, leg_x - 6.0, leg_y, it->second);
    svg << std::format(
        "  <text x=\"{:.1f}\" y=\"{:.1f}\" fill=\"#cdd6f4\" font-size=\"11\" "
        "text-anchor=\"end\" dominant-baseline=\"middle\">{}</text>\n",
        leg_x - 26.0, leg_y, it->first);
    leg_y += 16.0;
  }

  // Axis labels
  // X label
  svg << std::format(
      "  <text x=\"{:.1f}\" y=\"{:.1f}\" fill=\"#cdd6f4\" font-size=\"13\" "
      "text-anchor=\"middle\">{}</text>\n",
      (px0 + px1) / 2.0, kH - 8.0, x_label);

  // Y label (rotated)
  svg << std::format(
      "  <text transform=\"rotate(-90)\" x=\"{:.1f}\" y=\"{:.1f}\" "
      "fill=\"#cdd6f4\" font-size=\"13\" text-anchor=\"middle\">{}</text>\n",
      -((py0 + py1) / 2.0), 14.0, y_label);

  // Title
  svg << std::format(
      "  <text x=\"{:.1f}\" y=\"{:.1f}\" fill=\"#cdd6f4\" font-size=\"15\" "
      "font-weight=\"bold\" text-anchor=\"middle\">{}</text>\n",
      (px0 + px1) / 2.0, py0 - 12.0, title);

  svg << "</svg>\n";
  return svg.str();
}

} // namespace SvgPlotter
