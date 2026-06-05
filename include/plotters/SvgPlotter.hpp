/**
 * @file SvgPlotter.hpp
 * @brief Scientifically precise SVG generator for distribution function
 * histograms.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 *
 * @details
 * Produces a self-contained SVG string from a `Histogram` object with
 * "nice" tick intervals, scientific notation, and colorblind-safe palettes.
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "plotters/PathFont.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace correlation::plotters {

/**
 * @brief Theme and layout configuration for the plot.
 */
struct PlotConfig {
  /** @brief Visualization themes for the generated SVG. */
  enum class Theme {
    Light, ///< Standard light theme for publications.
    Dark   ///< Modern dark theme for UI integration (e.g., Catppuccin-esque).
  };

  Theme theme = Theme::Light; ///< Current theme selection.
  double width = 1200.0;      ///< SVG canvas width (px).
  double height = 900.0;      ///< SVG canvas height (px).
  bool show_grid = true;      ///< Whether to render background grid lines.
  bool show_markers = false;  ///< Whether to render data point markers (dots).

  // Publication settings
  double font_scale = 1.0;     ///< Multiplier for all font sizes
  double line_width = 3.0;     ///< Data line stroke width
  bool show_legend = true;     ///< Toggle legend visibility

  /** @brief Color palette selections */
  enum class Palette {
    OkabeIto,  ///< Colorblind-safe palette
    Grayscale, ///< B&W printing friendly
    Viridis    ///< Perceptually uniform
  };
  Palette palette = Palette::OkabeIto;

  /** @brief Standard publication sizes */
  enum class PresetSize {
    Default,
    SingleColumn,
    DoubleColumn,
    Presentation
  };
  PresetSize preset_size = PresetSize::Default;

  double effective_width() const {
    switch (preset_size) {
      case PresetSize::SingleColumn: return 1050.0;
      case PresetSize::DoubleColumn: return 2100.0;
      case PresetSize::Presentation: return 3000.0;
      default: return width;
    }
  }

  double effective_height() const {
    switch (preset_size) {
      case PresetSize::SingleColumn: return 788.0;
      case PresetSize::DoubleColumn: return 1575.0;
      case PresetSize::Presentation: return 2250.0;
      default: return height;
    }
  }

  /** @return Hex color string for the plot background. */
  std::string bg_color() const { return (theme == Theme::Light) ? "#FFFFFF" : "#1e1e2e"; }
  /** @return Hex color string for axes and ticks. */
  std::string axis_color() const { return (theme == Theme::Light) ? "#000000" : "#cdd6f4"; }
  /** @return Hex color string for grid lines. */
  std::string grid_color() const { return (theme == Theme::Light) ? "#e0e0e0" : "#45475a"; }
  /** @return Hex color string for labels and titles. */
  std::string text_color() const { return (theme == Theme::Light) ? "#333333" : "#a6adc8"; }
};

/**
 * @brief Information about the mouse hover position.
 */
struct HoverInfo {
  bool active = false;
  double mouse_x = -1.0;
  double mouse_y = -1.0;
  double widget_width = 0.0;
  double widget_height = 0.0;
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
    "#000000", // Black
};

/// Grayscale palette for B&W printing.
static const std::vector<std::string> kGrayscale = {
    "#000000", "#404040", "#808080", "#B0B0B0", "#D0D0D0"
};

/// Viridis perceptually uniform palette.
static const std::vector<std::string> kViridis = {
    "#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"
};

/**
 * @brief Retrieves a color from the selected palette.
 * @param i Index.
 * @param pal Selected palette.
 * @return Hex color string.
 */
inline std::string color(std::size_t i, PlotConfig::Palette pal) {
  switch (pal) {
    case PlotConfig::Palette::Grayscale: return kGrayscale[i % kGrayscale.size()];
    case PlotConfig::Palette::Viridis:   return kViridis[i % kViridis.size()];
    default:                              return kColors[i % kColors.size()];
  }
}

/**
 * @brief Maps a value from data space to SVG coordinate space.
 * @param v The coordinate in data space.
 * @param data_min Minimum data value.
 * @param data_max Maximum data value.
 * @param svg_min Target SVG start coordinate.
 * @param svg_max Target SVG end coordinate.
 * @return The mapped SVG coordinate.
 */
inline double mapValue(double v, double data_min, double data_max, double svg_min, double svg_max) {
  if (std::abs(data_max - data_min) < 1e-15)
    return (svg_min + svg_max) / 2.0;
  return svg_min + (v - data_min) / (data_max - data_min) * (svg_max - svg_min);
}

/**
 * @brief Logic for generating "nice" human-readable tick intervals.
 */
struct NiceScale {
  double min;                ///< Starting tick value.
  double max;                ///< Ending tick value.
  double spacing;            ///< Calculated distance between ticks.
  std::vector<double> ticks; ///< Generated tick locations.

  /**
   * @brief Calculates nice intervals for a range.
   * @param actual_min Measured data minimum.
   * @param actual_max Measured data maximum.
   * @param max_ticks Target number of ticks.
   */
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
  /**
   * @brief Rounds a range value to a "nice" human-readable number.
   * @param range The value range for the axis.
   * @param round Whether to perform aggressive rounding.
   * @return The rounded "nice" value.
   */
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
 *
 * Automatically converts "x 10^y" to "x 10^y" format using Unicode superscript
 * characters to avoid complex SVG text positioning issues.
 *
 * @param v The number to format.
 * @return A UTF-8 string ready for SVG path rendering.
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
      if (c == '-')
        res += "⁻";
      else if (c == '0')
        res += "⁰";
      else if (c == '1')
        res += "¹";
      else if (c == '2')
        res += "²";
      else if (c == '3')
        res += "³";
      else if (c == '4')
        res += "⁴";
      else if (c == '5')
        res += "⁵";
      else if (c == '6')
        res += "⁶";
      else if (c == '7')
        res += "⁷";
      else if (c == '8')
        res += "⁸";
      else if (c == '9')
        res += "⁹";
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
 * @brief Renders text as a filled SVG path using the Roboto outline font.
 * Uses evenodd fill-rule to render the font's inner holes properly.
 */
inline std::string renderTextAsPath(const std::string &text, double x, double y, double size, const std::string &anchor, const std::string &color) {
  std::string path_d = Roboto::instance().render(text, x, y, size, anchor);
  return std::format("  <path d=\"{}\" fill=\"{}\" fill-rule=\"evenodd\" stroke=\"none\"/>\n", path_d, color);
}

/**
 * @brief Renders a `Histogram` as a self-contained SVG string.
 *
 * @param hist     The Histogram to render.
 * @param config   Optional plot configuration (theme, size, etc.).
 * @param hover    Optional hover interaction info.
 * @returns        A complete SVG document as `std::string`.
 */
inline std::string renderHistogramAsSvg(const correlation::analysis::Histogram &hist, const PlotConfig &config = {}, const HoverInfo &hover = {}) {
  std::string title = hist.title.empty() ? "Histogram" : hist.title;
  std::string x_label = hist.x_label.empty() ? "x" : hist.x_label;
  std::string y_label = hist.y_label.empty() ? "y" : hist.y_label;

  // Add units if available
  if (!hist.x_unit.empty())
    x_label += std::format(" ({})", hist.x_unit);
  if (!hist.y_unit.empty()) {
    std::string y_unit = hist.y_unit;
    // Standardize Å^-1 to Å⁻¹ for path rendering
    if (y_unit == "Å^-1")
      y_unit = "Å⁻¹";
    y_label += std::format(" ({})", y_unit);
  }

  // ---- Layout configuration --------------------------------------------
  const double kW = config.effective_width();
  const double kH = config.effective_height();
  // Margins (enough room for labels)
  const double kLeft = 100.0;
  const double kRight = 40.0;
  const double kTop = 50.0;
  const double kBot = 90.0;

  const double px0 = kLeft;
  const double px1 = kW - kRight;
  const double py0 = kTop;
  const double py1 = kH - kBot;

  // ---- Data source -----------------------------------------------------
  const auto &partials = hist.smoothed_partials.empty() ? hist.partials : hist.smoothed_partials;
  const auto &xs = hist.bins;

  if (xs.empty() || partials.empty()) {
    std::string no_data_path = Roboto::instance().render("No data available", kW / 2.0, kH / 2.0 + 8.0, 24 * config.font_scale, "middle");
    return std::format("<svg width='{0:.0f}' height='{1:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {0:.0f} {1:.0f}\">"
                       "<rect width=\"100%\" height=\"100%\" fill=\"{2}\"/>"
                       "<path d=\"{3}\" fill=\"{4}\" fill-rule=\"evenodd\" stroke=\"none\"/></svg>",
                       kW, kH, config.bg_color(), no_data_path, config.text_color());
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

  detail::NiceScale xScale(raw_x_min, raw_x_max, 11);
  detail::NiceScale yScale(raw_y_min, raw_y_max, 8);

  // ---- Build SVG -------------------------------------------------------
  std::ostringstream svg;
  svg << std::format("<svg width='{:.0f}' height='{:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\">\n", kW,
                     kH, kW, kH);
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
    svg << renderTextAsPath(detail::fmtScientific(yv), px0 - 15.0, spy + 7.0, 20.0 * config.font_scale, "end", config.text_color());
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
    svg << renderTextAsPath(detail::fmtScientific(xv), spx, py1 + 25.0, 20.0 * config.font_scale, "middle", config.text_color());
  }

  // Draw axis border
  svg << std::format("  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" "
                     "fill=\"none\" stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                     px0, py0, px1 - px0, py1 - py0, config.axis_color());

  // Emphasis lines for y=0 or y=1
  std::vector<double> emphasis_values;
  if (title.find("g(r)") != std::string::npos) {
    emphasis_values.push_back(1.0);
  } else if (title.find("G(r)") != std::string::npos) {
    emphasis_values.push_back(0.0);
  } else if (title.find("S(Q)") != std::string::npos || title.find("S(q)") != std::string::npos) {
    emphasis_values.push_back(0.0);
    emphasis_values.push_back(1.0);
  }

  for (double focus_y : emphasis_values) {
    if (focus_y >= yScale.min && focus_y <= yScale.max) {
      double spy = detail::mapValue(focus_y, yScale.min, yScale.max, py1, py0);
      std::string extra = (std::abs(focus_y - 1.0) < 1e-6) ? " stroke-dasharray=\"5,5\"" : "";
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.8\"{}/>\n",
                         px0, spy, px1, spy, config.axis_color(), extra);
    }
  }

  // Polylines
  std::vector<std::pair<std::string, std::string>> legend;
  std::size_t ci = 0;
  for (const auto &[key, ys] : partials) {
    const std::string col = detail::color(ci++, config.palette);
    svg << std::format("  <polyline fill=\"none\" stroke=\"{}\" stroke-width=\"{:.1f}\" stroke-linejoin=\"round\" points=\"",
                       col, config.line_width);
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
  if (config.show_legend) {
    double lx = px1 - 15.0;
    double ly = py0 + 25.0;
    for (auto it = legend.rbegin(); it != legend.rend(); ++it) {
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"4.0\"/>\n",
                         lx - 40.0, ly, lx - 10.0, ly, it->second);
      svg << renderTextAsPath(it->first, lx - 45.0, ly + 6.0, 18.0 * config.font_scale, "end", config.text_color());
      ly += 28.0;
    }
  }

  // Titles/Labels
  svg << renderTextAsPath(x_label, (px0 + px1) / 2.0, py1 + 75.0, 28.0 * config.font_scale, "middle", config.text_color());

  // Y label rotated
  svg << std::format("  <g transform=\"translate({:.1f}, {:.1f}) rotate(-90)\">\n", 40.0, (py0 + py1) / 2.0);
  svg << renderTextAsPath(y_label, 0.0, 0.0, 28.0 * config.font_scale, "middle", config.text_color());
  svg << "  </g>\n";

  // Hover tracking interaction
  if (hover.active && hover.widget_width > 0.0 && hover.widget_height > 0.0) {
    double svg_aspect = kW / kH;
    double widget_aspect = hover.widget_width / hover.widget_height;
    double scale = 1.0;
    double dx = 0.0;
    double dy = 0.0;

    if (widget_aspect > svg_aspect) {
      scale = hover.widget_height / kH;
      dx = (hover.widget_width - kW * scale) / 2.0;
    } else {
      scale = hover.widget_width / kW;
      dy = (hover.widget_height - kH * scale) / 2.0;
    }

    double sx = (hover.mouse_x - dx) / scale;
    double sy = (hover.mouse_y - dy) / scale;

    if (sx >= px0 && sx <= px1) {
      // Map to data space
      double x_val = xScale.min + (sx - px0) / (px1 - px0) * (xScale.max - xScale.min);
      
      // Find nearest index in bins
      auto it = std::lower_bound(xs.begin(), xs.end(), x_val);
      std::size_t idx = 0;
      if (it == xs.end()) {
        idx = xs.size() - 1;
      } else if (it == xs.begin()) {
        idx = 0;
      } else {
        double d1 = *it - x_val;
        double d2 = x_val - *(it - 1);
        if (d2 < d1) {
          idx = std::distance(xs.begin(), it - 1);
        } else {
          idx = std::distance(xs.begin(), it);
        }
      }

      double target_x = xs[idx];
      double sx_data = detail::mapValue(target_x, xScale.min, xScale.max, px0, px1);

      // Draw vertical guide line
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.5\" stroke-dasharray=\"4,4\"/>\n",
                         sx_data, py0, sx_data, py1, config.axis_color());

      // Collect values and draw markers on curves
      std::vector<std::tuple<std::string, double, std::string>> hover_values; // name, value, color
      std::size_t ci = 0;
      for (const auto &[key, ys] : partials) {
        if (idx < ys.size()) {
          const std::string col = detail::color(ci++, config.palette);
          double y_val = ys[idx];
          double sy_data = detail::mapValue(y_val, yScale.min, yScale.max, py1, py0);

          // Bullet marker
          svg << std::format("  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"6\" fill=\"{}\" stroke=\"{}\" stroke-width=\"2\"/>\n",
                             sx_data, sy_data, col, config.bg_color());
          
          hover_values.push_back({key, y_val, col});
        }
      }

      // Draw Tooltip Box
      if (!hover_values.empty()) {
        double tooltip_w = 200.0;
        double tooltip_h = 35.0 + 22.0 * hover_values.size();
        
        // Tooltip position (flip sides depending on cursor location)
        double tx = (sx_data < kW / 2.0) ? sx_data + 15.0 : sx_data - tooltip_w - 15.0;
        double ty = std::clamp(sy - 30.0, py0 + 10.0, py1 - tooltip_h - 10.0);

        std::string card_bg = (config.theme == PlotConfig::Theme::Light) ? "#FFFFFF" : "#181825";
        std::string card_border = (config.theme == PlotConfig::Theme::Light) ? "#dddddd" : "#45475a";
        std::string text_col = (config.theme == PlotConfig::Theme::Light) ? "#333333" : "#cdd6f4";

        svg << std::format("  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" rx=\"6\" "
                           "fill=\"{}\" fill-opacity=\"0.92\" stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                           tx, ty, tooltip_w, tooltip_h, card_bg, card_border);

        // Header: x value (with unit or pure label)
        std::string x_unit_str = hist.x_unit;
        std::string header_txt = std::format("{} = {:.4f}{}", x_label, target_x, x_unit_str.empty() ? "" : " " + x_unit_str);
        // Clean up title (remove parenthesis unit if present, e.g. "r (Å)" to "r")
        auto paren = x_label.find(" (");
        if (paren != std::string::npos) {
          header_txt = std::format("{} = {:.4f}", x_label.substr(0, paren), target_x);
        }
        
        svg << renderTextAsPath(header_txt, tx + 12.0, ty + 20.0, 14.0 * config.font_scale, "start", text_col);

        double cur_y = ty + 42.0;
        for (const auto &[name, val, col] : hover_values) {
          // Color swatch dot
          svg << std::format("  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"5\" fill=\"{}\"/>\n",
                             tx + 18.0, cur_y - 4.0, col);
          // Label and value
          std::string line_txt = std::format("{}: {:.4f}", name, val);
          svg << renderTextAsPath(line_txt, tx + 30.0, cur_y, 13.0 * config.font_scale, "start", text_col);
          cur_y += 22.0;
        }
      }
    }
  }

  svg << "</svg>\n";
  return svg.str();
}

// ============================================================================
// Comparison overlay renderer
// ============================================================================

/**
 * @brief A labeled histogram for comparison rendering.
 */
struct LabeledHistogram {
  std::string label;                            ///< Run / dataset label.
  const correlation::analysis::Histogram *hist; ///< Pointer to histogram data.
};

/**
 * @brief Renders multiple histograms overlaid on a single plot.
 *
 * Each entry in @p datasets provides a display label and a pointer to the
 * Histogram data.  The first dataset defines the axes and labels; subsequent
 * datasets are drawn on the same grid with distinct colors.
 *
 * @param datasets  Vector of labeled histogram pointers.
 * @param partial_key  Which partial to plot (e.g. "Total", "Si-O").
 * @param config    Optional plot configuration (theme, size, etc.).
 * @param hover     Optional hover interaction info.
 * @returns         A complete SVG document as `std::string`.
 */
inline std::string renderComparisonSvg(const std::vector<LabeledHistogram> &datasets,
                                       const std::string &partial_key = "Total", const PlotConfig &config = {}, const HoverInfo &hover = {}) {
  if (datasets.empty())
    return "";

  // Use the first histogram for titles and axis labels.
  const auto &ref = *datasets.front().hist;
  std::string title = ref.title.empty() ? "Comparison" : ref.title;
  std::string x_label = ref.x_label.empty() ? "x" : ref.x_label;
  std::string y_label = ref.y_label.empty() ? "y" : ref.y_label;

  if (!ref.x_unit.empty())
    x_label += std::format(" ({})", ref.x_unit);
  if (!ref.y_unit.empty()) {
    std::string y_unit = ref.y_unit;
    if (y_unit == "Å^-1")
      y_unit = "Å⁻¹";
    y_label += std::format(" ({})", y_unit);
  }

  const double kW = config.effective_width();
  const double kH = config.effective_height();
  const double kLeft = 100.0;
  const double kRight = 40.0;
  const double kTop = 50.0;
  const double kBot = 90.0;

  const double px0 = kLeft;
  const double px1 = kW - kRight;
  const double py0 = kTop;
  const double py1 = kH - kBot;

  // Compute global ranges across all datasets.
  double raw_x_min = 1e30, raw_x_max = -1e30;
  double raw_y_min = 0.0, raw_y_max = -1e30;

  for (const auto &ds : datasets) {
    const auto &partials = ds.hist->smoothed_partials.empty() ? ds.hist->partials : ds.hist->smoothed_partials;
    auto it = partials.find(partial_key);
    if (it == partials.end())
      continue;
    const auto &xs = ds.hist->bins;
    const auto &ys = it->second;
    if (!xs.empty()) {
      raw_x_min = std::min(raw_x_min, xs.front());
      raw_x_max = std::max(raw_x_max, xs.back());
    }
    for (double v : ys) {
      raw_y_max = std::max(raw_y_max, v);
      raw_y_min = std::min(raw_y_min, v);
    }
  }

  if (raw_x_max <= raw_x_min || raw_y_max <= raw_y_min) {
    std::string no_data_path = Roboto::instance().render("No comparison data", kW / 2.0, kH / 2.0 + 8.0, 24 * config.font_scale, "middle");
    return std::format("<svg width='{}' height='{}' xmlns=\"http://www.w3.org/2000/svg\" "
                       "viewBox=\"0 0 {} {}\">"
                       "<rect width=\"100%\" height=\"100%\" fill=\"{}\"/>"
                       "<path d=\"{}\" fill=\"{}\" fill-rule=\"evenodd\" stroke=\"none\"/></svg>",
                       kW, kH, kW, kH, config.bg_color(), no_data_path, config.text_color());
  }

  double y_padding = (raw_y_max - raw_y_min) * 0.05;
  raw_y_max += y_padding;

  detail::NiceScale xScale(raw_x_min, raw_x_max, 11);
  detail::NiceScale yScale(raw_y_min, raw_y_max, 8);

  std::ostringstream svg;
  svg << std::format("<svg width='{:.0f}' height='{:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\">\n",
                     kW, kH, kW, kH);
  svg << std::format("  <rect width=\"100%\" height=\"100%\" fill=\"{}\" rx=\"6\"/>\n", config.bg_color());

  // Grid, ticks, and axes — same as single-plot renderer.
  for (double yv : yScale.ticks) {
    double spy = detail::mapValue(yv, yScale.min, yScale.max, py1, py0);
    if (config.show_grid) {
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"0.8\" stroke-dasharray=\"3,3\"/>\n",
                         px0, spy, px1, spy, config.grid_color());
    }
    svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                       px0 - 8.0, spy, px0, spy, config.axis_color());
    svg << renderTextAsPath(detail::fmtScientific(yv), px0 - 15.0, spy + 7.0, 20.0 * config.font_scale, "end", config.text_color());
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
    svg << renderTextAsPath(detail::fmtScientific(xv), spx, py1 + 25.0, 20.0 * config.font_scale, "middle", config.text_color());
  }

  // Axis border
  svg << std::format("  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" "
                     "fill=\"none\" stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                     px0, py0, px1 - px0, py1 - py0, config.axis_color());

  // Polylines — one per dataset
  std::vector<std::pair<std::string, std::string>> legend;
  std::size_t ci = 0;
  for (const auto &ds : datasets) {
    const auto &partials = ds.hist->smoothed_partials.empty() ? ds.hist->partials : ds.hist->smoothed_partials;
    auto it = partials.find(partial_key);
    if (it == partials.end())
      continue;

    const auto &xs = ds.hist->bins;
    const auto &ys = it->second;
    const std::string col = detail::color(ci++, config.palette);

    svg << std::format("  <polyline fill=\"none\" stroke=\"{}\" stroke-width=\"{:.1f}\" stroke-linejoin=\"round\" points=\"",
                       col, config.line_width);
    std::size_t n = std::min(xs.size(), ys.size());
    for (std::size_t i = 0; i < n; ++i) {
      double sx = detail::mapValue(xs[i], xScale.min, xScale.max, px0, px1);
      double sy = detail::mapValue(ys[i], yScale.min, yScale.max, py1, py0);
      svg << std::format("{:.2f},{:.2f} ", sx, sy);
    }
    svg << "\" />\n";
    legend.push_back({ds.label, col});
  }

  // Legend (top right)
  if (config.show_legend) {
    double lx = px1 - 15.0;
    double ly = py0 + 25.0;
    for (auto it = legend.rbegin(); it != legend.rend(); ++it) {
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"4.0\"/>\n",
                         lx - 40.0, ly, lx - 10.0, ly, it->second);
      svg << renderTextAsPath(it->first, lx - 45.0, ly + 6.0, 18.0 * config.font_scale, "end", config.text_color());
      ly += 28.0;
    }
  }

  // X-axis label
  svg << renderTextAsPath(x_label, (px0 + px1) / 2.0, py1 + 75.0, 28.0 * config.font_scale, "middle", config.text_color());

  // Y-axis label (rotated)
  svg << std::format("  <g transform=\"translate({:.1f}, {:.1f}) rotate(-90)\">\n", 40.0, (py0 + py1) / 2.0);
  svg << renderTextAsPath(y_label, 0.0, 0.0, 28.0 * config.font_scale, "middle", config.text_color());
  svg << "  </g>\n";

  // Hover tracking interaction for comparison plot
  if (hover.active && hover.widget_width > 0.0 && hover.widget_height > 0.0) {
    double svg_aspect = kW / kH;
    double widget_aspect = hover.widget_width / hover.widget_height;
    double scale = 1.0;
    double dx = 0.0;
    double dy = 0.0;

    if (widget_aspect > svg_aspect) {
      scale = hover.widget_height / kH;
      dx = (hover.widget_width - kW * scale) / 2.0;
    } else {
      scale = hover.widget_width / kW;
      dy = (hover.widget_height - kH * scale) / 2.0;
    }

    double sx = (hover.mouse_x - dx) / scale;
    double sy = (hover.mouse_y - dy) / scale;

    if (sx >= px0 && sx <= px1) {
      // Map to data space
      double x_val = xScale.min + (sx - px0) / (px1 - px0) * (xScale.max - xScale.min);
      
      // We assume first dataset's bins for the X coordinates
      const auto &xs = datasets.front().hist->bins;
      if (!xs.empty()) {
        auto it = std::lower_bound(xs.begin(), xs.end(), x_val);
        std::size_t idx = 0;
        if (it == xs.end()) {
          idx = xs.size() - 1;
        } else if (it == xs.begin()) {
          idx = 0;
        } else {
          double d1 = *it - x_val;
          double d2 = x_val - *(it - 1);
          if (d2 < d1) {
            idx = std::distance(xs.begin(), it - 1);
          } else {
            idx = std::distance(xs.begin(), it);
          }
        }

        double target_x = xs[idx];
        double sx_data = detail::mapValue(target_x, xScale.min, xScale.max, px0, px1);

        // Draw vertical guide line
        svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"1.5\" stroke-dasharray=\"4,4\"/>\n",
                           sx_data, py0, sx_data, py1, config.axis_color());

        // Collect values from each dataset and draw markers
        std::vector<std::tuple<std::string, double, std::string>> hover_values; // label, value, color
        std::size_t ci = 0;
        for (const auto &ds : datasets) {
          const auto &partials = ds.hist->smoothed_partials.empty() ? ds.hist->partials : ds.hist->smoothed_partials;
          auto pit = partials.find(partial_key);
          if (pit != partials.end() && idx < pit->second.size()) {
            const std::string col = detail::color(ci++, config.palette);
            double y_val = pit->second[idx];
            double sy_data = detail::mapValue(y_val, yScale.min, yScale.max, py1, py0);

            // Bullet marker
            svg << std::format("  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"6\" fill=\"{}\" stroke=\"{}\" stroke-width=\"2\"/>\n",
                               sx_data, sy_data, col, config.bg_color());
            
            hover_values.push_back({ds.label, y_val, col});
          } else {
            // Keep the color index aligned
            ci++;
          }
        }

        // Draw Tooltip Box
        if (!hover_values.empty()) {
          double tooltip_w = 200.0;
          double tooltip_h = 35.0 + 22.0 * hover_values.size();
          
          double tx = (sx_data < kW / 2.0) ? sx_data + 15.0 : sx_data - tooltip_w - 15.0;
          double ty = std::clamp(sy - 30.0, py0 + 10.0, py1 - tooltip_h - 10.0);

          std::string card_bg = (config.theme == PlotConfig::Theme::Light) ? "#FFFFFF" : "#181825";
          std::string card_border = (config.theme == PlotConfig::Theme::Light) ? "#dddddd" : "#45475a";
          std::string text_col = (config.theme == PlotConfig::Theme::Light) ? "#333333" : "#cdd6f4";

          svg << std::format("  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" rx=\"6\" "
                             "fill=\"{}\" fill-opacity=\"0.92\" stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                             tx, ty, tooltip_w, tooltip_h, card_bg, card_border);

          // Header
          std::string x_unit_str = ref.x_unit;
          std::string header_txt = std::format("{} = {:.4f}{}", x_label, target_x, x_unit_str.empty() ? "" : " " + x_unit_str);
          auto paren = x_label.find(" (");
          if (paren != std::string::npos) {
            header_txt = std::format("{} = {:.4f}", x_label.substr(0, paren), target_x);
          }

          svg << std::format("  <text x=\"{:.1f}\" y=\"{:.1f}\" font-family=\"'Inter', 'Roboto', sans-serif\" "
                             "font-size=\"14\" font-weight=\"bold\" fill=\"{}\">{}</text>\n",
                             tx + 12.0, ty + 24.0, text_col, header_txt);

          double cur_y = ty + 46.0;
          for (const auto &[name, val, col] : hover_values) {
            svg << std::format("  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"5\" fill=\"{}\"/>\n",
                               tx + 18.0, cur_y - 5.0, col);
            svg << std::format("  <text x=\"{:.1f}\" y=\"{:.1f}\" font-family=\"'Inter', 'Roboto', sans-serif\" "
                               "font-size=\"13\" fill=\"{}\">{}: {:.4f}</text>\n",
                               tx + 30.0, cur_y, text_col, name, val);
            cur_y += 22.0;
          }
        }
      }
    }
  }

  svg << "</svg>\n";
  return svg.str();
}

} // namespace correlation::plotters
