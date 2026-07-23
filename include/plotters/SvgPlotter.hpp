/**
 * @file SvgPlotter.hpp
 * @brief Zero-dependency SVG plotter for liquid and amorphous data.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "math/Precision.hpp"
#include "plotters/PathFont.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <map>
#include <ranges>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace correlation::plotters {

/**
 * @brief Theme and layout configuration for the plot.
 */
struct PlotConfig {
  /** @brief Visualization themes for the generated SVG. */
  enum class Theme : std::uint8_t {
    Light, ///< Standard light theme for publications.
    Dark   ///< Modern dark theme for UI integration (e.g., Catppuccin-esque).
  };

  Theme theme = Theme::Light;                 ///< Current theme selection.
  real_t width = static_cast<real_t>(1200.0); ///< SVG canvas width (px).
  real_t height = static_cast<real_t>(900.0); ///< SVG canvas height (px).
  bool show_grid = true;                      ///< Whether to render background grid lines.
  bool show_markers = false;                  ///< Whether to render data point markers (dots).
  bool fill_area = false;                     ///< Whether to render matching gradient fills under the curves.

  // Publication settings
  real_t font_scale = static_cast<real_t>(1.0); ///< Multiplier for all font sizes
  real_t line_width = static_cast<real_t>(3.0); ///< Data line stroke width
  bool show_legend = true;                      ///< Toggle legend visibility
  bool use_native_text = false;                 ///< Use standard SVG &lt;text&gt; elements instead of Hershey paths.

  /** @brief Color palette selections */
  enum class Palette : std::uint8_t {
    OkabeIto,  ///< Colorblind-safe palette
    Grayscale, ///< B&W printing friendly
    Viridis    ///< Perceptually uniform
  };
  Palette palette = Palette::OkabeIto;

  /** @brief Standard publication sizes */
  enum class PresetSize : std::uint8_t { Default, SingleColumn, DoubleColumn, Presentation };
  PresetSize preset_size = PresetSize::Default;

  real_t effective_width() const {
    switch (preset_size) {
    case PresetSize::SingleColumn:
      return static_cast<real_t>(1050.0);
    case PresetSize::DoubleColumn:
      return static_cast<real_t>(2100.0);
    case PresetSize::Presentation:
      return static_cast<real_t>(3000.0);
    default:
      return width;
    }
  }

  real_t effective_height() const {
    switch (preset_size) {
    case PresetSize::SingleColumn:
      return static_cast<real_t>(788.0);
    case PresetSize::DoubleColumn:
      return static_cast<real_t>(1575.0);
    case PresetSize::Presentation:
      return static_cast<real_t>(2250.0);
    default:
      return height;
    }
  }

  /** @return Hex color string for the plot background. */
  std::string bg_color() const { return (theme == Theme::Light) ? "#FFFFFF" : "#1e1e2e"; }
  /** @return Hex color string for axes and ticks. */
  std::string axis_color() const { return (theme == Theme::Light) ? "#000000" : "#cdd6f4"; }
  /** @return Hex color string for grid lines. */
  std::string grid_color() const { return (theme == Theme::Light) ? "#808080" : "#45475a"; }
  /** @return Hex color string for labels and titles. */
  std::string text_color() const { return (theme == Theme::Light) ? "#333333" : "#a6adc8"; }
};

/**
 * @brief Information about the mouse hover position.
 */
struct HoverInfo {
  bool active = false;
  real_t mouse_x = static_cast<real_t>(-1.0);
  real_t mouse_y = static_cast<real_t>(-1.0);
  real_t widget_width = static_cast<real_t>(0.0);
  real_t widget_height = static_cast<real_t>(0.0);
};

enum class TextAnchor : std::uint8_t { Start, Middle, End };

inline std::string renderTextAsPath(const std::string &text, real_t x_pos, real_t y_pos, real_t size, TextAnchor anchor,
                                    const std::string &color, bool use_native_text);

// ============================================================================
// Internal helpers
// ============================================================================
/**
 * @brief A labeled histogram for comparison rendering.
 */
struct LabeledHistogram {
  std::string label;                            ///< Run / dataset label.
  const correlation::analysis::Histogram *hist; ///< Pointer to histogram data.
};

namespace detail {

/// Okabe-Ito colorblind-safe palette.
constexpr std::array<std::string_view, 8> kColors = {
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
constexpr std::array<std::string_view, 5> kGrayscale = {"#000000", "#404040", "#808080", "#B0B0B0", "#D0D0D0"};

/// Viridis perceptually uniform palette.
constexpr std::array<std::string_view, 5> kViridis = {"#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"};

/**
 * @brief Retrieves a color from the selected palette.
 * @param index Index.
 * @param pal Selected palette.
 * @return Hex color string.
 */
inline std::string color(std::size_t index, PlotConfig::Palette pal) {
  switch (pal) {
  case PlotConfig::Palette::Grayscale:
    return std::string(kGrayscale.at(index % kGrayscale.size()));
  case PlotConfig::Palette::Viridis:
    return std::string(kViridis.at(index % kViridis.size()));
  default:
    return std::string(kColors.at(index % kColors.size()));
  }
}

/**
 * @brief Maps a value from data space to SVG coordinate space.
 * @param value The coordinate in data space.
 * @param data_min Minimum data value.
 * @param data_max Maximum data value.
 * @param svg_min Target SVG start coordinate.
 * @param svg_max Target SVG end coordinate.
 * @return The mapped SVG coordinate.
 */
inline real_t mapValue(real_t value, real_t data_min, real_t data_max, real_t svg_min, real_t svg_max) {
  if (std::abs(data_max - data_min) < static_cast<real_t>(1e-15)) {
    return (svg_min + svg_max) / static_cast<real_t>(2.0);
  }
  return svg_min + (value - data_min) / (data_max - data_min) * (svg_max - svg_min);
}

struct DataRange {
  real_t min = static_cast<real_t>(0.0);
  real_t max = static_cast<real_t>(0.0);
};

/**
 * @brief Logic for generating "nice" human-readable tick intervals.
 */
struct NiceScale {
  real_t min = static_cast<real_t>(0.0);     ///< Starting tick value.
  real_t max = static_cast<real_t>(0.0);     ///< Ending tick value.
  real_t spacing = static_cast<real_t>(0.0); ///< Calculated distance between ticks.
  std::vector<real_t> ticks;                 ///< Generated tick locations.

  NiceScale() = default;

  /**
   * @brief Calculates nice intervals for a range.
   * @param range The measured data range (min and max).
   * @param max_ticks Target number of ticks.
   */
  explicit NiceScale(const DataRange &range, int max_ticks = 6) {
    real_t actual_min = range.min;
    real_t actual_max = range.max;
    if (std::abs(actual_max - actual_min) < static_cast<real_t>(1e-12)) {
      min = actual_min - static_cast<real_t>(0.5);
      max = actual_min + static_cast<real_t>(0.5);
      spacing = static_cast<real_t>(0.1);
    } else {
      real_t range_val = niceNum(actual_max - actual_min, false);
      spacing = niceNum(range_val / static_cast<real_t>(max_ticks - 1), true);
      min = std::floor(actual_min / spacing) * spacing;
      max = std::ceil(actual_max / spacing) * spacing;
    }

    real_t range_span = max - min;
    int num_ticks = static_cast<int>(std::round(range_span / spacing)) + 1;
    for (int idx = 0; idx < num_ticks; ++idx) {
      real_t value = min + static_cast<real_t>(idx) * spacing;
      ticks.push_back(value);
    }
  }

private:
  /**
   * @brief Rounds a range value to a "nice" human-readable number.
   * @param range The value range for the axis.
   * @param round Whether to perform aggressive rounding.
   * @return The rounded "nice" value.
   */
  static real_t niceNum(real_t range, bool round) {
    real_t exponent = std::floor(std::log10(range));
    real_t fraction = range / std::pow(static_cast<real_t>(10.0), exponent);
    real_t nice_fraction = static_cast<real_t>(0.0);

    if (round) {
      if (fraction < static_cast<real_t>(1.5)) {
        nice_fraction = static_cast<real_t>(1.0);
      } else if (fraction < static_cast<real_t>(3.0)) {
        nice_fraction = static_cast<real_t>(2.0);
      } else if (fraction < static_cast<real_t>(7.0)) {
        nice_fraction = static_cast<real_t>(5.0);
      } else {
        nice_fraction = static_cast<real_t>(10.0);
      }
    } else {
      if (fraction <= static_cast<real_t>(1.0)) {
        nice_fraction = static_cast<real_t>(1.0);
      } else if (fraction <= static_cast<real_t>(2.0)) {
        nice_fraction = static_cast<real_t>(2.0);
      } else if (fraction <= static_cast<real_t>(5.0)) {
        nice_fraction = static_cast<real_t>(5.0);
      } else {
        nice_fraction = static_cast<real_t>(10.0);
      }
    }
    return nice_fraction * std::pow(static_cast<real_t>(10.0), exponent);
  }
};

/**
 * @brief Formats a number for SVG display, using scientific notation if needed.
 *
 * Automatically converts "x 10^y" to "x 10^y" format using Unicode superscript
 * characters to avoid complex SVG text positioning issues.
 *
 * @param value The number to format.
 * @return A UTF-8 string ready for SVG path rendering.
 */
inline std::string fmtScientific(real_t value) {
  real_t abs_value = std::abs(value);
  if (abs_value < static_cast<real_t>(1e-12)) {
    return "0";
  }

  // Use scientific notation for very small/large values
  if (abs_value < static_cast<real_t>(0.001) || abs_value >= static_cast<real_t>(10000.0)) {
    int exponent = static_cast<int>(std::floor(std::log10(abs_value)));
    real_t fraction = value / std::pow(static_cast<real_t>(10.0), static_cast<real_t>(exponent));
    std::string res = std::format("{:.1f}×10", fraction);
    std::string exp_s = std::to_string(exponent);
    for (char chr : exp_s) {
      if (chr == '-') {
        res += "⁻";
      } else if (chr == '0') {
        res += "⁰";
      } else if (chr == '1') {
        res += "¹";
      } else if (chr == '2') {
        res += "²";
      } else if (chr == '3') {
        res += "³";
      } else if (chr == '4') {
        res += "⁴";
      } else if (chr == '5') {
        res += "⁵";
      } else if (chr == '6') {
        res += "⁶";
      } else if (chr == '7') {
        res += "⁷";
      } else if (chr == '8') {
        res += "⁸";
      } else if (chr == '9') {
        res += "⁹";
      }
    }
    return res;
  }

  // Otherwise use simple formatting
  std::string str = std::format("{:.2f}", value);
  auto dot = str.find('.');
  if (dot != std::string::npos) {
    std::size_t last = str.find_last_not_of('0');
    if (last != std::string::npos && last > dot) {
      str = str.substr(0, last + 1);
    } else if (last == dot) {
      str = str.substr(0, dot);
    }
  }
  return str;
}

struct SvgHistogramRenderer {
  const correlation::analysis::Histogram *hist = nullptr;
  const PlotConfig *config = nullptr;
  const HoverInfo *hover = nullptr;
  const std::map<std::string, real_t> *weights = nullptr;

  std::string title;
  std::string x_label;
  std::string y_label;

  real_t kW = static_cast<real_t>(0.0);
  real_t kH = static_cast<real_t>(0.0);
  real_t px0 = static_cast<real_t>(0.0);
  real_t px1 = static_cast<real_t>(0.0);
  real_t py0 = static_cast<real_t>(0.0);
  real_t py1 = static_cast<real_t>(0.0);

  std::map<std::string, std::vector<real_t>> partials;
  const std::vector<real_t> *xs = nullptr;

  detail::NiceScale xScale;
  detail::NiceScale yScale;

  std::ostringstream svg;

  SvgHistogramRenderer(const correlation::analysis::Histogram &histogram, const PlotConfig &cfg, const HoverInfo &hov,
                       const std::map<std::string, real_t> &component_weights)
      : hist(&histogram), config(&cfg), hover(&hov), weights(&component_weights), xs(&histogram.bins) {
    initializeLayoutAndLabels();
    filterPartials();
  }

  void initializeLayoutAndLabels() {
    title = hist->title.empty() ? "Histogram" : hist->title;
    x_label = hist->x_label.empty() ? "x" : hist->x_label;
    y_label = hist->y_label.empty() ? "y" : hist->y_label;

    if (!hist->x_unit.empty()) {
      x_label += std::format(" ({})", hist->x_unit);
    }
    if (!hist->y_unit.empty()) {
      std::string y_unit = hist->y_unit;
      if (y_unit == "Å^-1") {
        y_unit = "Å⁻¹";
      }
      y_label += std::format(" ({})", y_unit);
    }

    kW = config->effective_width();
    kH = config->effective_height();
    const real_t kLeft = static_cast<real_t>(100.0);
    const real_t kRight = static_cast<real_t>(40.0);
    const real_t kTop = static_cast<real_t>(50.0);
    const real_t kBot = static_cast<real_t>(90.0);

    px0 = kLeft;
    px1 = kW - kRight;
    py0 = kTop;
    py1 = kH - kBot;
  }

  void filterPartials() {
    const auto &raw_partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;
    auto total_it = raw_partials.find("Total");
    if (total_it != raw_partials.end()) {
      partials["Total"] = total_it->second;
    }

    std::vector<std::pair<std::string, real_t>> candidates;
    for (const auto &[key, value] : raw_partials) {
      if (key == "Total" || key.starts_with("Frequency_")) {
        continue;
      }
      real_t score = static_cast<real_t>(0.0);
      if (weights != nullptr && !weights->empty()) {
        auto wit = weights->find(key);
        if (wit != weights->end()) {
          score = wit->second;
        }
      } else {
        for (real_t val : value) {
          score += std::abs(val);
        }
      }
      candidates.emplace_back(key, score);
    }

    std::ranges::sort(candidates, [](const auto &lhs, const auto &rhs) { return lhs.second > rhs.second; });

    std::size_t limit = std::min(candidates.size(), std::size_t(10));
    for (std::size_t i = 0; i < limit; ++i) {
      const std::string &key = candidates.at(i).first;
      auto iter = raw_partials.find(key);
      if (iter != raw_partials.end()) {
        partials[key] = iter->second;
      }
    }
  }

  bool hasNoData() const { return xs == nullptr || xs->empty() || partials.empty(); }

  std::string renderNoData() const {
    if (config->use_native_text) {
      return std::format(
          "<svg width='{0:.0f}' height='{1:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {0:.0f} {1:.0f}\">"
          "<rect width=\"100%\" height=\"100%\" fill=\"{2}\"/>"
          "<text x=\"{3:.1f}\" y=\"{4:.1f}\" font-family=\"'Outfit', 'Plus Jakarta Sans', 'Inter', 'Roboto', 'Helvetica Neue', sans-serif\" font-size=\"{5:.1f}\" text-anchor=\"middle\" fill=\"{6}\">No data available</text></svg>",
          kW, kH, config->bg_color(), kW / static_cast<real_t>(2.0),
          kH / static_cast<real_t>(2.0) + static_cast<real_t>(8.0), static_cast<real_t>(24.0) * config->font_scale,
          config->text_color());
    }
    std::string no_data_path = Roboto::instance().render("No data available", kW / static_cast<real_t>(2.0),
                                                         kH / static_cast<real_t>(2.0) + static_cast<real_t>(8.0),
                                                         static_cast<real_t>(24.0) * config->font_scale, "middle");
    return std::format(
        "<svg width='{0:.0f}' height='{1:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {0:.0f} {1:.0f}\">"
        "<rect width=\"100%\" height=\"100%\" fill=\"{2}\"/>"
        "<path d=\"{3}\" fill=\"{4}\" fill-rule=\"evenodd\" stroke=\"none\"/></svg>",
        kW, kH, config->bg_color(), no_data_path, config->text_color());
  }

  void computeScales() {
    if (xs == nullptr || xs->empty()) {
      return;
    }
    real_t raw_x_min = xs->front();
    real_t raw_x_max = xs->back();
    real_t raw_y_min = static_cast<real_t>(0.0);
    real_t raw_y_max = static_cast<real_t>(0.0);
    for (const auto &[key, value] : partials) {
      for (real_t val : value) {
        raw_y_max = std::max(raw_y_max, val);
        raw_y_min = std::min(raw_y_min, val);
      }
    }
    real_t y_padding = (raw_y_max - raw_y_min) * static_cast<real_t>(0.05);
    raw_y_max += y_padding;

    xScale = detail::NiceScale(detail::DataRange{.min = raw_x_min, .max = raw_x_max}, 11);
    yScale = detail::NiceScale(detail::DataRange{.min = raw_y_min, .max = raw_y_max}, 8);
  }

  void writeHeader() {
    svg << std::format(
        "<svg width='{:.0f}' height='{:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\" shape-rendering=\"geometricPrecision\" text-rendering=\"geometricPrecision\">\n",
        kW, kH, kW, kH);
    svg << "  <defs>\n"
        << "    <filter id=\"tooltip-shadow\" x=\"-10%\" y=\"-10%\" width=\"120%\" height=\"120%\">\n"
        << "      <feDropShadow dx=\"2\" dy=\"4\" stdDeviation=\"4\" flood-color=\"#000000\" flood-opacity=\"0.15\"/>\n"
        << "    </filter>\n";

    std::size_t color_idx = 0;
    for (const auto &[key, value] : partials) {
      std::string col = color(color_idx, config->palette);
      std::string grad_id = std::format("area-grad-{}", color_idx);
      svg << std::format("    <linearGradient id=\"{}\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">\n"
                         "      <stop offset=\"0%\" stop-color=\"{}\" stop-opacity=\"0.35\"/>\n"
                         "      <stop offset=\"100%\" stop-color=\"{}\" stop-opacity=\"0.0\"/>\n"
                         "    </linearGradient>\n",
                         grad_id, col, col);
      color_idx++;
    }
    svg << "  </defs>\n";
    svg << std::format("  <rect width=\"100%\" height=\"100%\" fill=\"{}\" rx=\"6\"/>\n", config->bg_color());
  }

  void drawGridAndAxes() {
    svg << "  <!-- Grid & Ticks -->\n";
    for (real_t y_val : yScale.ticks) {
      real_t spy = mapValue(y_val, yScale.min, yScale.max, py1, py0);
      if (config->show_grid) {
        svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"0.8\" stroke-dasharray=\"3,3\"/>\n",
                           px0, spy, px1, spy, config->grid_color());
      }
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                         px0 - static_cast<real_t>(8.0), spy, px0, spy, config->axis_color());
      svg << renderTextAsPath(fmtScientific(y_val), px0 - static_cast<real_t>(15.0), spy + static_cast<real_t>(7.0),
                              static_cast<real_t>(20.0) * config->font_scale, TextAnchor::End, config->text_color(),
                              config->use_native_text);
    }

    for (real_t x_val : xScale.ticks) {
      real_t spx = mapValue(x_val, xScale.min, xScale.max, px0, px1);
      if (config->show_grid) {
        svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"0.8\" stroke-dasharray=\"3,3\"/>\n",
                           spx, py0, spx, py1, config->grid_color());
      }
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                         spx, py1, spx, py1 + static_cast<real_t>(8.0), config->axis_color());
      svg << renderTextAsPath(fmtScientific(x_val), spx, py1 + static_cast<real_t>(25.0),
                              static_cast<real_t>(20.0) * config->font_scale, TextAnchor::Middle, config->text_color(),
                              config->use_native_text);
    }

    svg << std::format("  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" "
                       "fill=\"none\" stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                       px0, py0, px1 - px0, py1 - py0, config->axis_color());
  }

  void drawEmphasisLines() {
    std::vector<real_t> emphasis_values;
    if (title.contains("g(r)")) {
      emphasis_values.push_back(static_cast<real_t>(1.0));
    } else if (title.contains("G(r)")) {
      emphasis_values.push_back(static_cast<real_t>(0.0));
    } else if (title.contains("S(Q)") || title.contains("S(q)")) {
      emphasis_values.push_back(static_cast<real_t>(0.0));
      emphasis_values.push_back(static_cast<real_t>(1.0));
    }

    for (real_t focus_y : emphasis_values) {
      if (focus_y >= yScale.min && focus_y <= yScale.max) {
        real_t spy = mapValue(focus_y, yScale.min, yScale.max, py1, py0);
        std::string extra = (std::abs(focus_y - static_cast<real_t>(1.0)) < static_cast<real_t>(1e-6))
                                ? " stroke-dasharray=\"5,5\""
                                : "";
        svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"1.8\"{}/>\n",
                           px0, spy, px1, spy, config->axis_color(), extra);
      }
    }
  }

  void drawAreaFills() {
    if (config == nullptr || !config->fill_area || xs == nullptr) {
      return;
    }
    std::size_t fill_ci = 0;
    for (const auto &[key, value] : partials) {
      std::string grad_id = std::format("area-grad-{}", fill_ci);
      std::size_t num_points = std::min(xs->size(), value.size());
      if (num_points > 1) {
        svg << std::format("  <polygon fill=\"url(#{})\" stroke=\"none\" points=\"", grad_id);
        real_t sx_start = mapValue(xs->at(0), xScale.min, xScale.max, px0, px1);
        svg << std::format("{:.2f},{:.2f} ", sx_start, py1);
        for (std::size_t point_idx = 0; point_idx < num_points; ++point_idx) {
          real_t screen_x = mapValue(xs->at(point_idx), xScale.min, xScale.max, px0, px1);
          real_t screen_y = mapValue(value.at(point_idx), yScale.min, yScale.max, py1, py0);
          real_t sy_clamped = std::min(screen_y, py1);
          svg << std::format("{:.2f},{:.2f} ", screen_x, sy_clamped);
        }
        real_t sx_end = mapValue(xs->at(num_points - 1), xScale.min, xScale.max, px0, px1);
        svg << std::format("{:.2f},{:.2f}", sx_end, py1);
        svg << "\" />\n";
      }
      fill_ci++;
    }
  }

  std::vector<std::pair<std::string, std::string>> drawPolylines() {
    std::vector<std::pair<std::string, std::string>> legend;
    if (xs == nullptr) {
      return legend;
    }
    std::size_t color_idx = 0;
    for (const auto &[key, value] : partials) {
      const std::string col = color(color_idx++, config->palette);
      svg << std::format(
          R"(  <polyline fill="none" stroke="{}" stroke-width="{:.1f}" stroke-linejoin="round" points=")", col,
          config->line_width);
      std::size_t num_points = std::min(xs->size(), value.size());
      for (std::size_t point_idx = 0; point_idx < num_points; ++point_idx) {
        real_t screen_x = mapValue(xs->at(point_idx), xScale.min, xScale.max, px0, px1);
        real_t screen_y = mapValue(value.at(point_idx), yScale.min, yScale.max, py1, py0);
        svg << std::format("{:.2f},{:.2f} ", screen_x, screen_y);
      }
      svg << "\" />\n";
      legend.emplace_back(key, col);
    }
    return legend;
  }

  void drawMarkers() {
    if (!config->show_markers || xs == nullptr) {
      return;
    }
    std::size_t marker_ci = 0;
    for (const auto &[key, value] : partials) {
      const std::string col = color(marker_ci++, config->palette);
      std::size_t num_points = std::min(xs->size(), value.size());
      for (std::size_t point_idx = 0; point_idx < num_points; ++point_idx) {
        real_t screen_x = mapValue(xs->at(point_idx), xScale.min, xScale.max, px0, px1);
        real_t screen_y = mapValue(value.at(point_idx), yScale.min, yScale.max, py1, py0);
        svg << std::format("  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"3.5\" fill=\"{}\" stroke=\"none\"/>\n", screen_x,
                           screen_y, col);
      }
    }
  }

  void drawLegend(const std::vector<std::pair<std::string, std::string>> &legend) {
    if (!config->show_legend) {
      return;
    }
    real_t legend_x = px1 - static_cast<real_t>(15.0);
    real_t legend_y = py0 + static_cast<real_t>(25.0);
    for (const auto &iter : std::views::reverse(legend)) {
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"4.0\"/>\n",
                         legend_x - static_cast<real_t>(40.0), legend_y, legend_x - static_cast<real_t>(10.0), legend_y,
                         iter.second);
      svg << renderTextAsPath(iter.first, legend_x - static_cast<real_t>(45.0), legend_y + static_cast<real_t>(6.0),
                              static_cast<real_t>(18.0) * config->font_scale, TextAnchor::End, config->text_color(),
                              config->use_native_text);
      legend_y += static_cast<real_t>(28.0);
    }
  }

  void drawTitlesAndLabels() {
    svg << renderTextAsPath(x_label, (px0 + px1) / static_cast<real_t>(2.0), py1 + static_cast<real_t>(75.0),
                            static_cast<real_t>(28.0) * config->font_scale, TextAnchor::Middle, config->text_color(),
                            config->use_native_text);

    svg << std::format("  <g transform=\"translate({:.1f}, {:.1f}) rotate(-90)\">\n", static_cast<real_t>(40.0),
                       (py0 + py1) / static_cast<real_t>(2.0));
    svg << renderTextAsPath(y_label, static_cast<real_t>(0.0), static_cast<real_t>(0.0),
                            static_cast<real_t>(28.0) * config->font_scale, TextAnchor::Middle, config->text_color(),
                            config->use_native_text);
    svg << "  </g>\n";
  }

  struct NearestPoint {
    std::size_t index = 0;
    std::string key;
  };

  NearestPoint findNearestPoint(real_t mouse_sx, real_t mouse_sy) const {
    real_t min_dist_sq = static_cast<real_t>(1e30);
    NearestPoint best;
    for (const auto &[key, value] : partials) {
      std::size_t num_points = std::min(xs->size(), value.size());
      for (std::size_t point_idx = 0; point_idx < num_points; ++point_idx) {
        real_t sx_i = mapValue(xs->at(point_idx), xScale.min, xScale.max, px0, px1);
        real_t sy_i = mapValue(value.at(point_idx), yScale.min, yScale.max, py1, py0);
        real_t dx_i = sx_i - mouse_sx;
        real_t dy_i = sy_i - mouse_sy;
        real_t dist_sq = dx_i * dx_i + dy_i * dy_i;
        if (dist_sq < min_dist_sq) {
          min_dist_sq = dist_sq;
          best.index = point_idx;
          best.key = key;
        }
      }
    }
    return best;
  }

  struct TooltipPosition {
    real_t sx_data;
    real_t snapped_sy_data;
    real_t target_x;
  };

  void drawTooltipBox(const TooltipPosition &pos, const std::string &best_key,
                      const std::vector<std::tuple<std::string, real_t, std::string>> &hover_values) {
    if (hover_values.empty()) {
      return;
    }
    real_t tooltip_w = static_cast<real_t>(200.0);
    real_t tooltip_h = static_cast<real_t>(35.0) + static_cast<real_t>(22.0) * static_cast<real_t>(hover_values.size());

    // Tooltip position (flip sides depending on cursor location)
    real_t tooltip_x = (pos.sx_data < kW / static_cast<real_t>(2.0))
                           ? pos.sx_data + static_cast<real_t>(15.0)
                           : pos.sx_data - tooltip_w - static_cast<real_t>(15.0);
    real_t tooltip_y = (pos.snapped_sy_data >= static_cast<real_t>(0.0))
                           ? pos.snapped_sy_data - tooltip_h / static_cast<real_t>(2.0)
                           : py0 + static_cast<real_t>(15.0);
    tooltip_y =
        std::max(py0 + static_cast<real_t>(8.0), std::min(py1 - tooltip_h - static_cast<real_t>(8.0), tooltip_y));

    std::string card_bg = (config->theme == PlotConfig::Theme::Light) ? "#FFFFFF" : "#181825";
    std::string card_border = (config->theme == PlotConfig::Theme::Light) ? "#dddddd" : "#45475a";
    std::string fg_color = (config->theme == PlotConfig::Theme::Light) ? "#333333" : "#cdd6f4";

    svg << std::format(
        "  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" rx=\"6\" "
        "fill=\"{}\" fill-opacity=\"0.92\" stroke=\"{}\" stroke-width=\"1.5\" filter=\"url(#tooltip-shadow)\"/>\n",
        tooltip_x, tooltip_y, tooltip_w, tooltip_h, card_bg, card_border);

    // Header: x value (with unit or pure label)
    std::string x_unit_str = hist->x_unit;
    std::string header_txt =
        std::format("{} = {:.4f}{}", x_label, pos.target_x, x_unit_str.empty() ? "" : " " + x_unit_str);
    // Clean up title (remove parenthesis unit if present, e.g. "r (Å)" to "r")
    auto paren = x_label.find(" (");
    if (paren != std::string::npos) {
      header_txt = std::format("{} = {:.4f}", x_label.substr(0, paren), pos.target_x);
    }

    svg << renderTextAsPath(header_txt, tooltip_x + static_cast<real_t>(12.0), tooltip_y + static_cast<real_t>(20.0),
                            static_cast<real_t>(14.0) * config->font_scale, TextAnchor::Start, fg_color,
                            config->use_native_text);

    real_t cur_y = tooltip_y + static_cast<real_t>(42.0);
    for (const auto &[name, val, col] : hover_values) {
      // Color swatch dot
      svg << std::format("  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"5\" fill=\"{}\"/>\n",
                         tooltip_x + static_cast<real_t>(18.0), cur_y - static_cast<real_t>(4.0), col);
      // Label and value
      std::string line_txt = std::format("{}: {:.4f}", name, val);
      if (name == best_key) {
        line_txt += " (nearest)";
      }
      svg << renderTextAsPath(line_txt, tooltip_x + static_cast<real_t>(30.0), cur_y,
                              static_cast<real_t>(13.0) * config->font_scale, TextAnchor::Start, fg_color,
                              config->use_native_text);
      cur_y += static_cast<real_t>(22.0);
    }
  }

  void drawHoverTooltip() {
    if (hover == nullptr || !hover->active || hover->widget_width <= static_cast<real_t>(0.0) ||
        hover->widget_height <= static_cast<real_t>(0.0) || xs == nullptr) {
      return;
    }
    real_t svg_aspect = kW / kH;
    real_t widget_aspect = hover->widget_width / hover->widget_height;
    real_t scale = static_cast<real_t>(1.0);
    real_t offset_x = static_cast<real_t>(0.0);
    real_t offset_y = static_cast<real_t>(0.0);

    if (widget_aspect > svg_aspect) {
      scale = hover->widget_height / kH;
      offset_x = (hover->widget_width - kW * scale) / static_cast<real_t>(2.0);
    } else {
      scale = hover->widget_width / kW;
      offset_y = (hover->widget_height - kH * scale) / static_cast<real_t>(2.0);
    }

    real_t mouse_sx = (hover->mouse_x - offset_x) / scale;
    real_t mouse_sy = (hover->mouse_y - offset_y) / scale;

    if (mouse_sx < px0 || mouse_sx > px1 || mouse_sy < py0 || mouse_sy > py1) {
      return;
    }

    NearestPoint best = findNearestPoint(mouse_sx, mouse_sy);

    std::size_t idx = best.index;
    real_t target_x = xs->at(idx);
    real_t sx_data = mapValue(target_x, xScale.min, xScale.max, px0, px1);

    // Draw vertical guide line
    svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"1.5\" stroke-dasharray=\"4,4\"/>\n",
                       sx_data, py0, sx_data, py1, config->axis_color());

    // Collect values and draw markers on curves
    std::vector<std::tuple<std::string, real_t, std::string>> hover_values; // name, value, color
    std::size_t color_idx = 0;
    real_t snapped_sy_data = static_cast<real_t>(-1.0);
    for (const auto &[key, value] : partials) {
      if (idx < value.size()) {
        const std::string col = color(color_idx++, config->palette);
        real_t y_val = value.at(idx);
        real_t sy_data = mapValue(y_val, yScale.min, yScale.max, py1, py0);

        if (key == best.key) {
          snapped_sy_data = sy_data;
        }

        // Bullet marker
        svg << std::format(
            "  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"6\" fill=\"{}\" stroke=\"{}\" stroke-width=\"2\"/>\n", sx_data,
            sy_data, col, config->bg_color());

        hover_values.emplace_back(key, y_val, col);
      }
    }

    // Draw horizontal guide line to the snapped data point of the closest curve
    if (snapped_sy_data >= py0 && snapped_sy_data <= py1) {
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.5\" stroke-dasharray=\"4,4\"/>\n",
                         px0, snapped_sy_data, sx_data, snapped_sy_data, config->axis_color());
    }

    drawTooltipBox({sx_data, snapped_sy_data, target_x}, best.key, hover_values);
  }

  std::string getResult() {
    svg << "</svg>\n";
    return svg.str();
  }
};

struct SvgComparisonRenderer {
  const std::vector<LabeledHistogram> *datasets = nullptr;
  const std::string *partial_key = nullptr;
  const PlotConfig *config = nullptr;
  const HoverInfo *hover = nullptr;

  real_t kWidth = static_cast<real_t>(0.0);
  real_t kHeight = static_cast<real_t>(0.0);
  real_t px0 = static_cast<real_t>(0.0);
  real_t px1 = static_cast<real_t>(0.0);
  real_t py0 = static_cast<real_t>(0.0);
  real_t py1 = static_cast<real_t>(0.0);

  real_t raw_x_min = static_cast<real_t>(1e30);
  real_t raw_x_max = static_cast<real_t>(-1e30);
  real_t raw_y_min = static_cast<real_t>(0.0);
  real_t raw_y_max = static_cast<real_t>(-1e30);

  detail::NiceScale xScale;
  detail::NiceScale yScale;

  std::string x_label;
  std::string y_label;

  std::ostringstream svg;
  std::vector<std::pair<std::string, std::string>> legend;

  SvgComparisonRenderer(const std::vector<LabeledHistogram> &datasets_ref, const std::string &partial_key_ref,
                        const PlotConfig &config_ref, const HoverInfo &hover_ref)
      : datasets(&datasets_ref), partial_key(&partial_key_ref), config(&config_ref), hover(&hover_ref) {}

  void initializeLabelsAndLayout() {
    const auto &ref = *datasets->front().hist;
    x_label = ref.x_label.empty() ? "x" : ref.x_label;
    y_label = ref.y_label.empty() ? "y" : ref.y_label;

    if (!ref.x_unit.empty()) {
      x_label += std::format(" ({})", ref.x_unit);
    }
    if (!ref.y_unit.empty()) {
      std::string y_unit = ref.y_unit;
      if (y_unit == "Å^-1") {
        y_unit = "Å⁻¹";
      }
      y_label += std::format(" ({})", y_unit);
    }

    kWidth = config->effective_width();
    kHeight = config->effective_height();
    const real_t kLeft = static_cast<real_t>(100.0);
    const real_t kRight = static_cast<real_t>(40.0);
    const real_t kTop = static_cast<real_t>(50.0);
    const real_t kBot = static_cast<real_t>(90.0);

    px0 = kLeft;
    px1 = kWidth - kRight;
    py0 = kTop;
    py1 = kHeight - kBot;
  }

  void computeGlobalRanges() {
    for (const auto &dataset : *datasets) {
      const auto &partials =
          dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
      auto iterator = partials.find(*partial_key);
      if (iterator == partials.end()) {
        continue;
      }
      const auto &x_values = dataset.hist->bins;
      const auto &y_values = iterator->second;
      if (!x_values.empty()) {
        raw_x_min = std::min(raw_x_min, x_values.front());
        raw_x_max = std::max(raw_x_max, x_values.back());
      }
      for (real_t y_value : y_values) {
        raw_y_max = std::max(raw_y_max, y_value);
        raw_y_min = std::min(raw_y_min, y_value);
      }
    }
  }

  bool hasValidRanges() const { return raw_x_max > raw_x_min && raw_y_max > raw_y_min; }

  std::string renderNoComparisonData() const {
    if (config->use_native_text) {
      return std::format(
          "<svg width='{0:.0f}' height='{1:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {0:.0f} {1:.0f}\">"
          "<rect width=\"100%\" height=\"100%\" fill=\"{2}\"/>"
          "<text x=\"{3:.1f}\" y=\"{4:.1f}\" font-family=\"'Outfit', 'Plus Jakarta Sans', 'Inter', 'Roboto', 'Helvetica Neue', sans-serif\" font-size=\"{5:.1f}\" text-anchor=\"middle\" fill=\"{6}\">No comparison data</text></svg>",
          kWidth, kHeight, config->bg_color(), kWidth / static_cast<real_t>(2.0),
          kHeight / static_cast<real_t>(2.0) + static_cast<real_t>(8.0), static_cast<real_t>(24.0) * config->font_scale,
          config->text_color());
    }
    std::string no_data_path = Roboto::instance().render("No comparison data", kWidth / static_cast<real_t>(2.0),
                                                         kHeight / static_cast<real_t>(2.0) + static_cast<real_t>(8.0),
                                                         static_cast<real_t>(24.0) * config->font_scale, "middle");
    return std::format(
        "<svg width='{0:.0f}' height='{1:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {0:.0f} {1:.0f}\">"
        "<rect width=\"100%\" height=\"100%\" fill=\"{2}\"/>"
        "<path d=\"{3}\" fill=\"{4}\" fill-rule=\"evenodd\" stroke=\"none\"/></svg>",
        kWidth, kHeight, config->bg_color(), no_data_path, config->text_color());
  }

  void computeScales() {
    real_t y_padding = (raw_y_max - raw_y_min) * static_cast<real_t>(0.05);
    raw_y_max += y_padding;

    xScale = detail::NiceScale(detail::DataRange{.min = raw_x_min, .max = raw_x_max}, 11);
    yScale = detail::NiceScale(detail::DataRange{.min = raw_y_min, .max = raw_y_max}, 8);
  }

  void writeHeaderAndDefs() {
    svg << std::format(
        "<svg width='{:.0f}' height='{:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\" shape-rendering=\"geometricPrecision\" text-rendering=\"geometricPrecision\">\n",
        kWidth, kHeight, kWidth, kHeight);
    svg << "  <defs>\n"
        << "    <filter id=\"tooltip-shadow\" x=\"-10%\" y=\"-10%\" width=\"120%\" height=\"120%\">\n"
        << "      <feDropShadow dx=\"2\" dy=\"4\" stdDeviation=\"4\" flood-color=\"#000000\" flood-opacity=\"0.15\"/>\n"
        << "    </filter>\n";

    // Create linear gradients for area fills in comparison
    std::size_t color_idx = 0;
    for (const auto &dataset : *datasets) {
      std::string col = detail::color(color_idx, config->palette);
      std::string grad_id = std::format("area-grad-{}", color_idx);
      svg << std::format("    <linearGradient id=\"{}\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">\n"
                         "      <stop offset=\"0%\" stop-color=\"{}\" stop-opacity=\"0.35\"/>\n"
                         "      <stop offset=\"100%\" stop-color=\"{}\" stop-opacity=\"0.0\"/>\n"
                         "    </linearGradient>\n",
                         grad_id, col, col);
      color_idx++;
    }
    svg << "  </defs>\n";
    svg << std::format("  <rect width=\"100%\" height=\"100%\" fill=\"{}\" rx=\"6\"/>\n", config->bg_color());
  }

  void drawGridAndAxes() {
    for (real_t y_value : yScale.ticks) {
      real_t spy = detail::mapValue(y_value, yScale.min, yScale.max, py1, py0);
      if (config->show_grid) {
        svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"0.8\" stroke-dasharray=\"3,3\"/>\n",
                           px0, spy, px1, spy, config->grid_color());
      }
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                         px0 - static_cast<real_t>(8.0), spy, px0, spy, config->axis_color());
      svg << renderTextAsPath(detail::fmtScientific(y_value), px0 - static_cast<real_t>(15.0),
                              spy + static_cast<real_t>(7.0), static_cast<real_t>(20.0) * config->font_scale,
                              TextAnchor::End, config->text_color(), config->use_native_text);
    }

    for (real_t x_value : xScale.ticks) {
      real_t spx = detail::mapValue(x_value, xScale.min, xScale.max, px0, px1);
      if (config->show_grid) {
        svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"0.8\" stroke-dasharray=\"3,3\"/>\n",
                           spx, py0, spx, py1, config->grid_color());
      }
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                         spx, py1, spx, py1 + static_cast<real_t>(8.0), config->axis_color());
      svg << renderTextAsPath(detail::fmtScientific(x_value), spx, py1 + static_cast<real_t>(25.0),
                              static_cast<real_t>(20.0) * config->font_scale, TextAnchor::Middle, config->text_color(),
                              config->use_native_text);
    }

    // Axis border
    svg << std::format("  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" "
                       "fill=\"none\" stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                       px0, py0, px1 - px0, py1 - py0, config->axis_color());
  }

  void drawAreaFills() {
    if (config->fill_area) {
      std::size_t fill_ci = 0;
      for (const auto &dataset : *datasets) {
        const auto &partials =
            dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
        auto iterator = partials.find(*partial_key);
        if (iterator != partials.end()) {
          const auto &x_values = dataset.hist->bins;
          const auto &y_values = iterator->second;
          std::string grad_id = std::format("area-grad-{}", fill_ci);
          std::size_t n_points = std::min(x_values.size(), y_values.size());
          if (n_points > 1) {
            svg << std::format("  <polygon fill=\"url(#{})\" stroke=\"none\" points=\"", grad_id);
            real_t sp_x_start = detail::mapValue(x_values[0], xScale.min, xScale.max, px0, px1);
            svg << std::format("{:.2f},{:.2f} ", sp_x_start, py1);
            for (std::size_t i = 0; i < n_points; ++i) {
              real_t sp_x = detail::mapValue(x_values[i], xScale.min, xScale.max, px0, px1);
              real_t sp_y = detail::mapValue(y_values[i], yScale.min, yScale.max, py1, py0);
              real_t sp_y_clamped = std::min(sp_y, py1);
              svg << std::format("{:.2f},{:.2f} ", sp_x, sp_y_clamped);
            }
            real_t sp_x_end = detail::mapValue(x_values[n_points - 1], xScale.min, xScale.max, px0, px1);
            svg << std::format("{:.2f},{:.2f}", sp_x_end, py1);
            svg << "\" />\n";
          }
        }
        fill_ci++;
      }
    }
  }

  void drawPolylines() {
    std::size_t color_idx = 0;
    for (const auto &dataset : *datasets) {
      const auto &partials =
          dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
      auto partial_iter = partials.find(*partial_key);
      if (partial_iter == partials.end()) {
        continue;
      }

      const auto &x_values = dataset.hist->bins;
      const auto &y_values = partial_iter->second;
      const std::string col = detail::color(color_idx++, config->palette);

      svg << std::format(
          R"(  <polyline fill="none" stroke="{}" stroke-width="{:.1f}" stroke-linejoin="round" points=")", col,
          config->line_width);
      std::size_t n_points = std::min(x_values.size(), y_values.size());
      for (std::size_t i = 0; i < n_points; ++i) {
        real_t sp_x = detail::mapValue(x_values[i], xScale.min, xScale.max, px0, px1);
        real_t sp_y = detail::mapValue(y_values[i], yScale.min, yScale.max, py1, py0);
        svg << std::format("{:.2f},{:.2f} ", sp_x, sp_y);
      }
      svg << "\" />\n";
      legend.emplace_back(dataset.label, col);
    }
  }

  void drawMarkers() {
    if (config->show_markers) {
      std::size_t marker_color_idx = 0;
      for (const auto &dataset : *datasets) {
        const auto &partials =
            dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
        auto partial_iter = partials.find(*partial_key);
        if (partial_iter != partials.end()) {
          const auto &x_values = dataset.hist->bins;
          const auto &y_values = partial_iter->second;
          const std::string col = detail::color(marker_color_idx++, config->palette);
          std::size_t n_points = std::min(x_values.size(), y_values.size());
          for (std::size_t i = 0; i < n_points; ++i) {
            real_t sp_x = detail::mapValue(x_values[i], xScale.min, xScale.max, px0, px1);
            real_t sp_y = detail::mapValue(y_values[i], yScale.min, yScale.max, py1, py0);
            svg << std::format("  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"3.5\" fill=\"{}\" stroke=\"none\"/>\n", sp_x,
                               sp_y, col);
          }
        } else {
          marker_color_idx++;
        }
      }
    }
  }

  void drawLegend() {
    if (config->show_legend) {
      real_t legend_x = px1 - static_cast<real_t>(15.0);
      real_t legend_y = py0 + static_cast<real_t>(25.0);
      for (auto [label, color] : std::views::reverse(legend)) {
        svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"4.0\"/>\n",
                           legend_x - static_cast<real_t>(40.0), legend_y, legend_x - static_cast<real_t>(10.0),
                           legend_y, color);
        svg << renderTextAsPath(label, legend_x - static_cast<real_t>(45.0), legend_y + static_cast<real_t>(6.0),
                                static_cast<real_t>(18.0) * config->font_scale, TextAnchor::End, config->text_color(),
                                config->use_native_text);
        legend_y += static_cast<real_t>(28.0);
      }
    }
  }

  void drawTitlesAndLabels() {
    // X-axis label
    svg << renderTextAsPath(x_label, (px0 + px1) / static_cast<real_t>(2.0), py1 + static_cast<real_t>(75.0),
                            static_cast<real_t>(28.0) * config->font_scale, TextAnchor::Middle, config->text_color(),
                            config->use_native_text);

    // Y-axis label (rotated)
    svg << std::format("  <g transform=\"translate({:.1f}, {:.1f}) rotate(-90)\">\n", static_cast<real_t>(40.0),
                       (py0 + py1) / static_cast<real_t>(2.0));
    svg << renderTextAsPath(y_label, static_cast<real_t>(0.0), static_cast<real_t>(0.0),
                            static_cast<real_t>(28.0) * config->font_scale, TextAnchor::Middle, config->text_color(),
                            config->use_native_text);
    svg << "  </g>\n";
  }

  struct NearestPoint {
    std::size_t index = 0;
    std::string label;
  };

  NearestPoint findNearestPoint(real_t sp_x, real_t sp_y) const {
    real_t min_dist_sq = static_cast<real_t>(1e30);
    NearestPoint best;

    for (const auto &dataset : *datasets) {
      const auto &partials =
          dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
      auto pit = partials.find(*partial_key);
      if (pit == partials.end()) {
        continue;
      }
      const auto &x_values = dataset.hist->bins;
      const auto &y_values = pit->second;
      std::size_t n_points = std::min(x_values.size(), y_values.size());
      for (std::size_t i = 0; i < n_points; ++i) {
        real_t sp_x_i = detail::mapValue(x_values[i], xScale.min, xScale.max, px0, px1);
        real_t sp_y_i = detail::mapValue(y_values[i], yScale.min, yScale.max, py1, py0);
        real_t dx_i = sp_x_i - sp_x;
        real_t dy_i = sp_y_i - sp_y;
        real_t dist_sq = dx_i * dx_i + dy_i * dy_i;
        if (dist_sq < min_dist_sq) {
          min_dist_sq = dist_sq;
          best.index = i;
          best.label = dataset.label;
        }
      }
    }
    return best;
  }

  std::vector<std::tuple<std::string, real_t, std::string>>
  collectHoverValuesAndDrawMarkers(std::size_t idx, const std::string &best_label, real_t sx_data,
                                   real_t &snapped_sy_data) {
    std::vector<std::tuple<std::string, real_t, std::string>> hover_values;
    std::size_t color_idx = 0;
    for (const auto &dataset : *datasets) {
      const auto &partials =
          dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
      auto pit = partials.find(*partial_key);
      if (pit != partials.end() && idx < pit->second.size()) {
        const std::string col = detail::color(color_idx++, config->palette);
        real_t y_val = pit->second[idx];
        real_t sy_data = detail::mapValue(y_val, yScale.min, yScale.max, py1, py0);

        if (dataset.label == best_label) {
          snapped_sy_data = sy_data;
        }

        // Bullet marker
        svg << std::format(
            "  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"6\" fill=\"{}\" stroke=\"{}\" stroke-width=\"2\"/>\n", sx_data,
            sy_data, col, config->bg_color());

        hover_values.emplace_back(dataset.label, y_val, col);
      } else {
        color_idx++;
      }
    }
    return hover_values;
  }

  struct TooltipPosition {
    real_t sx_data = static_cast<real_t>(0.0);
    real_t snapped_sy_data = static_cast<real_t>(0.0);
    real_t target_x = static_cast<real_t>(0.0);
  };

  void drawTooltipBox(const TooltipPosition &pos, const std::string &best_label,
                      const std::vector<std::tuple<std::string, real_t, std::string>> &hover_values) {
    if (hover_values.empty()) {
      return;
    }
    real_t tooltip_w = static_cast<real_t>(200.0);
    real_t tooltip_h = static_cast<real_t>(35.0) + static_cast<real_t>(22.0) * static_cast<real_t>(hover_values.size());

    real_t tooltip_x = (pos.sx_data < kWidth / static_cast<real_t>(2.0))
                           ? pos.sx_data + static_cast<real_t>(15.0)
                           : pos.sx_data - tooltip_w - static_cast<real_t>(15.0);
    real_t tooltip_y = (pos.snapped_sy_data >= static_cast<real_t>(0.0))
                           ? pos.snapped_sy_data - tooltip_h / static_cast<real_t>(2.0)
                           : py0 + static_cast<real_t>(15.0);
    tooltip_y =
        std::max(py0 + static_cast<real_t>(8.0), std::min(py1 - tooltip_h - static_cast<real_t>(8.0), tooltip_y));

    std::string card_bg = (config->theme == PlotConfig::Theme::Light) ? "#FFFFFF" : "#181825";
    std::string card_border = (config->theme == PlotConfig::Theme::Light) ? "#dddddd" : "#45475a";
    std::string fg_color = (config->theme == PlotConfig::Theme::Light) ? "#333333" : "#cdd6f4";

    svg << std::format(
        "  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" rx=\"6\" "
        "fill=\"{}\" fill-opacity=\"0.92\" stroke=\"{}\" stroke-width=\"1.5\" filter=\"url(#tooltip-shadow)\"/>\n",
        tooltip_x, tooltip_y, tooltip_w, tooltip_h, card_bg, card_border);

    // Header
    const auto &ref = *datasets->front().hist;
    std::string x_unit_str = ref.x_unit;
    std::string header_txt =
        std::format("{} = {:.4f}{}", x_label, pos.target_x, x_unit_str.empty() ? "" : " " + x_unit_str);
    auto paren = x_label.find(" (");
    if (paren != std::string::npos) {
      header_txt = std::format("{} = {:.4f}", x_label.substr(0, paren), pos.target_x);
    }

    svg << renderTextAsPath(header_txt, tooltip_x + static_cast<real_t>(12.0), tooltip_y + static_cast<real_t>(24.0),
                            static_cast<real_t>(14.0) * config->font_scale, TextAnchor::Start, fg_color,
                            config->use_native_text);

    real_t cur_y = tooltip_y + static_cast<real_t>(46.0);
    for (const auto &[name, val, col] : hover_values) {
      svg << std::format("  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"5\" fill=\"{}\"/>\n",
                         tooltip_x + static_cast<real_t>(18.0), cur_y - static_cast<real_t>(5.0), col);
      std::string line_txt = std::format("{}: {:.4f}", name, val);
      if (name == best_label) {
        line_txt += " (nearest)";
      }
      svg << renderTextAsPath(line_txt, tooltip_x + static_cast<real_t>(30.0), cur_y,
                              static_cast<real_t>(13.0) * config->font_scale, TextAnchor::Start, fg_color,
                              config->use_native_text);
      cur_y += static_cast<real_t>(22.0);
    }
  }

  void drawHoverTooltip() {
    if (!hover->active || hover->widget_width <= static_cast<real_t>(0.0) ||
        hover->widget_height <= static_cast<real_t>(0.0)) {
      return;
    }
    real_t svg_aspect = kWidth / kHeight;
    real_t widget_aspect = hover->widget_width / hover->widget_height;
    real_t scale = static_cast<real_t>(1.0);
    real_t offset_x = static_cast<real_t>(0.0);
    real_t offset_y = static_cast<real_t>(0.0);

    if (widget_aspect > svg_aspect) {
      scale = hover->widget_height / kHeight;
      offset_x = (hover->widget_width - kWidth * scale) / static_cast<real_t>(2.0);
    } else {
      scale = hover->widget_width / kWidth;
      offset_y = (hover->widget_height - kHeight * scale) / static_cast<real_t>(2.0);
    }

    real_t sp_x = (hover->mouse_x - offset_x) / scale;
    real_t sp_y = (hover->mouse_y - offset_y) / scale;

    if (sp_x < px0 || sp_x > px1 || sp_y < py0 || sp_y > py1) {
      return;
    }

    auto [best_idx, best_label] = findNearestPoint(sp_x, sp_y);

    const auto &x_bins = datasets->front().hist->bins;
    if (x_bins.empty() || best_idx >= x_bins.size()) {
      return;
    }

    real_t target_x = x_bins[best_idx];
    real_t sx_data = detail::mapValue(target_x, xScale.min, xScale.max, px0, px1);

    // Draw vertical guide line
    svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"1.5\" stroke-dasharray=\"4,4\"/>\n",
                       sx_data, py0, sx_data, py1, config->axis_color());

    real_t snapped_sy_data = static_cast<real_t>(-1.0);
    auto hover_values = collectHoverValuesAndDrawMarkers(best_idx, best_label, sx_data, snapped_sy_data);

    // Draw horizontal guide line to the snapped data point of the closest dataset
    if (snapped_sy_data >= py0 && snapped_sy_data <= py1) {
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.5\" stroke-dasharray=\"4,4\"/>\n",
                         px0, snapped_sy_data, sx_data, snapped_sy_data, config->axis_color());
    }

    drawTooltipBox({sx_data, snapped_sy_data, target_x}, best_label, hover_values);
  }

  std::string getResult() {
    svg << "</svg>\n";
    return svg.str();
  }
};

} // namespace detail

// ============================================================================
// Public API
// ============================================================================

/**
 * @brief Renders text as a filled SVG path using the Roboto outline font.
 * Uses evenodd fill-rule to render the font's inner holes properly.
 */
inline std::string renderTextAsPath(const std::string &text, real_t x_pos, real_t y_pos, real_t size, TextAnchor anchor,
                                    const std::string &color, bool use_native_text = false) {
  std::string anchor_str = "start";
  if (anchor == TextAnchor::Middle) {
    anchor_str = "middle";
  } else if (anchor == TextAnchor::End) {
    anchor_str = "end";
  }

  if (use_native_text) {
    return std::format(
        "  <text x=\"{:.1f}\" y=\"{:.1f}\" font-family=\"'Outfit', 'Plus Jakarta Sans', 'Inter', 'Roboto', 'Helvetica Neue', sans-serif\" font-size=\"{:.1f}\" text-anchor=\"{}\" fill=\"{}\">{}</text>\n",
        x_pos, y_pos, size, anchor_str, color, text);
  }
  std::string path_d = Roboto::instance().render(text, static_cast<double>(x_pos), static_cast<double>(y_pos),
                                                 static_cast<double>(size), anchor_str);
  return std::format("  <path d=\"{}\" fill=\"{}\" fill-rule=\"evenodd\" stroke=\"none\"/>\n", path_d, color);
}

/**
 * @brief Renders a `Histogram` as a self-contained SVG string.
 *
 * @param hist     The Histogram to render.
 * @param config   Optional plot configuration (theme, size, etc.).
 * @param hover    Optional hover interaction info.
 * @param weights  Optional weights for each partial component to prioritize rendering.
 * @returns        A complete SVG document as `std::string`.
 */
inline std::string renderHistogramAsSvg(const correlation::analysis::Histogram &hist, const PlotConfig &config = {},
                                        const HoverInfo &hover = {},
                                        const std::map<std::string, real_t> &weights = {}) {
  detail::SvgHistogramRenderer renderer(hist, config, hover, weights);
  if (renderer.hasNoData()) {
    return renderer.renderNoData();
  }
  renderer.computeScales();
  renderer.writeHeader();
  renderer.drawGridAndAxes();
  renderer.drawEmphasisLines();
  renderer.drawAreaFills();
  auto legend = renderer.drawPolylines();
  renderer.drawMarkers();
  renderer.drawLegend(legend);
  renderer.drawTitlesAndLabels();
  renderer.drawHoverTooltip();
  return renderer.getResult();
}

// ============================================================================
// Comparison overlay renderer
// ============================================================================

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
                                       const std::string &partial_key = "Total", const PlotConfig &config = {},
                                       const HoverInfo &hover = {}) {
  if (datasets.empty()) {
    return "";
  }
  detail::SvgComparisonRenderer renderer(datasets, partial_key, config, hover);
  renderer.initializeLabelsAndLayout();
  renderer.computeGlobalRanges();
  if (!renderer.hasValidRanges()) {
    return renderer.renderNoComparisonData();
  }
  renderer.computeScales();
  renderer.writeHeaderAndDefs();
  renderer.drawGridAndAxes();
  renderer.drawAreaFills();
  renderer.drawPolylines();
  renderer.drawMarkers();
  renderer.drawLegend();
  renderer.drawTitlesAndLabels();
  renderer.drawHoverTooltip();
  return renderer.getResult();
}

} // namespace correlation::plotters
