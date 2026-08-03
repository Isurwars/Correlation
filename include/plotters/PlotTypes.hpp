/**
 * @file PlotTypes.hpp
 * @brief Common types, configuration structures, and palette/formatting utilities for SVG plotting.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "math/Precision.hpp"
#include "plotters/PathFont.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <span>
#include <string>
#include <string_view>
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
  bool show_difference_curve = false;         ///< Toggle overlaid absolute difference curve Y_diff = Y_ref - Y_target.

  // Publication settings
  real_t font_scale = static_cast<real_t>(1.0); ///< Multiplier for all font sizes
  real_t line_width = static_cast<real_t>(3.0); ///< Data line stroke width
  bool show_legend = true;                      ///< Toggle legend visibility
  bool use_native_text = false;                 ///< Use standard SVG &lt;text&gt; elements instead of Hershey paths.

  /** @brief Color palette selections */
  enum class Palette : std::uint8_t {
    OkabeIto,  ///< Colorblind-safe palette
    Grayscale, ///< B&W printing friendly
    Viridis,   ///< Perceptually uniform
    Magma,     ///< Sequential dark to bright purple/yellow
    Heatmap,   ///< Hot thermal gradient (Black-Red-Yellow-White)
    Rainbow,   ///< Standard Jet/HSV rainbow gradient
    Turbo,     ///< Google Turbo perceptually smooth rainbow
    Plasma,    ///< Sequential purple to yellow/orange
    Inferno,   ///< Sequential dark to bright yellow
    Cividis    ///< Colorblind-optimized perceptually uniform
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

/**
 * @enum TextAnchor
 * @brief Horizontal alignment / anchor options for SVG text rendering.
 */
enum class TextAnchor : std::uint8_t {
  Start,  ///< Left-aligned / start-anchored text.
  Middle, ///< Center-aligned / middle-anchored text.
  End     ///< Right-aligned / end-anchored text.
};

/**
 * @brief Custom rendering style for individual comparison curves.
 */
struct CurveStyle {
  std::string color_hex;     ///< Custom hex color string (e.g. "#E63946"). Empty uses default palette.
  float stroke_width = 2.0F; ///< Stroke line thickness in px.
  int dash_style = 0;        ///< Dash pattern style (0: Solid, 1: Dashed, 2: Dotted).
  bool visible = true;       ///< Visibility toggle flag.
};

/**
 * @brief A labeled histogram for comparison rendering.
 */
struct LabeledHistogram {
  std::string label;                            ///< Run / dataset label.
  const correlation::analysis::Histogram *hist; ///< Pointer to histogram data.
  CurveStyle style{};                           ///< Per-curve custom style settings.
  bool is_difference = false;                   ///< Flag indicating if this is an overlaid difference curve.
};

/**
 * @brief Samples a 1D dataset y(x) at point x using linear interpolation with boundary clamping to [x_min, x_max].
 */
inline real_t sampleHistogramClamped(const std::vector<real_t> &x_bins, const std::vector<real_t> &y_vals,
                                     real_t x_val) {
  if (x_bins.empty() || y_vals.empty()) {
    return static_cast<real_t>(0.0);
  }
  if (x_val <= x_bins.front()) {
    return y_vals.front();
  }
  if (x_val >= x_bins.back()) {
    return y_vals.back();
  }
  auto it_idx = std::lower_bound(x_bins.begin(), x_bins.end(), x_val);
  std::size_t idx = std::distance(x_bins.begin(), it_idx);
  if (idx == 0) {
    return y_vals[0];
  }
  real_t x_0 = x_bins[idx - 1];
  real_t x_1 = x_bins[idx];
  real_t y_0 = y_vals[idx - 1];
  real_t y_1 = y_vals[idx];
  if (std::abs(x_1 - x_0) < static_cast<real_t>(1e-12)) {
    return y_0;
  }
  return y_0 + (x_val - x_0) / (x_1 - x_0) * (y_1 - y_0);
}

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
        "  <text x=\"{:.1f}\" y=\"{:.1f}\" font-family=\"'Outfit', 'Plus Jakarta Sans', 'Inter', 'Roboto', 'Helvetica Neue', "
        "sans-serif\" font-size=\"{:.1f}\" text-anchor=\"{}\" fill=\"{}\">{}</text>\n",
        x_pos, y_pos, size, anchor_str, color, text);
  }
  std::string path_d = Roboto::instance().render(TextRenderParameters{
      .text = text,
      .start_x = x_pos,
      .start_y = y_pos,
      .font_size = size,
      .anchor = anchor_str,
  });
  return std::format("  <path d=\"{}\" fill=\"{}\" fill-rule=\"evenodd\" stroke=\"none\"/>\n", path_d, color);
}

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

struct ColorStop {
  float t;
  uint8_t r, g, b;
};

inline std::string interpolateColorStops(float position, std::span<const ColorStop> stops) {
  if (stops.empty()) {
    return "#000000";
  }
  if (position <= stops.front().t) {
    return std::format("#{:02X}{:02X}{:02X}", static_cast<unsigned int>(stops.front().r),
                       static_cast<unsigned int>(stops.front().g), static_cast<unsigned int>(stops.front().b));
  }
  if (position >= stops.back().t) {
    return std::format("#{:02X}{:02X}{:02X}", static_cast<unsigned int>(stops.back().r),
                       static_cast<unsigned int>(stops.back().g), static_cast<unsigned int>(stops.back().b));
  }
  for (std::size_t i = 0; i < stops.size() - 1; ++i) {
    if (position >= stops[i].t && position <= stops[i + 1].t) {
      float factor = (position - stops[i].t) / (stops[i + 1].t - stops[i].t);
      auto red = static_cast<unsigned int>(
          std::round(std::lerp(static_cast<float>(stops[i].r), static_cast<float>(stops[i + 1].r), factor)));
      auto green = static_cast<unsigned int>(
          std::round(std::lerp(static_cast<float>(stops[i].g), static_cast<float>(stops[i + 1].g), factor)));
      auto blue = static_cast<unsigned int>(
          std::round(std::lerp(static_cast<float>(stops[i].b), static_cast<float>(stops[i + 1].b), factor)));
      return std::format("#{:02X}{:02X}{:02X}", red, green, blue);
    }
  }
  return std::format("#{:02X}{:02X}{:02X}", static_cast<unsigned int>(stops.back().r),
                     static_cast<unsigned int>(stops.back().g), static_cast<unsigned int>(stops.back().b));
}

inline std::string sampleContinuousColormap(float position, PlotConfig::Palette pal) {
  switch (pal) {
  case PlotConfig::Palette::Magma: {
    constexpr std::array<ColorStop, 5> kMagmaStops = {{{.t = 0.00F, .r = 0, .g = 0, .b = 4},
                                                       {.t = 0.25F, .r = 81, .g = 18, .b = 124},
                                                       {.t = 0.50F, .r = 183, .g = 55, .b = 121},
                                                       {.t = 0.75F, .r = 252, .g = 137, .b = 97},
                                                       {.t = 1.00F, .r = 252, .g = 253, .b = 191}}};
    return interpolateColorStops(position, kMagmaStops);
  }
  case PlotConfig::Palette::Heatmap: {
    constexpr std::array<ColorStop, 4> kHeatmapStops = {{{.t = 0.00F, .r = 0, .g = 0, .b = 0},
                                                         {.t = 0.33F, .r = 255, .g = 0, .b = 0},
                                                         {.t = 0.66F, .r = 255, .g = 255, .b = 0},
                                                         {.t = 1.00F, .r = 255, .g = 255, .b = 255}}};
    return interpolateColorStops(position, kHeatmapStops);
  }
  case PlotConfig::Palette::Rainbow: {
    constexpr std::array<ColorStop, 5> kRainbowStops = {{{.t = 0.00F, .r = 0, .g = 0, .b = 255},
                                                         {.t = 0.25F, .r = 0, .g = 255, .b = 255},
                                                         {.t = 0.50F, .r = 0, .g = 255, .b = 0},
                                                         {.t = 0.75F, .r = 255, .g = 255, .b = 0},
                                                         {.t = 1.00F, .r = 255, .g = 0, .b = 0}}};
    return interpolateColorStops(position, kRainbowStops);
  }
  case PlotConfig::Palette::Turbo: {
    constexpr std::array<ColorStop, 5> kTurboStops = {{{.t = 0.00F, .r = 48, .g = 18, .b = 59},
                                                       {.t = 0.25F, .r = 26, .g = 228, .b = 182},
                                                       {.t = 0.50F, .r = 162, .g = 252, .b = 60},
                                                       {.t = 0.75F, .r = 251, .g = 128, .b = 34},
                                                       {.t = 1.00F, .r = 122, .g = 4, .b = 3}}};
    return interpolateColorStops(position, kTurboStops);
  }
  case PlotConfig::Palette::Plasma: {
    constexpr std::array<ColorStop, 5> kPlasmaStops = {{{.t = 0.00F, .r = 13, .g = 8, .b = 135},
                                                        {.t = 0.25F, .r = 106, .g = 0, .b = 168},
                                                        {.t = 0.50F, .r = 177, .g = 42, .b = 144},
                                                        {.t = 0.75F, .r = 225, .g = 100, .b = 98},
                                                        {.t = 1.00F, .r = 252, .g = 166, .b = 54}}};
    return interpolateColorStops(position, kPlasmaStops);
  }
  case PlotConfig::Palette::Inferno: {
    constexpr std::array<ColorStop, 5> kInfernoStops = {{{.t = 0.00F, .r = 0, .g = 0, .b = 4},
                                                         {.t = 0.25F, .r = 87, .g = 9, .b = 107},
                                                         {.t = 0.50F, .r = 187, .g = 55, .b = 84},
                                                         {.t = 0.75F, .r = 249, .g = 142, .b = 9},
                                                         {.t = 1.00F, .r = 252, .g = 255, .b = 164}}};
    return interpolateColorStops(position, kInfernoStops);
  }
  case PlotConfig::Palette::Cividis: {
    constexpr std::array<ColorStop, 5> kCividisStops = {{{.t = 0.00F, .r = 0, .g = 32, .b = 81},
                                                         {.t = 0.25F, .r = 58, .g = 71, .b = 108},
                                                         {.t = 0.50F, .r = 107, .g = 112, .b = 116},
                                                         {.t = 0.75F, .r = 161, .g = 156, .b = 114},
                                                         {.t = 1.00F, .r = 254, .g = 254, .b = 98}}};
    return interpolateColorStops(position, kCividisStops);
  }
  default:
    return "#000000";
  }
}

/**
 * @brief Retrieves a color from the selected palette.
 * @param index Index of curve.
 * @param total_count Total count of curves for rank interpolation.
 * @param pal Selected palette.
 * @return Hex color string.
 */
inline std::string color(std::size_t index, std::size_t total_count, PlotConfig::Palette pal) {
  switch (pal) {
  case PlotConfig::Palette::Grayscale:
    return std::string(kGrayscale.at(index % kGrayscale.size()));
  case PlotConfig::Palette::Viridis:
    return std::string(kViridis.at(index % kViridis.size()));
  case PlotConfig::Palette::OkabeIto:
    return std::string(kColors.at(index % kColors.size()));
  default: {
    float position = (total_count <= 1) ? 0.5F : static_cast<float>(index) / static_cast<float>(total_count - 1);
    return sampleContinuousColormap(position, pal);
  }
  }
}

/**
 * @brief Legacy overload retrieving color with default count.
 */
inline std::string color(std::size_t index, PlotConfig::Palette pal) { return color(index, 8, pal); }

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

/**
 * @struct DataRange
 * @brief Represents the minimum and maximum scalar bounds of a dataset.
 */
struct DataRange {
  real_t min = static_cast<real_t>(0.0); ///< Minimum value in data space.
  real_t max = static_cast<real_t>(0.0); ///< Maximum value in data space.
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
 */
inline std::string fmtScientific(real_t value) {
  real_t abs_value = std::abs(value);
  if (abs_value < static_cast<real_t>(1e-12)) {
    return "0";
  }

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

} // namespace detail
} // namespace correlation::plotters
