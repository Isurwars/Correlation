/**
 * @file PdfPlotter.cpp
 * @brief PDF generator implementation for distribution function histograms using pdfgen.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "plotters/PdfPlotter.hpp"

#include "analysis/DistributionFunctions.hpp"
#include "math/Precision.hpp"
#include "plotters/PlotTypes.hpp"
#include "plotters/SvgPlotter.hpp"
#include "plotters/pdfgen.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <map>
#include <numbers>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace correlation::plotters {

namespace {

// Constants
constexpr real_t PI = std::numbers::pi_v<real_t>;
constexpr uint32_t PDF_TRANSPARENT_COLOR = 0xff000000U;
constexpr uint32_t PDF_COLOR_BLACK = 0x00000000U;

[[nodiscard]] constexpr uint32_t pdfRgb(uint32_t red, uint32_t green, uint32_t blue) noexcept {
  return ((red & 0xffU) << 16U) | ((green & 0xffU) << 8U) | (blue & 0xffU);
}

// Extract color from string "#RRGGBB"
[[nodiscard]] uint32_t parseHexColor(std::string_view hex) {
  if (hex.length() == 7 && hex[0] == '#') {
    unsigned int red = 0;
    unsigned int green = 0;
    unsigned int blue = 0;
    std::from_chars(hex.data() + 1, hex.data() + 3, red, 16);
    std::from_chars(hex.data() + 3, hex.data() + 5, green, 16);
    std::from_chars(hex.data() + 5, hex.data() + 7, blue, 16);
    return pdfRgb(red, green, blue);
  }
  return PDF_COLOR_BLACK;
}

struct BlendParams {
  uint32_t fg;
  uint32_t bg;
};

// Blend foreground color with background color using a given opacity to simulate transparency
[[nodiscard]] uint32_t blendColor(BlendParams colors, real_t opacity) {
  auto const r_fg = static_cast<real_t>((colors.fg >> 16U) & 0xffU);
  auto const g_fg = static_cast<real_t>((colors.fg >> 8U) & 0xffU);
  auto const b_fg = static_cast<real_t>(colors.fg & 0xffU);

  auto const r_bg = static_cast<real_t>((colors.bg >> 16U) & 0xffU);
  auto const g_bg = static_cast<real_t>((colors.bg >> 8U) & 0xffU);
  auto const b_bg = static_cast<real_t>(colors.bg & 0xffU);

  auto const red =
      static_cast<uint32_t>(r_fg * opacity + r_bg * (static_cast<real_t>(1.0) - opacity));
  auto const green =
      static_cast<uint32_t>(g_fg * opacity + g_bg * (static_cast<real_t>(1.0) - opacity));
  auto const blue =
      static_cast<uint32_t>(b_fg * opacity + b_bg * (static_cast<real_t>(1.0) - opacity));

  return pdfRgb(red, green, blue);
}

// Format numbers in standard ASCII scientific notation for PDF compatibility
[[nodiscard]] std::string fmtScientificPdf(real_t value) {
  real_t const abs_value = std::abs(value);
  if (abs_value < static_cast<real_t>(1e-12)) {
    return "0";
  }

  if (abs_value < static_cast<real_t>(0.001) || abs_value >= static_cast<real_t>(10000.0)) {
    auto const exponent = static_cast<int>(std::floor(std::log10(abs_value)));
    real_t const fraction =
        value / std::pow(static_cast<real_t>(10.0), static_cast<real_t>(exponent));
    return std::format("{:.1f}x10^{}", fraction, exponent);
  }

  std::string str = std::format("{:.2f}", value);
  auto const dot = str.find('.');
  if (dot != std::string::npos) {
    std::size_t const last = str.find_last_not_of('0');
    if (last != std::string::npos && last > dot) {
      str = str.substr(0, last + 1);
    } else if (last == dot) {
      str = str.substr(0, dot);
    }
  }
  return str;
}

// Sanitize unit strings to use ASCII superscripts instead of Unicode
[[nodiscard]] std::string sanitizeUnitPdf(const std::string &unit) {
  std::string clean = unit;
  auto const replace_all = [](std::string &str, const std::string &search_for,
                              const std::string &replacement) {
    size_t start_pos = 0;
    while ((start_pos = str.find(search_for, start_pos)) != std::string::npos) {
      str.replace(start_pos, search_for.length(), replacement);
      start_pos += replacement.length();
    }
  };
  replace_all(clean, "⁻¹", "^-1");
  replace_all(clean, "⁻", "^-");
  replace_all(clean, "⁰", "^0");
  replace_all(clean, "¹", "^1");
  replace_all(clean, "²", "^2");
  replace_all(clean, "³", "^3");
  replace_all(clean, "⁴", "^4");
  replace_all(clean, "⁵", "^5");
  replace_all(clean, "⁶", "^6");
  replace_all(clean, "⁷", "^7");
  replace_all(clean, "⁸", "^8");
  replace_all(clean, "⁹", "^9");
  replace_all(clean, "×", "x");
  return clean;
}

struct PdfPoint {
  real_t x;
  real_t y;
};

// Helper to draw text with specified alignment (TextAnchor) and optional rotation angle
void drawPdfText(pdf_doc *pdf, struct pdf_object *page, const std::string &text, PdfPoint pos,
                 real_t size, TextAnchor anchor, uint32_t color,
                 real_t angle = static_cast<real_t>(0.0)) {
  std::string const sanitized = sanitizeUnitPdf(text);
  float width = 0.0F;
  pdf_get_font_text_width(pdf, "Helvetica", sanitized.c_str(), static_cast<float>(size), &width);

  auto offset_factor = static_cast<real_t>(0.0);
  if (anchor == TextAnchor::Middle) {
    offset_factor = static_cast<real_t>(0.5);
  } else if (anchor == TextAnchor::End) {
    offset_factor = static_cast<real_t>(1.0);
  }

  real_t const offset_x = -offset_factor * static_cast<real_t>(width) * std::cos(angle);
  real_t const offset_y = -offset_factor * static_cast<real_t>(width) * std::sin(angle);
  real_t const x_aligned = pos.x + offset_x;
  real_t const y_aligned = pos.y + offset_y;

  if (angle != static_cast<real_t>(0.0)) {
    pdf_add_text_rotate(pdf, page, sanitized.c_str(), static_cast<float>(size),
                        static_cast<float>(x_aligned), static_cast<float>(y_aligned),
                        static_cast<float>(angle), color);
  } else {
    pdf_add_text(pdf, page, sanitized.c_str(), static_cast<float>(size),
                 static_cast<float>(x_aligned), static_cast<float>(y_aligned), color);
  }
}

struct PdfHistogramRenderer {
  const correlation::analysis::Histogram *hist;
  const PlotConfig *config;
  pdf_doc *pdf;
  struct pdf_object *page;

  real_t canvas_width = static_cast<real_t>(0.0);
  real_t canvas_height = static_cast<real_t>(0.0);
  real_t px0 = static_cast<real_t>(0.0);
  real_t px1 = static_cast<real_t>(0.0);
  real_t py0 = static_cast<real_t>(0.0);
  real_t py1 = static_cast<real_t>(0.0);

  uint32_t bg_col = 0;
  uint32_t axis_col = 0;
  uint32_t grid_col = 0;
  uint32_t text_col = 0;

  real_t min_x = static_cast<real_t>(0.0);
  real_t max_x = static_cast<real_t>(0.0);
  real_t min_y = static_cast<real_t>(0.0);
  real_t max_y = static_cast<real_t>(0.0);

  detail::NiceScale xScale;
  detail::NiceScale yScale;
  std::map<std::string, std::vector<real_t>> partials;

  PdfHistogramRenderer(const correlation::analysis::Histogram &histogram, const PlotConfig &cfg,
                       pdf_doc *pdf_doc_ptr, struct pdf_object *pdf_page)
      : hist(&histogram), config(&cfg), pdf(pdf_doc_ptr), page(pdf_page),
        canvas_width(cfg.effectiveWidth()), canvas_height(cfg.effectiveHeight()),
        px0(static_cast<real_t>(100.0)), px1(cfg.effectiveWidth() - static_cast<real_t>(40.0)),
        py0(static_cast<real_t>(50.0)), py1(cfg.effectiveHeight() - static_cast<real_t>(90.0)),
        bg_col(parseHexColor(cfg.bgColor())), axis_col(parseHexColor(cfg.axisColor())),
        grid_col(parseHexColor(cfg.gridColor())), text_col(parseHexColor(cfg.textColor())) {}

  [[nodiscard]] real_t toPdfY(real_t y_svg) const { return canvas_height - y_svg; }

  bool prepareData() {
    min_x = hist->bins.front();
    max_x = hist->bins.back();
    min_y = static_cast<real_t>(0.0);
    max_y = static_cast<real_t>(0.0);

    partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;

    for (const auto &[key, vals] : partials) {
      for (real_t const val : vals) {
        max_y = std::max(max_y, val);
        min_y = std::min(min_y, val);
      }
    }

    if (max_y == min_y) {
      max_y += static_cast<real_t>(1.0);
    }
    if (max_x == min_x) {
      max_x += static_cast<real_t>(1.0);
    }

    real_t const y_padding = (max_y - min_y) * static_cast<real_t>(0.05);
    max_y += y_padding;

    xScale = detail::NiceScale(
        detail::DataRange{
            .min = min_x,
            .max = max_x,
        },
        11);
    yScale = detail::NiceScale(
        detail::DataRange{
            .min = min_y,
            .max = max_y,
        },
        8);
    return true;
  }

  void drawGridAndTicks() const {
    for (real_t const y_val : yScale.ticks) {
      real_t const spy = detail::mapValue(y_val, yScale.min, yScale.max, py1, py0);
      real_t const pdf_y = toPdfY(spy);
      if (config->show_grid) {
        pdf_add_line(pdf, page, static_cast<float>(px0), static_cast<float>(pdf_y),
                     static_cast<float>(px1), static_cast<float>(pdf_y), 0.8F, grid_col);
      }
      pdf_add_line(pdf, page, static_cast<float>(px0 - static_cast<real_t>(8.0)),
                   static_cast<float>(pdf_y), static_cast<float>(px0), static_cast<float>(pdf_y),
                   1.5F, axis_col);
      drawPdfText(pdf, page, fmtScientificPdf(y_val),
                  {
                      .x = px0 - static_cast<real_t>(15.0),
                      .y = pdf_y - static_cast<real_t>(5.0) * config->font_scale,
                  },
                  static_cast<real_t>(20.0) * config->font_scale, TextAnchor::End, text_col);
    }

    for (real_t const x_val : xScale.ticks) {
      real_t const spx = detail::mapValue(x_val, xScale.min, xScale.max, px0, px1);
      real_t const pdf_y = toPdfY(py1);
      if (config->show_grid) {
        pdf_add_line(pdf, page, static_cast<float>(spx), static_cast<float>(toPdfY(py0)),
                     static_cast<float>(spx), static_cast<float>(pdf_y), 0.8F, grid_col);
      }
      pdf_add_line(pdf, page, static_cast<float>(spx),
                   static_cast<float>(pdf_y - static_cast<real_t>(8.0)), static_cast<float>(spx),
                   static_cast<float>(pdf_y), 1.5F, axis_col);
      drawPdfText(pdf, page, fmtScientificPdf(x_val),
                  {.x = spx, .y = pdf_y - static_cast<real_t>(25.0) * config->font_scale},
                  static_cast<real_t>(20.0) * config->font_scale, TextAnchor::Middle, text_col);
    }

    pdf_add_rectangle(pdf, page, static_cast<float>(px0), static_cast<float>(toPdfY(py1)),
                      static_cast<float>(px1 - px0), static_cast<float>(py1 - py0), 1.5F, axis_col);
  }

  void drawEmphasisLines() const {
    std::vector<real_t> emphasis_values;
    std::string const title_str = hist->title;
    if (title_str.contains("g(r)")) {
      emphasis_values.push_back(static_cast<real_t>(1.0));
    } else if (title_str.contains("G(r)")) {
      emphasis_values.push_back(static_cast<real_t>(0.0));
    } else if (title_str.contains("S(Q)") || title_str.contains("S(q)")) {
      emphasis_values.push_back(static_cast<real_t>(0.0));
      emphasis_values.push_back(static_cast<real_t>(1.0));
    }

    for (real_t const focus_y : emphasis_values) {
      if (focus_y >= yScale.min && focus_y <= yScale.max) {
        real_t const spy = detail::mapValue(focus_y, yScale.min, yScale.max, py1, py0);
        real_t const pdf_y = toPdfY(spy);
        pdf_add_line(pdf, page, static_cast<float>(px0), static_cast<float>(pdf_y),
                     static_cast<float>(px1), static_cast<float>(pdf_y), 1.8F, axis_col);
      }
    }
  }

  void drawAreaFills() const {
    if (config->fill_area) {
      std::size_t fill_ci = 0;
      for (const auto &[key, vals] : partials) {
        uint32_t const col = parseHexColor(detail::color(fill_ci++, config->palette));
        uint32_t const shade_col = blendColor(
            {
                .fg = col,
                .bg = bg_col,
            },
            static_cast<real_t>(0.20));
        std::size_t const num_points = std::min(hist->bins.size(), vals.size());
        if (num_points > 1) {
          std::vector<float> x_coords;
          std::vector<float> y_coords;
          x_coords.reserve(num_points + 2);
          y_coords.reserve(num_points + 2);

          real_t const sx_start =
              detail::mapValue(hist->bins.front(), xScale.min, xScale.max, px0, px1);
          x_coords.push_back(static_cast<float>(sx_start));
          y_coords.push_back(static_cast<float>(toPdfY(py1)));

          for (std::size_t i = 0; i < num_points; ++i) {
            real_t const screen_x =
                detail::mapValue(hist->bins[i], xScale.min, xScale.max, px0, px1);
            real_t const screen_y = detail::mapValue(vals[i], yScale.min, yScale.max, py1, py0);
            real_t const sy_clamped = std::min(screen_y, py1);

            x_coords.push_back(static_cast<float>(screen_x));
            y_coords.push_back(static_cast<float>(toPdfY(sy_clamped)));
          }

          real_t const sx_end =
              detail::mapValue(hist->bins[num_points - 1], xScale.min, xScale.max, px0, px1);
          x_coords.push_back(static_cast<float>(sx_end));
          y_coords.push_back(static_cast<float>(toPdfY(py1)));

          pdf_add_filled_polygon(pdf, page, x_coords.data(), y_coords.data(),
                                 static_cast<int>(x_coords.size()), 0.0F, shade_col);
        }
      }
    }
  }

  [[nodiscard]] std::vector<std::pair<std::string, uint32_t>> drawPolylines() const {
    std::vector<std::pair<std::string, uint32_t>> legend_items;
    std::size_t color_idx = 0;
    for (const auto &[key, vals] : partials) {
      uint32_t const col = parseHexColor(detail::color(color_idx++, config->palette));
      legend_items.emplace_back(key, col);

      std::size_t const num_points = std::min(hist->bins.size(), vals.size());
      for (size_t i = 1; i < num_points; ++i) {
        real_t const x_1 = detail::mapValue(hist->bins[i - 1], xScale.min, xScale.max, px0, px1);
        real_t const y_1 = detail::mapValue(vals[i - 1], yScale.min, yScale.max, py1, py0);
        real_t const x_2 = detail::mapValue(hist->bins[i], xScale.min, xScale.max, px0, px1);
        real_t const y_2 = detail::mapValue(vals[i], yScale.min, yScale.max, py1, py0);

        pdf_add_line(pdf, page, static_cast<float>(x_1), static_cast<float>(toPdfY(y_1)),
                     static_cast<float>(x_2), static_cast<float>(toPdfY(y_2)),
                     static_cast<float>(config->line_width), col);
      }
    }
    return legend_items;
  }

  void drawMarkers() const {
    if (config->show_markers) {
      std::size_t marker_ci = 0;
      for (const auto &[key, vals] : partials) {
        uint32_t const col = parseHexColor(detail::color(marker_ci++, config->palette));
        std::size_t const num_points = std::min(hist->bins.size(), vals.size());
        for (std::size_t i = 0; i < num_points; ++i) {
          real_t const screen_x = detail::mapValue(hist->bins[i], xScale.min, xScale.max, px0, px1);
          real_t const screen_y = detail::mapValue(vals[i], yScale.min, yScale.max, py1, py0);
          pdf_add_circle(pdf, page, static_cast<float>(screen_x),
                         static_cast<float>(toPdfY(screen_y)),
                         static_cast<float>(config->marker_size), 0.0F, col, col);
        }
      }
    }
  }

  void drawLegend(const std::vector<std::pair<std::string, uint32_t>> &legend_items) const {
    if (config->show_legend) {
      real_t const legend_x = px1 - static_cast<real_t>(15.0);
      real_t legend_y = py0 + static_cast<real_t>(25.0);
      for (const auto &[label, col] : std::views::reverse(legend_items)) {
        real_t const pdf_y = toPdfY(legend_y);
        pdf_add_line(pdf, page, static_cast<float>(legend_x - static_cast<real_t>(40.0)),
                     static_cast<float>(pdf_y),
                     static_cast<float>(legend_x - static_cast<real_t>(10.0)),
                     static_cast<float>(pdf_y), 4.0F, col);
        drawPdfText(pdf, page, label,
                    {
                        .x = legend_x - static_cast<real_t>(45.0),
                        .y = pdf_y - static_cast<real_t>(5.0) * config->font_scale,
                    },
                    static_cast<real_t>(18.0) * config->font_scale, TextAnchor::End, text_col);
        legend_y += static_cast<real_t>(28.0);
      }
    }
  }

  void drawTitlesAndLabels() const {
    std::string const title_txt = hist->title.empty() ? "Histogram" : hist->title;
    std::string x_label = hist->x_label.empty() ? "x" : hist->x_label;
    std::string y_label = hist->y_label.empty() ? "y" : hist->y_label;

    if (!hist->x_unit.empty()) {
      x_label += " (" + hist->x_unit + ")";
    }
    if (!hist->y_unit.empty()) {
      std::string y_unit = hist->y_unit;
      if (y_unit == "Å^-1") {
        y_unit = "Å⁻¹";
      }
      y_label += " (" + y_unit + ")";
    }

    drawPdfText(
        pdf, page, title_txt,
        {.x = canvas_width / static_cast<real_t>(2.0), .y = toPdfY(static_cast<real_t>(30.0))},
        static_cast<real_t>(24.0) * config->font_scale, TextAnchor::Middle, text_col);
    drawPdfText(
        pdf, page, x_label,
        {.x = (px0 + px1) / static_cast<real_t>(2.0), .y = toPdfY(py1 + static_cast<real_t>(75.0))},
        static_cast<real_t>(28.0) * config->font_scale, TextAnchor::Middle, text_col);
    drawPdfText(
        pdf, page, y_label,
        {.x = static_cast<real_t>(40.0), .y = toPdfY((py0 + py1) / static_cast<real_t>(2.0))},
        static_cast<real_t>(28.0) * config->font_scale, TextAnchor::Middle, text_col,
        PI / static_cast<real_t>(2.0));
  }
};

struct PdfComparisonRenderer {
  const std::vector<LabeledHistogram> *datasets;
  const std::string *query_key;
  const PlotConfig *config;
  pdf_doc *pdf;
  struct pdf_object *page;

  real_t canvas_width = static_cast<real_t>(0.0);
  real_t canvas_height = static_cast<real_t>(0.0);
  real_t px0 = static_cast<real_t>(0.0);
  real_t px1 = static_cast<real_t>(0.0);
  real_t py0 = static_cast<real_t>(0.0);
  real_t py1 = static_cast<real_t>(0.0);

  uint32_t bg_col = 0;
  uint32_t axis_col = 0;
  uint32_t grid_col = 0;
  uint32_t text_col = 0;

  real_t raw_x_min = static_cast<real_t>(1e30);
  real_t raw_x_max = static_cast<real_t>(-1e30);
  real_t raw_y_min = static_cast<real_t>(0.0);
  real_t raw_y_max = static_cast<real_t>(-1e30);

  detail::NiceScale xScale;
  detail::NiceScale yScale;
  std::vector<std::pair<std::string, uint32_t>> legend_items;

  PdfComparisonRenderer(const std::vector<LabeledHistogram> &datasets_ref, const std::string &key,
                        const PlotConfig &cfg, pdf_doc *pdf_doc_ptr, struct pdf_object *pdf_page)
      : datasets(&datasets_ref), query_key(&key), config(&cfg), pdf(pdf_doc_ptr), page(pdf_page),
        canvas_width(cfg.effectiveWidth()), canvas_height(cfg.effectiveHeight()),
        px0(static_cast<real_t>(100.0)), px1(cfg.effectiveWidth() - static_cast<real_t>(40.0)),
        py0(static_cast<real_t>(50.0)), py1(cfg.effectiveHeight() - static_cast<real_t>(90.0)),
        bg_col(parseHexColor(cfg.bgColor())), axis_col(parseHexColor(cfg.axisColor())),
        grid_col(parseHexColor(cfg.gridColor())), text_col(parseHexColor(cfg.textColor())) {}

  [[nodiscard]] real_t toPdfY(real_t y_svg) const { return canvas_height - y_svg; }

  bool prepareData() {
    for (const auto &dataset : *datasets) {
      const auto &partials = dataset.hist->smoothed_partials.empty()
                                 ? dataset.hist->partials
                                 : dataset.hist->smoothed_partials;
      auto const iterator = partials.find(*query_key);
      if (iterator == partials.end()) {
        continue;
      }
      const auto &x_values = dataset.hist->bins;
      const auto &y_values = iterator->second;
      if (!x_values.empty()) {
        raw_x_min = std::min(raw_x_min, x_values.front());
        raw_x_max = std::max(raw_x_max, x_values.back());
      }
      for (real_t const y_value : y_values) {
        raw_y_max = std::max(raw_y_max, y_value);
        raw_y_min = std::min(raw_y_min, y_value);
      }
    }

    if (raw_x_max <= raw_x_min || raw_y_max <= raw_y_min) {
      return false;
    }

    real_t const y_padding = (raw_y_max - raw_y_min) * static_cast<real_t>(0.05);
    raw_y_max += y_padding;

    xScale = detail::NiceScale(
        detail::DataRange{
            .min = raw_x_min,
            .max = raw_x_max,
        },
        11);
    yScale = detail::NiceScale(
        detail::DataRange{
            .min = raw_y_min,
            .max = raw_y_max,
        },
        8);
    return true;
  }

  void drawGridAndTicks() const {
    for (real_t const y_val : yScale.ticks) {
      real_t const spy = detail::mapValue(y_val, yScale.min, yScale.max, py1, py0);
      real_t const pdf_y = toPdfY(spy);
      if (config->show_grid) {
        pdf_add_line(pdf, page, static_cast<float>(px0), static_cast<float>(pdf_y),
                     static_cast<float>(px1), static_cast<float>(pdf_y), 0.8F, grid_col);
      }
      pdf_add_line(pdf, page, static_cast<float>(px0 - static_cast<real_t>(8.0)),
                   static_cast<float>(pdf_y), static_cast<float>(px0), static_cast<float>(pdf_y),
                   1.5F, axis_col);
      drawPdfText(pdf, page, fmtScientificPdf(y_val),
                  {
                      .x = px0 - static_cast<real_t>(15.0),
                      .y = pdf_y - static_cast<real_t>(5.0) * config->font_scale,
                  },
                  static_cast<real_t>(20.0) * config->font_scale, TextAnchor::End, text_col);
    }

    for (real_t const x_val : xScale.ticks) {
      real_t const spx = detail::mapValue(x_val, xScale.min, xScale.max, px0, px1);
      real_t const pdf_y = toPdfY(py1);
      if (config->show_grid) {
        pdf_add_line(pdf, page, static_cast<float>(spx), static_cast<float>(toPdfY(py0)),
                     static_cast<float>(spx), static_cast<float>(pdf_y), 0.8F, grid_col);
      }
      pdf_add_line(pdf, page, static_cast<float>(spx),
                   static_cast<float>(pdf_y - static_cast<real_t>(8.0)), static_cast<float>(spx),
                   static_cast<float>(pdf_y), 1.5F, axis_col);
      drawPdfText(pdf, page, fmtScientificPdf(x_val),
                  {.x = spx, .y = pdf_y - static_cast<real_t>(25.0) * config->font_scale},
                  static_cast<real_t>(20.0) * config->font_scale, TextAnchor::Middle, text_col);
    }

    pdf_add_rectangle(pdf, page, static_cast<float>(px0), static_cast<float>(toPdfY(py1)),
                      static_cast<float>(px1 - px0), static_cast<float>(py1 - py0), 1.5F, axis_col);
  }

  void drawAreaFills() const {
    if (config->fill_area) {
      std::size_t fill_ci = 0;
      for (const auto &dataset : *datasets) {
        const auto &partials = dataset.hist->smoothed_partials.empty()
                                   ? dataset.hist->partials
                                   : dataset.hist->smoothed_partials;
        auto const iterator = partials.find(*query_key);
        if (iterator != partials.end()) {
          const auto &x_values = dataset.hist->bins;
          const auto &y_values = iterator->second;
          uint32_t const col = parseHexColor(detail::color(fill_ci, config->palette));
          uint32_t const shade_col = blendColor(
              {
                  .fg = col,
                  .bg = bg_col,
              },
              static_cast<real_t>(0.20));
          std::size_t const n_points = std::min(x_values.size(), y_values.size());
          if (n_points > 1) {
            std::vector<float> x_coords;
            std::vector<float> y_coords;
            x_coords.reserve(n_points + 2);
            y_coords.reserve(n_points + 2);

            real_t const sx_start = detail::mapValue(x_values[0], xScale.min, xScale.max, px0, px1);
            x_coords.push_back(static_cast<float>(sx_start));
            y_coords.push_back(static_cast<float>(toPdfY(py1)));

            for (std::size_t i = 0; i < n_points; ++i) {
              real_t const sp_x = detail::mapValue(x_values[i], xScale.min, xScale.max, px0, px1);
              real_t const sp_y = detail::mapValue(y_values[i], yScale.min, yScale.max, py1, py0);
              real_t const sp_y_clamped = std::min(sp_y, py1);

              x_coords.push_back(static_cast<float>(sp_x));
              y_coords.push_back(static_cast<float>(toPdfY(sp_y_clamped)));
            }

            real_t const sx_end =
                detail::mapValue(x_values[n_points - 1], xScale.min, xScale.max, px0, px1);
            x_coords.push_back(static_cast<float>(sx_end));
            y_coords.push_back(static_cast<float>(toPdfY(py1)));

            pdf_add_filled_polygon(pdf, page, x_coords.data(), y_coords.data(),
                                   static_cast<int>(x_coords.size()), 0.0F, shade_col);
          }
        }
        fill_ci++;
      }
    }
  }

  void drawPolylines() {
    std::size_t color_idx = 0;
    for (const auto &dataset : *datasets) {
      const auto &partials = dataset.hist->smoothed_partials.empty()
                                 ? dataset.hist->partials
                                 : dataset.hist->smoothed_partials;
      auto const partial_iter = partials.find(*query_key);
      if (partial_iter == partials.end()) {
        continue;
      }

      const auto &x_values = dataset.hist->bins;
      const auto &y_values = partial_iter->second;
      uint32_t const col = parseHexColor(detail::color(color_idx++, config->palette));
      legend_items.emplace_back(dataset.label, col);

      std::size_t const n_points = std::min(x_values.size(), y_values.size());
      for (size_t i = 1; i < n_points; ++i) {
        real_t const x_1 = detail::mapValue(x_values[i - 1], xScale.min, xScale.max, px0, px1);
        real_t const y_1 = detail::mapValue(y_values[i - 1], yScale.min, yScale.max, py1, py0);
        real_t const x_2 = detail::mapValue(x_values[i], xScale.min, xScale.max, px0, px1);
        real_t const y_2 = detail::mapValue(y_values[i], yScale.min, yScale.max, py1, py0);

        pdf_add_line(pdf, page, static_cast<float>(x_1), static_cast<float>(toPdfY(y_1)),
                     static_cast<float>(x_2), static_cast<float>(toPdfY(y_2)),
                     static_cast<float>(config->line_width), col);
      }
    }
  }

  void drawMarkers() const {
    if (config->show_markers) {
      std::size_t marker_color_idx = 0;
      for (const auto &dataset : *datasets) {
        const auto &partials = dataset.hist->smoothed_partials.empty()
                                   ? dataset.hist->partials
                                   : dataset.hist->smoothed_partials;
        auto const partial_iter = partials.find(*query_key);
        if (partial_iter != partials.end()) {
          const auto &x_values = dataset.hist->bins;
          const auto &y_values = partial_iter->second;
          uint32_t const col = parseHexColor(detail::color(marker_color_idx++, config->palette));
          std::size_t const n_points = std::min(x_values.size(), y_values.size());
          for (std::size_t i = 0; i < n_points; ++i) {
            real_t const sp_x = detail::mapValue(x_values[i], xScale.min, xScale.max, px0, px1);
            real_t const sp_y = detail::mapValue(y_values[i], yScale.min, yScale.max, py1, py0);
            pdf_add_circle(pdf, page, static_cast<float>(sp_x), static_cast<float>(toPdfY(sp_y)),
                           static_cast<float>(config->marker_size), 0.0F, col, col);
          }
        } else {
          marker_color_idx++;
        }
      }
    }
  }

  void drawLegend() const {
    if (config->show_legend) {
      real_t const legend_x = px1 - static_cast<real_t>(15.0);
      real_t legend_y = py0 + static_cast<real_t>(25.0);
      for (const auto &[label, col] : std::views::reverse(legend_items)) {
        real_t const pdf_y = toPdfY(legend_y);
        pdf_add_line(pdf, page, static_cast<float>(legend_x - static_cast<real_t>(40.0)),
                     static_cast<float>(pdf_y),
                     static_cast<float>(legend_x - static_cast<real_t>(10.0)),
                     static_cast<float>(pdf_y), 4.0F, col);
        drawPdfText(pdf, page, label,
                    {
                        .x = legend_x - static_cast<real_t>(45.0),
                        .y = pdf_y - static_cast<real_t>(5.0) * config->font_scale,
                    },
                    static_cast<real_t>(18.0) * config->font_scale, TextAnchor::End, text_col);
        legend_y += static_cast<real_t>(28.0);
      }
    }
  }

  void drawTitlesAndLabels() const {
    const auto &ref = *datasets->front().hist;
    std::string const title_txt = *query_key + " Comparison";
    std::string x_label = ref.x_label.empty() ? "x" : ref.x_label;
    std::string y_label = ref.y_label.empty() ? "y" : ref.y_label;

    if (!ref.x_unit.empty()) {
      x_label += " (" + ref.x_unit + ")";
    }
    if (!ref.y_unit.empty()) {
      std::string y_unit = ref.y_unit;
      if (y_unit == "Å^-1") {
        y_unit = "Å⁻¹";
      }
      y_label += " (" + y_unit + ")";
    }

    drawPdfText(
        pdf, page, title_txt,
        {.x = canvas_width / static_cast<real_t>(2.0), .y = toPdfY(static_cast<real_t>(30.0))},
        static_cast<real_t>(24.0) * config->font_scale, TextAnchor::Middle, text_col);
    drawPdfText(
        pdf, page, x_label,
        {.x = (px0 + px1) / static_cast<real_t>(2.0), .y = toPdfY(py1 + static_cast<real_t>(75.0))},
        static_cast<real_t>(28.0) * config->font_scale, TextAnchor::Middle, text_col);
    drawPdfText(
        pdf, page, y_label,
        {.x = static_cast<real_t>(40.0), .y = toPdfY((py0 + py1) / static_cast<real_t>(2.0))},
        static_cast<real_t>(28.0) * config->font_scale, TextAnchor::Middle, text_col,
        PI / static_cast<real_t>(2.0));
  }
};

} // namespace

void renderHistogramAsPdf(const correlation::analysis::Histogram &hist, const std::string &filepath,
                          const PlotConfig &config) {
  real_t const canvas_width = config.effectiveWidth();
  real_t const canvas_height = config.effectiveHeight();

  // Setup pdfgen
  pdf_info const info = {};
  pdf_doc *pdf =
      pdf_create(static_cast<float>(canvas_width), static_cast<float>(canvas_height), &info);
  struct pdf_object *page = pdf_append_page(pdf);

  uint32_t const bg_col = parseHexColor(config.bgColor());
  uint32_t const text_col = parseHexColor(config.textColor());

  // Background
  pdf_add_filled_rectangle(pdf, page, 0.0F, 0.0F, static_cast<float>(canvas_width),
                           static_cast<float>(canvas_height), 0.0F, bg_col, PDF_TRANSPARENT_COLOR);

  pdf_set_font(pdf, "Helvetica");

  const auto to_pdf_y = [&](real_t y_svg) { return canvas_height - y_svg; };

  if (hist.bins.empty() || (hist.smoothed_partials.empty() && hist.partials.empty())) {
    drawPdfText(pdf, page, "No data available",
                {
                    .x = canvas_width / static_cast<real_t>(2.0),
                    .y = to_pdf_y(canvas_height / static_cast<real_t>(2.0)),
                },
                static_cast<real_t>(24.0) * config.font_scale, TextAnchor::Middle, text_col);
    pdf_save(pdf, filepath.c_str());
    pdf_destroy(pdf);
    return;
  }

  PdfHistogramRenderer renderer(hist, config, pdf, page);
  if (renderer.prepareData()) {
    renderer.drawGridAndTicks();
    renderer.drawEmphasisLines();
    renderer.drawAreaFills();
    auto const legend_items = renderer.drawPolylines();
    renderer.drawMarkers();
    renderer.drawLegend(legend_items);
    renderer.drawTitlesAndLabels();
  }

  pdf_save(pdf, filepath.c_str());
  pdf_destroy(pdf);
}

void renderComparisonPdf(const std::vector<LabeledHistogram> &datasets,
                         const ComparisonQuery &query, const PlotConfig &config) {
  if (datasets.empty()) {
    return;
  }

  real_t const canvas_width = config.effectiveWidth();
  real_t const canvas_height = config.effectiveHeight();

  pdf_info const info = {};
  pdf_doc *pdf =
      pdf_create(static_cast<float>(canvas_width), static_cast<float>(canvas_height), &info);
  struct pdf_object *page = pdf_append_page(pdf);

  uint32_t const bg_col = parseHexColor(config.bgColor());
  uint32_t const text_col = parseHexColor(config.textColor());

  // Background
  pdf_add_filled_rectangle(pdf, page, 0.0F, 0.0F, static_cast<float>(canvas_width),
                           static_cast<float>(canvas_height), 0.0F, bg_col, PDF_TRANSPARENT_COLOR);

  pdf_set_font(pdf, "Helvetica");

  const auto to_pdf_y = [&](real_t y_svg) { return canvas_height - y_svg; };

  PdfComparisonRenderer renderer(datasets, query.key, config, pdf, page);
  if (renderer.prepareData()) {
    renderer.drawGridAndTicks();
    renderer.drawAreaFills();
    renderer.drawPolylines();
    renderer.drawMarkers();
    renderer.drawLegend();
    renderer.drawTitlesAndLabels();
  } else {
    drawPdfText(pdf, page, "No comparison data",
                {
                    .x = canvas_width / static_cast<real_t>(2.0),
                    .y = to_pdf_y(canvas_height / static_cast<real_t>(2.0)),
                },
                static_cast<real_t>(24.0) * config.font_scale, TextAnchor::Middle, text_col);
  }

  pdf_save(pdf, query.filepath.c_str());
  pdf_destroy(pdf);
}

} // namespace correlation::plotters
