/**
 * @file PdfPlotter.hpp
 * @brief PDF generator for distribution function histograms using pdfgen.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "plotters/SvgPlotter.hpp" // For PlotConfig, NiceScale, LabeledHistogram, color
#include "plotters/pdfgen.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <map>
#include <string>
#include <vector>

namespace correlation::plotters {

namespace detail {

// Constants
constexpr double kPi = 3.14159265358979323846;

// Extract color from string "#RRGGBB"
inline uint32_t parseHexColor(const std::string &hex) {
  if (hex.length() == 7 && hex[0] == '#') {
    uint32_t red = std::stoul(hex.substr(1, 2), nullptr, 16);
    uint32_t green = std::stoul(hex.substr(3, 2), nullptr, 16);
    uint32_t blue = std::stoul(hex.substr(5, 2), nullptr, 16);
    return PDF_RGB(red, green, blue);
  }
  return PDF_BLACK;
}

struct BlendParams {
  uint32_t fg;
  uint32_t bg;
};

// Blend foreground color with background color using a given opacity to simulate transparency
inline uint32_t blendColor(BlendParams colors, double opacity) {
  double r_fg = ((colors.fg >> 16) & 0xff);
  double g_fg = ((colors.fg >> 8) & 0xff);
  double b_fg = (colors.fg & 0xff);

  double r_bg = ((colors.bg >> 16) & 0xff);
  double g_bg = ((colors.bg >> 8) & 0xff);
  double b_bg = (colors.bg & 0xff);

  uint32_t red = static_cast<uint32_t>(r_fg * opacity + r_bg * (1.0 - opacity));
  uint32_t green = static_cast<uint32_t>(g_fg * opacity + g_bg * (1.0 - opacity));
  uint32_t blue = static_cast<uint32_t>(b_fg * opacity + b_bg * (1.0 - opacity));

  return PDF_RGB(red, green, blue);
}

// Format numbers in standard ASCII scientific notation for PDF compatibility
inline std::string fmtScientificPdf(double value) {
  double abs_value = std::abs(value);
  if (abs_value < 1e-12) {
    return "0";
  }

  if (abs_value < 0.001 || abs_value >= 10000.0) {
    int exponent = static_cast<int>(std::floor(std::log10(abs_value)));
    double fraction = value / std::pow(10.0, exponent);
    return std::format("{:.1f}x10^{}", fraction, exponent);
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

// Sanitize unit strings to use ASCII superscripts instead of Unicode
inline std::string sanitizeUnitPdf(const std::string &unit) {
  std::string clean = unit;
  auto replaceAll = [](std::string &str, const std::string &search_for, const std::string &replacement) {
    size_t start_pos = 0;
    while ((start_pos = str.find(search_for, start_pos)) != std::string::npos) {
      str.replace(start_pos, search_for.length(), replacement);
      start_pos += replacement.length();
    }
  };
  replaceAll(clean, "⁻¹", "^-1");
  replaceAll(clean, "⁻", "^-");
  replaceAll(clean, "⁰", "^0");
  replaceAll(clean, "¹", "^1");
  replaceAll(clean, "²", "^2");
  replaceAll(clean, "³", "^3");
  replaceAll(clean, "⁴", "^4");
  replaceAll(clean, "⁵", "^5");
  replaceAll(clean, "⁶", "^6");
  replaceAll(clean, "⁷", "^7");
  replaceAll(clean, "⁸", "^8");
  replaceAll(clean, "⁹", "^9");
  replaceAll(clean, "×", "x");
  return clean;
}

struct PdfPoint {
  double x;
  double y;
};

// Helper to draw text with specified alignment (TextAnchor) and optional rotation angle
inline void drawPdfText(pdf_doc *pdf, struct pdf_object *page, const std::string &text, PdfPoint pos, double size,
                        TextAnchor anchor, uint32_t color, double angle = 0.0) {
  std::string sanitized = sanitizeUnitPdf(text);
  float width = 0.0F;
  pdf_get_font_text_width(pdf, "Helvetica", sanitized.c_str(), static_cast<float>(size), &width);

  double offset_factor = 0.0;
  if (anchor == TextAnchor::Middle) {
    offset_factor = 0.5;
  } else if (anchor == TextAnchor::End) {
    offset_factor = 1.0;
  }

  double offset_x = -offset_factor * width * std::cos(angle);
  double offset_y = -offset_factor * width * std::sin(angle);
  double x_aligned = pos.x + offset_x;
  double y_aligned = pos.y + offset_y;

  if (angle != 0.0) {
    pdf_add_text_rotate(pdf, page, sanitized.c_str(), static_cast<float>(size), static_cast<float>(x_aligned),
                        static_cast<float>(y_aligned), static_cast<float>(angle), color);
  } else {
    pdf_add_text(pdf, page, sanitized.c_str(), static_cast<float>(size), static_cast<float>(x_aligned),
                 static_cast<float>(y_aligned), color);
  }
}

struct PdfHistogramRenderer {
  const correlation::analysis::Histogram *hist;
  const PlotConfig *config;
  pdf_doc *pdf;
  struct pdf_object *page;

  double canvas_width = 0.0;
  double canvas_height = 0.0;
  double px0 = 0.0;
  double px1 = 0.0;
  double py0 = 0.0;
  double py1 = 0.0;

  uint32_t bg_col = 0;
  uint32_t axis_col = 0;
  uint32_t grid_col = 0;
  uint32_t text_col = 0;

  double min_x = 0.0;
  double max_x = 0.0;
  double min_y = 0.0;
  double max_y = 0.0;

  detail::NiceScale xScale;
  detail::NiceScale yScale;
  std::map<std::string, std::vector<real_t>> partials;

  PdfHistogramRenderer(const correlation::analysis::Histogram &histogram, const PlotConfig &cfg, pdf_doc *pdf_doc_ptr,
                       struct pdf_object *pdf_page)
      : hist(&histogram),
        config(&cfg),
        pdf(pdf_doc_ptr),
        page(pdf_page),
        canvas_width(cfg.effective_width()),
        canvas_height(cfg.effective_height()),
        px0(100.0),
        px1(cfg.effective_width() - 40.0),
        py0(50.0),
        py1(cfg.effective_height() - 90.0),
        bg_col(detail::parseHexColor(cfg.bg_color())),
        axis_col(detail::parseHexColor(cfg.axis_color())),
        grid_col(detail::parseHexColor(cfg.grid_color())),
        text_col(detail::parseHexColor(cfg.text_color())) {}

  double toPdfY(double y_svg) const {
    return canvas_height - y_svg;
  }

  bool prepareData() {
    min_x = hist->bins.front();
    max_x = hist->bins.back();
    min_y = 0.0;
    max_y = 0.0;

    partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;

    for (const auto &[key, vals] : partials) {
      for (real_t val : vals) {
        max_y = std::max(max_y, static_cast<double>(val));
        min_y = std::min(min_y, static_cast<double>(val));
      }
    }

    if (max_y == min_y) {
      max_y += 1.0;
    }
    if (max_x == min_x) {
      max_x += 1.0;
    }

    double y_padding = (max_y - min_y) * 0.05;
    max_y += y_padding;

    xScale = detail::NiceScale(detail::DataRange{.min = min_x, .max = max_x}, 11);
    yScale = detail::NiceScale(detail::DataRange{.min = min_y, .max = max_y}, 8);
    return true;
  }

  void drawGridAndTicks() const {
    for (double y_val : yScale.ticks) {
      double spy = detail::mapValue(y_val, yScale.min, yScale.max, py1, py0);
      double pdf_y = toPdfY(spy);
      if (config->show_grid) {
        pdf_add_line(pdf, page, static_cast<float>(px0), static_cast<float>(pdf_y), static_cast<float>(px1),
                     static_cast<float>(pdf_y), 0.8F, grid_col);
      }
      pdf_add_line(pdf, page, static_cast<float>(px0 - 8.0), static_cast<float>(pdf_y), static_cast<float>(px0),
                   static_cast<float>(pdf_y), 1.5F, axis_col);
      detail::drawPdfText(pdf, page, detail::fmtScientificPdf(y_val), {px0 - 15.0, pdf_y - 5.0 * config->font_scale},
                          20.0 * config->font_scale, TextAnchor::End, text_col);
    }

    for (double x_val : xScale.ticks) {
      double spx = detail::mapValue(x_val, xScale.min, xScale.max, px0, px1);
      double pdf_y = toPdfY(py1);
      if (config->show_grid) {
        pdf_add_line(pdf, page, static_cast<float>(spx), static_cast<float>(toPdfY(py0)), static_cast<float>(spx),
                     static_cast<float>(pdf_y), 0.8F, grid_col);
      }
      pdf_add_line(pdf, page, static_cast<float>(spx), static_cast<float>(pdf_y - 8.0), static_cast<float>(spx),
                   static_cast<float>(pdf_y), 1.5F, axis_col);
      detail::drawPdfText(pdf, page, detail::fmtScientificPdf(x_val), {spx, pdf_y - 25.0 * config->font_scale},
                          20.0 * config->font_scale, TextAnchor::Middle, text_col);
    }

    pdf_add_rectangle(pdf, page, static_cast<float>(px0), static_cast<float>(toPdfY(py1)), static_cast<float>(px1 - px0),
                      static_cast<float>(py1 - py0), 1.5F, axis_col);
  }

  void drawEmphasisLines() const {
    std::vector<double> emphasis_values;
    std::string title_str = hist->title;
    if (title_str.contains("g(r)")) {
      emphasis_values.push_back(1.0);
    } else if (title_str.contains("G(r)")) {
      emphasis_values.push_back(0.0);
    } else if (title_str.contains("S(Q)") || title_str.contains("S(q)")) {
      emphasis_values.push_back(0.0);
      emphasis_values.push_back(1.0);
    }

    for (double focus_y : emphasis_values) {
      if (focus_y >= yScale.min && focus_y <= yScale.max) {
        double spy = detail::mapValue(focus_y, yScale.min, yScale.max, py1, py0);
        double pdf_y = toPdfY(spy);
        pdf_add_line(pdf, page, static_cast<float>(px0), static_cast<float>(pdf_y), static_cast<float>(px1),
                     static_cast<float>(pdf_y), 1.8F, axis_col);
      }
    }
  }

  void drawAreaFills() const {
    if (config->fill_area) {
      std::size_t fill_ci = 0;
      for (const auto &[key, vals] : partials) {
        uint32_t col = detail::parseHexColor(detail::color(fill_ci++, config->palette));
        uint32_t shade_col = detail::blendColor({.fg = col, .bg = bg_col}, 0.20);
        std::size_t num_points = std::min(hist->bins.size(), vals.size());
        if (num_points > 1) {
          std::vector<float> x_coords;
          std::vector<float> y_coords;
          x_coords.reserve(num_points + 2);
          y_coords.reserve(num_points + 2);

          double sx_start = detail::mapValue(hist->bins.front(), xScale.min, xScale.max, px0, px1);
          x_coords.push_back(static_cast<float>(sx_start));
          y_coords.push_back(static_cast<float>(toPdfY(py1)));

          for (std::size_t i = 0; i < num_points; ++i) {
            double screen_x = detail::mapValue(hist->bins[i], xScale.min, xScale.max, px0, px1);
            double screen_y = detail::mapValue(vals[i], yScale.min, yScale.max, py1, py0);
            double sy_clamped = std::min(screen_y, py1);

            x_coords.push_back(static_cast<float>(screen_x));
            y_coords.push_back(static_cast<float>(toPdfY(sy_clamped)));
          }

          double sx_end = detail::mapValue(hist->bins[num_points - 1], xScale.min, xScale.max, px0, px1);
          x_coords.push_back(static_cast<float>(sx_end));
          y_coords.push_back(static_cast<float>(toPdfY(py1)));

          pdf_add_filled_polygon(pdf, page, x_coords.data(), y_coords.data(), static_cast<int>(x_coords.size()), 0.0F,
                                 shade_col);
        }
      }
    }
  }

  std::vector<std::pair<std::string, uint32_t>> drawPolylines() const {
    std::vector<std::pair<std::string, uint32_t>> legend_items;
    std::size_t color_idx = 0;
    for (const auto &[key, vals] : partials) {
      uint32_t col = detail::parseHexColor(detail::color(color_idx++, config->palette));
      legend_items.emplace_back(key, col);

      std::size_t num_points = std::min(hist->bins.size(), vals.size());
      for (size_t i = 1; i < num_points; ++i) {
        double x_1 = detail::mapValue(hist->bins[i - 1], xScale.min, xScale.max, px0, px1);
        double y_1 = detail::mapValue(vals[i - 1], yScale.min, yScale.max, py1, py0);
        double x_2 = detail::mapValue(hist->bins[i], xScale.min, xScale.max, px0, px1);
        double y_2 = detail::mapValue(vals[i], yScale.min, yScale.max, py1, py0);

        pdf_add_line(pdf, page, static_cast<float>(x_1), static_cast<float>(toPdfY(y_1)), static_cast<float>(x_2),
                     static_cast<float>(toPdfY(y_2)), static_cast<float>(config->line_width), col);
      }
    }
    return legend_items;
  }

  void drawMarkers() const {
    if (config->show_markers) {
      std::size_t marker_ci = 0;
      for (const auto &[key, vals] : partials) {
        uint32_t col = detail::parseHexColor(detail::color(marker_ci++, config->palette));
        std::size_t num_points = std::min(hist->bins.size(), vals.size());
        for (std::size_t i = 0; i < num_points; ++i) {
          double screen_x = detail::mapValue(hist->bins[i], xScale.min, xScale.max, px0, px1);
          double screen_y = detail::mapValue(vals[i], yScale.min, yScale.max, py1, py0);
          pdf_add_circle(pdf, page, static_cast<float>(screen_x), static_cast<float>(toPdfY(screen_y)), 3.5F, 0.0F, col,
                         col);
        }
      }
    }
  }

  void drawLegend(const std::vector<std::pair<std::string, uint32_t>> &legend_items) const {
    if (config->show_legend) {
      double legend_x = px1 - 15.0;
      double legend_y = py0 + 25.0;
      for (auto it = legend_items.rbegin(); it != legend_items.rend(); ++it) {
        const auto &[label, col] = *it;
        double pdf_y = toPdfY(legend_y);
        pdf_add_line(pdf, page, static_cast<float>(legend_x - 40.0), static_cast<float>(pdf_y),
                     static_cast<float>(legend_x - 10.0), static_cast<float>(pdf_y), 4.0F, col);
        detail::drawPdfText(pdf, page, label, {legend_x - 45.0, pdf_y - 5.0 * config->font_scale}, 18.0 * config->font_scale,
                            TextAnchor::End, text_col);
        legend_y += 28.0;
      }
    }
  }

  void drawTitlesAndLabels() const {
    std::string title_txt = hist->title.empty() ? "Histogram" : hist->title;
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

    detail::drawPdfText(pdf, page, title_txt, {canvas_width / 2.0, toPdfY(30.0)}, 24.0 * config->font_scale, TextAnchor::Middle,
                        text_col);
    detail::drawPdfText(pdf, page, x_label, {(px0 + px1) / 2.0, toPdfY(py1 + 75.0)}, 28.0 * config->font_scale,
                        TextAnchor::Middle, text_col);
    detail::drawPdfText(pdf, page, y_label, {40.0, toPdfY((py0 + py1) / 2.0)}, 28.0 * config->font_scale, TextAnchor::Middle,
                        text_col, detail::kPi / 2.0);
  }
};

struct PdfComparisonRenderer {
  const std::vector<LabeledHistogram> *datasets;
  const std::string *query_key;
  const PlotConfig *config;
  pdf_doc *pdf;
  struct pdf_object *page;

  double canvas_width = 0.0;
  double canvas_height = 0.0;
  double px0 = 0.0;
  double px1 = 0.0;
  double py0 = 0.0;
  double py1 = 0.0;

  uint32_t bg_col = 0;
  uint32_t axis_col = 0;
  uint32_t grid_col = 0;
  uint32_t text_col = 0;

  double raw_x_min = 1e30;
  double raw_x_max = -1e30;
  double raw_y_min = 0.0;
  double raw_y_max = -1e30;

  detail::NiceScale xScale;
  detail::NiceScale yScale;
  std::vector<std::pair<std::string, uint32_t>> legend_items;

  PdfComparisonRenderer(const std::vector<LabeledHistogram> &datasets_ref, const std::string &key, const PlotConfig &cfg,
                        pdf_doc *pdf_doc_ptr, struct pdf_object *pdf_page)
      : datasets(&datasets_ref),
        query_key(&key),
        config(&cfg),
        pdf(pdf_doc_ptr),
        page(pdf_page),
        canvas_width(cfg.effective_width()),
        canvas_height(cfg.effective_height()),
        px0(100.0),
        px1(cfg.effective_width() - 40.0),
        py0(50.0),
        py1(cfg.effective_height() - 90.0),
        bg_col(detail::parseHexColor(cfg.bg_color())),
        axis_col(detail::parseHexColor(cfg.axis_color())),
        grid_col(detail::parseHexColor(cfg.grid_color())),
        text_col(detail::parseHexColor(cfg.text_color())) {}

  double toPdfY(double y_svg) const {
    return canvas_height - y_svg;
  }

  bool prepareData() {
    for (const auto &dataset : *datasets) {
      const auto &partials =
          dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
      auto iterator = partials.find(*query_key);
      if (iterator == partials.end()) {
        continue;
      }
      const auto &x_values = dataset.hist->bins;
      const auto &y_values = iterator->second;
      if (!x_values.empty()) {
        raw_x_min = std::min(raw_x_min, static_cast<double>(x_values.front()));
        raw_x_max = std::max(raw_x_max, static_cast<double>(x_values.back()));
      }
      for (real_t y_value : y_values) {
        raw_y_max = std::max(raw_y_max, static_cast<double>(y_value));
        raw_y_min = std::min(raw_y_min, static_cast<double>(y_value));
      }
    }

    if (raw_x_max <= raw_x_min || raw_y_max <= raw_y_min) {
      return false;
    }

    double y_padding = (raw_y_max - raw_y_min) * 0.05;
    raw_y_max += y_padding;

    xScale = detail::NiceScale(detail::DataRange{.min = raw_x_min, .max = raw_x_max}, 11);
    yScale = detail::NiceScale(detail::DataRange{.min = raw_y_min, .max = raw_y_max}, 8);
    return true;
  }

  void drawGridAndTicks() const {
    for (double y_val : yScale.ticks) {
      double spy = detail::mapValue(y_val, yScale.min, yScale.max, py1, py0);
      double pdf_y = toPdfY(spy);
      if (config->show_grid) {
        pdf_add_line(pdf, page, static_cast<float>(px0), static_cast<float>(pdf_y), static_cast<float>(px1),
                     static_cast<float>(pdf_y), 0.8F, grid_col);
      }
      pdf_add_line(pdf, page, static_cast<float>(px0 - 8.0), static_cast<float>(pdf_y), static_cast<float>(px0),
                   static_cast<float>(pdf_y), 1.5F, axis_col);
      detail::drawPdfText(pdf, page, detail::fmtScientificPdf(y_val), {px0 - 15.0, pdf_y - 5.0 * config->font_scale},
                          20.0 * config->font_scale, TextAnchor::End, text_col);
    }

    for (double x_val : xScale.ticks) {
      double spx = detail::mapValue(x_val, xScale.min, xScale.max, px0, px1);
      double pdf_y = toPdfY(py1);
      if (config->show_grid) {
        pdf_add_line(pdf, page, static_cast<float>(spx), static_cast<float>(toPdfY(py0)), static_cast<float>(spx),
                     static_cast<float>(pdf_y), 0.8F, grid_col);
      }
      pdf_add_line(pdf, page, static_cast<float>(spx), static_cast<float>(pdf_y - 8.0), static_cast<float>(spx),
                   static_cast<float>(pdf_y), 1.5F, axis_col);
      detail::drawPdfText(pdf, page, detail::fmtScientificPdf(x_val), {spx, pdf_y - 25.0 * config->font_scale},
                          20.0 * config->font_scale, TextAnchor::Middle, text_col);
    }

    pdf_add_rectangle(pdf, page, static_cast<float>(px0), static_cast<float>(toPdfY(py1)), static_cast<float>(px1 - px0),
                      static_cast<float>(py1 - py0), 1.5F, axis_col);
  }

  void drawAreaFills() const {
    if (config->fill_area) {
      std::size_t fill_ci = 0;
      for (const auto &dataset : *datasets) {
        const auto &partials =
            dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
        auto iterator = partials.find(*query_key);
        if (iterator != partials.end()) {
          const auto &x_values = dataset.hist->bins;
          const auto &y_values = iterator->second;
          uint32_t col = detail::parseHexColor(detail::color(fill_ci, config->palette));
          uint32_t shade_col = detail::blendColor({.fg = col, .bg = bg_col}, 0.20);
          std::size_t n_points = std::min(x_values.size(), y_values.size());
          if (n_points > 1) {
            std::vector<float> x_coords;
            std::vector<float> y_coords;
            x_coords.reserve(n_points + 2);
            y_coords.reserve(n_points + 2);

            double sx_start = detail::mapValue(x_values[0], xScale.min, xScale.max, px0, px1);
            x_coords.push_back(static_cast<float>(sx_start));
            y_coords.push_back(static_cast<float>(toPdfY(py1)));

            for (std::size_t i = 0; i < n_points; ++i) {
              double sp_x = detail::mapValue(x_values[i], xScale.min, xScale.max, px0, px1);
              double sp_y = detail::mapValue(y_values[i], yScale.min, yScale.max, py1, py0);
              double sp_y_clamped = std::min(sp_y, py1);

              x_coords.push_back(static_cast<float>(sp_x));
              y_coords.push_back(static_cast<float>(toPdfY(sp_y_clamped)));
            }

            double sx_end = detail::mapValue(x_values[n_points - 1], xScale.min, xScale.max, px0, px1);
            x_coords.push_back(static_cast<float>(sx_end));
            y_coords.push_back(static_cast<float>(toPdfY(py1)));

            pdf_add_filled_polygon(pdf, page, x_coords.data(), y_coords.data(), static_cast<int>(x_coords.size()), 0.0F,
                                   shade_col);
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
      auto partial_iter = partials.find(*query_key);
      if (partial_iter == partials.end()) {
        continue;
      }

      const auto &x_values = dataset.hist->bins;
      const auto &y_values = partial_iter->second;
      uint32_t col = detail::parseHexColor(detail::color(color_idx++, config->palette));
      legend_items.emplace_back(dataset.label, col);

      std::size_t n_points = std::min(x_values.size(), y_values.size());
      for (size_t i = 1; i < n_points; ++i) {
        double x_1 = detail::mapValue(x_values[i - 1], xScale.min, xScale.max, px0, px1);
        double y_1 = detail::mapValue(y_values[i - 1], yScale.min, yScale.max, py1, py0);
        double x_2 = detail::mapValue(x_values[i], xScale.min, xScale.max, px0, px1);
        double y_2 = detail::mapValue(y_values[i], yScale.min, yScale.max, py1, py0);

        pdf_add_line(pdf, page, static_cast<float>(x_1), static_cast<float>(toPdfY(y_1)), static_cast<float>(x_2),
                     static_cast<float>(toPdfY(y_2)), static_cast<float>(config->line_width), col);
      }
    }
  }

  void drawMarkers() const {
    if (config->show_markers) {
      std::size_t marker_color_idx = 0;
      for (const auto &dataset : *datasets) {
        const auto &partials =
            dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
        auto partial_iter = partials.find(*query_key);
        if (partial_iter != partials.end()) {
          const auto &x_values = dataset.hist->bins;
          const auto &y_values = partial_iter->second;
          uint32_t col = detail::parseHexColor(detail::color(marker_color_idx++, config->palette));
          std::size_t n_points = std::min(x_values.size(), y_values.size());
          for (std::size_t i = 0; i < n_points; ++i) {
            double sp_x = detail::mapValue(x_values[i], xScale.min, xScale.max, px0, px1);
            double sp_y = detail::mapValue(y_values[i], yScale.min, yScale.max, py1, py0);
            pdf_add_circle(pdf, page, static_cast<float>(sp_x), static_cast<float>(toPdfY(sp_y)), 3.5F, 0.0F, col,
                           col);
          }
        } else {
          marker_color_idx++;
        }
      }
    }
  }

  void drawLegend() const {
    if (config->show_legend) {
      double legend_x = px1 - 15.0;
      double legend_y = py0 + 25.0;
      for (auto it = legend_items.rbegin(); it != legend_items.rend(); ++it) {
        const auto &[label, col] = *it;
        double pdf_y = toPdfY(legend_y);
        pdf_add_line(pdf, page, static_cast<float>(legend_x - 40.0), static_cast<float>(pdf_y),
                     static_cast<float>(legend_x - 10.0), static_cast<float>(pdf_y), 4.0F, col);
        detail::drawPdfText(pdf, page, label, {legend_x - 45.0, pdf_y - 5.0 * config->font_scale}, 18.0 * config->font_scale,
                            TextAnchor::End, text_col);
        legend_y += 28.0;
      }
    }
  }

  void drawTitlesAndLabels() const {
    const auto &ref = *datasets->front().hist;
    std::string title_txt = *query_key + " Comparison";
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

    detail::drawPdfText(pdf, page, title_txt, {canvas_width / 2.0, toPdfY(30.0)}, 24.0 * config->font_scale, TextAnchor::Middle,
                        text_col);
    detail::drawPdfText(pdf, page, x_label, {(px0 + px1) / 2.0, toPdfY(py1 + 75.0)}, 28.0 * config->font_scale,
                        TextAnchor::Middle, text_col);
    detail::drawPdfText(pdf, page, y_label, {40.0, toPdfY((py0 + py1) / 2.0)}, 28.0 * config->font_scale, TextAnchor::Middle,
                        text_col, detail::kPi / 2.0);
  }
};

} // namespace detail

inline void renderHistogramAsPdf(const correlation::analysis::Histogram &hist, const std::string &filepath,
                                 const PlotConfig &config) {
  double canvas_width = config.effective_width();
  double canvas_height = config.effective_height();

  // Setup pdfgen
  pdf_info info = {};
  pdf_doc *pdf = pdf_create(static_cast<float>(canvas_width), static_cast<float>(canvas_height), &info);
  struct pdf_object *page = pdf_append_page(pdf);

  uint32_t bg_col = detail::parseHexColor(config.bg_color());
  uint32_t text_col = detail::parseHexColor(config.text_color());

  // Background
  pdf_add_filled_rectangle(pdf, page, 0.0F, 0.0F, static_cast<float>(canvas_width), static_cast<float>(canvas_height), 0.0F, bg_col,
                           PDF_TRANSPARENT);

  pdf_set_font(pdf, "Helvetica");

  auto toPdfY = [&](double y_svg) { return canvas_height - y_svg; };

  if (hist.bins.empty() || (hist.smoothed_partials.empty() && hist.partials.empty())) {
    detail::drawPdfText(pdf, page, "No data available", {canvas_width / 2.0, toPdfY(canvas_height / 2.0)}, 24.0 * config.font_scale,
                        TextAnchor::Middle, text_col);
    pdf_save(pdf, filepath.c_str());
    pdf_destroy(pdf);
    return;
  }

  detail::PdfHistogramRenderer renderer(hist, config, pdf, page);
  if (renderer.prepareData()) {
    renderer.drawGridAndTicks();
    renderer.drawEmphasisLines();
    renderer.drawAreaFills();
    auto legend_items = renderer.drawPolylines();
    renderer.drawMarkers();
    renderer.drawLegend(legend_items);
    renderer.drawTitlesAndLabels();
  }

  pdf_save(pdf, filepath.c_str());
  pdf_destroy(pdf);
}

struct ComparisonQuery {
  std::string key;
  std::string filepath;
};

inline void renderComparisonPdf(const std::vector<LabeledHistogram> &datasets, const ComparisonQuery &query,
                                const PlotConfig &config) {
  if (datasets.empty()) {
    return;
  }

  double canvas_width = config.effective_width();
  double canvas_height = config.effective_height();

  pdf_info info = {};
  pdf_doc *pdf = pdf_create(static_cast<float>(canvas_width), static_cast<float>(canvas_height), &info);
  struct pdf_object *page = pdf_append_page(pdf);

  uint32_t bg_col = detail::parseHexColor(config.bg_color());
  uint32_t text_col = detail::parseHexColor(config.text_color());

  // Background
  pdf_add_filled_rectangle(pdf, page, 0.0F, 0.0F, static_cast<float>(canvas_width), static_cast<float>(canvas_height), 0.0F, bg_col,
                           PDF_TRANSPARENT);

  pdf_set_font(pdf, "Helvetica");

  auto toPdfY = [&](double y_svg) { return canvas_height - y_svg; };

  detail::PdfComparisonRenderer renderer(datasets, query.key, config, pdf, page);
  if (renderer.prepareData()) {
    renderer.drawGridAndTicks();
    renderer.drawAreaFills();
    renderer.drawPolylines();
    renderer.drawMarkers();
    renderer.drawLegend();
    renderer.drawTitlesAndLabels();
  } else {
    detail::drawPdfText(pdf, page, "No comparison data", {canvas_width / 2.0, toPdfY(canvas_height / 2.0)}, 24.0 * config.font_scale,
                        TextAnchor::Middle, text_col);
  }

  pdf_save(pdf, query.filepath.c_str());
  pdf_destroy(pdf);
}

} // namespace correlation::plotters
