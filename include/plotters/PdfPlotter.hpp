/**
 * @file PdfPlotter.hpp
 * @brief PDF generator for distribution function histograms using pdfgen.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "plotters/SvgPlotter.hpp" // For PlotConfig
#include "plotters/pdfgen.h"

#include <algorithm>
#include <map>
#include <string>
#include <vector>

namespace correlation::plotters {

namespace detail {
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

struct PdfPadding {
  double left = 0.0;
  double right = 0.0;
  double top = 0.0;
  double bottom = 0.0;
};

struct PdfDimensions {
  double width = 0.0;
  double height = 0.0;
};

inline void drawPdfAxes(pdf_doc *pdf, struct pdf_object *page, const PdfPadding &padding, const PdfDimensions &dims,
                        uint32_t text_col) {
  // Y-axis
  pdf_add_line(pdf, page, static_cast<float>(padding.left), static_cast<float>(dims.height - padding.bottom),
               static_cast<float>(padding.left), static_cast<float>(padding.top), 2.0F, text_col);
  // X-axis
  pdf_add_line(pdf, page, static_cast<float>(padding.left), static_cast<float>(dims.height - padding.bottom),
               static_cast<float>(dims.width - padding.right), static_cast<float>(dims.height - padding.bottom), 2.0F,
               text_col);
}
} // namespace detail

inline void renderHistogramAsPdf(const correlation::analysis::Histogram &hist, const std::string &filepath,
                                 const PlotConfig &config) {

  // Setup pdfgen
  pdf_info info = {};
  pdf_doc *pdf = pdf_create(static_cast<float>(config.width), static_cast<float>(config.height), &info);
  struct pdf_object *page = pdf_append_page(pdf);

  uint32_t bg_col = (config.theme == PlotConfig::Theme::Dark) ? PDF_RGB(30, 30, 46) : PDF_RGB(255, 255, 255);
  uint32_t text_col = (config.theme == PlotConfig::Theme::Dark) ? PDF_RGB(205, 214, 244) : PDF_RGB(33, 33, 33);

  // Background
  pdf_add_rectangle(pdf, page, 0.0F, 0.0F, static_cast<float>(config.width), static_cast<float>(config.height), 0.0F,
                    bg_col);

  double pad_left = 80.0 * config.font_scale;
  double pad_right = 40.0;
  double pad_top = 80.0 * config.font_scale;
  double pad_bottom = 60.0 * config.font_scale;

  detail::drawPdfAxes(pdf, page,
                      detail::PdfPadding{.left = pad_left, .right = pad_right, .top = pad_top, .bottom = pad_bottom},
                      detail::PdfDimensions{.width = config.width, .height = config.height}, text_col);
  pdf_set_font(pdf, "Helvetica");
  pdf_add_text(pdf, page, hist.title.c_str(), static_cast<float>(20.0 * config.font_scale),
               static_cast<float>(config.width / 2.0 - 50.0), 40.0F, text_col);
  pdf_add_text(pdf, page, hist.x_label.c_str(), static_cast<float>(14.0 * config.font_scale),
               static_cast<float>(config.width / 2.0 - 20.0), static_cast<float>(config.height - 20.0), text_col);

  if (hist.bins.empty() || hist.partials.empty()) {
    pdf_save(pdf, filepath.c_str());
    pdf_destroy(pdf);
    return;
  }

  double min_x = hist.bins.front();
  double max_x = hist.bins.back();
  double max_y = 0.0;
  double min_y = 0.0;

  // We don't need getPalette, we will use color() directly
  const auto &partials = hist.smoothed_partials.empty() ? hist.partials : hist.smoothed_partials;

  for (const auto &[key, vals] : partials) {
    for (double val : vals) {
      max_y = std::max(max_y, val);
      min_y = std::min(min_y, val);
    }
  }

  if (max_y == min_y) {
    max_y += 1.0;
  }
  if (max_x == min_x) {
    max_x += 1.0;
  }

  auto mapX = [&](double x_val) {
    return pad_left + (x_val - min_x) / (max_x - min_x) * (config.width - pad_left - pad_right);
  };
  auto mapY = [&](double y_val) {
    return config.height - pad_bottom - (y_val - min_y) / (max_y - min_y) * (config.height - pad_top - pad_bottom);
  };

  int color_idx = 0;
  for (const auto &[key, vals] : partials) {
    uint32_t col = detail::parseHexColor(detail::color(color_idx, config.palette));
    for (size_t i = 1; i < hist.bins.size(); ++i) {
      double x_1 = mapX(hist.bins[i - 1]);
      double y_1 = mapY(vals[i - 1]);
      double x_2 = mapX(hist.bins[i]);
      double y_2 = mapY(vals[i]);
      pdf_add_line(pdf, page, static_cast<float>(x_1), static_cast<float>(y_1), static_cast<float>(x_2),
                   static_cast<float>(y_2), static_cast<float>(config.line_width), col);
    }
    color_idx++;
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
  // Simple fallback for comparison PDF
  if (datasets.empty()) {
    return;
  }
  renderHistogramAsPdf(*datasets[0].hist, query.filepath, config);
}

} // namespace correlation::plotters
