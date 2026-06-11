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

#include <format>
#include <map>
#include <string>
#include <vector>

namespace correlation::plotters {

namespace detail {
// Extract color from string "#RRGGBB"
inline uint32_t parseHexColor(const std::string &hex) {
  if (hex.length() == 7 && hex[0] == '#') {
    uint32_t r = std::stoul(hex.substr(1, 2), nullptr, 16);
    uint32_t g = std::stoul(hex.substr(3, 2), nullptr, 16);
    uint32_t b = std::stoul(hex.substr(5, 2), nullptr, 16);
    return PDF_RGB(r, g, b);
  }
  return PDF_BLACK;
}

inline void drawPdfAxes(pdf_doc *pdf, struct pdf_object *page, double pad_left, double pad_right, double pad_top, double pad_bottom, double width, double height, uint32_t text_col) {
  // Y-axis
  pdf_add_line(pdf, page, pad_left, height - pad_bottom, pad_left, pad_top, 2.0, text_col);
  // X-axis
  pdf_add_line(pdf, page, pad_left, height - pad_bottom, width - pad_right, height - pad_bottom, 2.0, text_col);
}
} // namespace detail

inline void renderHistogramAsPdf(const correlation::analysis::Histogram &hist, const std::string &filepath,
                                 const PlotConfig &config) {

  // Setup pdfgen
  pdf_info info = {};
  pdf_doc *pdf = pdf_create(config.width, config.height, &info);
  struct pdf_object *page = pdf_append_page(pdf);
  
  uint32_t bg_col = (config.theme == PlotConfig::Theme::Dark) ? PDF_RGB(30, 30, 46) : PDF_RGB(255, 255, 255);
  uint32_t text_col = (config.theme == PlotConfig::Theme::Dark) ? PDF_RGB(205, 214, 244) : PDF_RGB(33, 33, 33);
  
  // Background
  pdf_add_rectangle(pdf, page, 0, 0, config.width, config.height, 0, bg_col);

  double pad_left = 80.0 * config.font_scale;
  double pad_right = 40.0;
  double pad_top = 80.0 * config.font_scale;
  double pad_bottom = 60.0 * config.font_scale;

  detail::drawPdfAxes(pdf, page, pad_left, pad_right, pad_top, pad_bottom, config.width, config.height, text_col);
  pdf_set_font(pdf, "Helvetica");
  pdf_add_text(pdf, page, hist.title.c_str(), 20 * config.font_scale, config.width/2.0 - 50.0, 40, text_col);
  pdf_add_text(pdf, page, hist.x_label.c_str(), 14 * config.font_scale, config.width/2.0 - 20.0, config.height - 20, text_col);

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
    for (double v : vals) {
      if (v > max_y) max_y = v;
      if (v < min_y) min_y = v;
    }
  }

  if (max_y == min_y) max_y += 1.0;
  if (max_x == min_x) max_x += 1.0;

  auto mapX = [&](double x) { return pad_left + (x - min_x) / (max_x - min_x) * (config.width - pad_left - pad_right); };
  auto mapY = [&](double y) { return config.height - pad_bottom - (y - min_y) / (max_y - min_y) * (config.height - pad_top - pad_bottom); };

  int color_idx = 0;
  for (const auto &[key, vals] : partials) {
    uint32_t col = detail::parseHexColor(detail::color(color_idx, config.palette));
    for (size_t i = 1; i < hist.bins.size(); ++i) {
      double x1 = mapX(hist.bins[i - 1]);
      double y1 = mapY(vals[i - 1]);
      double x2 = mapX(hist.bins[i]);
      double y2 = mapY(vals[i]);
      pdf_add_line(pdf, page, x1, y1, x2, y2, config.line_width, col);
    }
    color_idx++;
  }

  pdf_save(pdf, filepath.c_str());
  pdf_destroy(pdf);
}

inline void renderComparisonPdf(const std::vector<LabeledHistogram> &datasets, const std::string &key,
                                const std::string &filepath, const PlotConfig &config) {
  // Simple fallback for comparison PDF
  if (datasets.empty()) return;
  renderHistogramAsPdf(*datasets[0].hist, filepath, config);
}

} // namespace correlation::plotters
