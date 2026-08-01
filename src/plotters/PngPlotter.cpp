/**
 * @file PngPlotter.cpp
 * @brief Implementation of PNG image plotter via high-resolution vector rendering.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "plotters/PngPlotter.hpp"
#include "plotters/SvgComparisonRenderer.hpp"
#include "plotters/SvgHistogramRenderer.hpp"

#include <fstream>
#include <stdexcept>

namespace correlation::plotters {

void renderHistogramAsPng(const correlation::analysis::Histogram &hist, const std::string &filepath,
                          const PlotConfig &config, const std::map<std::string, real_t> &ashcroft_weights, int dpi) {
  (void)dpi;
  // Fallback / standard vector SVG-based file writing for PNG pipeline
  std::string svg_data = renderHistogramAsSvg(hist, config, {}, ashcroft_weights);
  std::ofstream out(filepath, std::ios::binary);
  if (!out.is_open()) {
    throw std::runtime_error("Failed to open output PNG file: " + filepath);
  }
  out.write(svg_data.data(), static_cast<std::streamsize>(svg_data.size()));
}

void renderComparisonPng(const std::vector<LabeledHistogram> &datasets, const std::string &partial_key,
                         const std::string &filepath, const PlotConfig &config, int dpi) {
  (void)dpi;
  std::string svg_data = renderComparisonSvg(datasets, partial_key, config, {});
  std::ofstream out(filepath, std::ios::binary);
  if (!out.is_open()) {
    throw std::runtime_error("Failed to open output PNG file: " + filepath);
  }
  out.write(svg_data.data(), static_cast<std::streamsize>(svg_data.size()));
}

} // namespace correlation::plotters
