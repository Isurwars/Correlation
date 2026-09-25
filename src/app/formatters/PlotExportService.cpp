/**
 * @file PlotExportService.cpp
 * @brief Implementation of decoupled plot export service.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/PlotExportService.hpp"
#include "plotters/PdfPlotter.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>

namespace correlation::app {

namespace {

[[nodiscard]] bool isPdfExtension(const std::string &filepath) {
  std::string ext = std::filesystem::path(filepath).extension().string();
  std::ranges::transform(ext, ext.begin(),
                         [](unsigned char chr) { return static_cast<char>(std::tolower(chr)); });
  return ext == ".pdf";
}

[[nodiscard]] std::expected<void, std::string> writeStringToFile(const std::string &filepath,
                                                                 const std::string &content) {
  std::ofstream out(filepath);
  if (!out.is_open()) {
    return std::unexpected("Failed to open file for writing: " + filepath);
  }
  out << content;
  if (!out.good()) {
    return std::unexpected("I/O error occurred while writing: " + filepath);
  }
  return {};
}

} // namespace

std::string PlotExportService::getComparisonKey(const correlation::analysis::Histogram *hist) {
  if (hist == nullptr) {
    return "Total";
  }
  const auto &partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;
  if (!partials.empty() && !partials.contains("Total")) {
    return partials.begin()->first;
  }
  return "Total";
}

std::expected<void, std::string>
PlotExportService::exportHistogram(const std::string &filepath,
                                   const correlation::analysis::Histogram &hist,
                                   const correlation::plotters::PlotConfig &config,
                                   const std::map<std::string, real_t> &ashcroft_weights) {
  if (isPdfExtension(filepath)) {
    try {
      correlation::plotters::renderHistogramAsPdf(hist, filepath, config);
      return {};
    } catch (const std::exception &ex) {
      return std::unexpected(std::string("PDF export failed: ") + ex.what());
    }
  }

  const std::string svg =
      correlation::plotters::renderHistogramAsSvg(hist, config, {}, ashcroft_weights);
  return writeStringToFile(filepath, svg);
}

std::expected<void, std::string> PlotExportService::exportComparison(
    const std::string &filepath, std::span<const correlation::plotters::LabeledHistogram> datasets,
    const std::string &comparison_key, const correlation::plotters::PlotConfig &config) {
  const std::vector<correlation::plotters::LabeledHistogram> dataset_vec(datasets.begin(),
                                                                         datasets.end());
  if (isPdfExtension(filepath)) {
    try {
      correlation::plotters::renderComparisonPdf(
          dataset_vec, {.key = comparison_key, .filepath = filepath}, config);
      return {};
    } catch (const std::exception &ex) {
      return std::unexpected(std::string("Comparison PDF export failed: ") + ex.what());
    }
  }

  const std::string svg =
      correlation::plotters::renderComparisonSvg(dataset_vec, comparison_key, config);
  return writeStringToFile(filepath, svg);
}

} // namespace correlation::app
