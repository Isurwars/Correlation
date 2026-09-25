/**
 * @file PlotExportService.hpp
 * @brief Decoupled plot export service for SVG and PDF generation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "plotters/SvgPlotter.hpp"

#include <expected>
#include <map>
#include <span>
#include <string>

namespace correlation::app {

/**
 * @class PlotExportService
 * @brief Stateless service executing plot file exports decoupled from UI frameworks.
 */
class PlotExportService {
public:
  /**
   * @brief Determines the comparison curve key for a given histogram (e.g. "Total" or first partial).
   * @param[in] hist Histogram to inspect.
   * @return Key string.
   */
  [[nodiscard]] static std::string getComparisonKey(const correlation::analysis::Histogram *hist);

  /**
   * @brief Exports a single histogram plot to file (SVG or PDF).
   * @param[in] filepath Target destination file path.
   * @param[in] hist Histogram to render.
   * @param[in] config Plot styling configuration.
   * @param[in] ashcroft_weights Element weights for S(Q) partials if applicable.
   * @return Expected void on success or error description string on failure.
   */
  [[nodiscard]] static std::expected<void, std::string>
  exportHistogram(const std::string &filepath, const correlation::analysis::Histogram &hist,
                  const correlation::plotters::PlotConfig &config,
                  const std::map<std::string, real_t> &ashcroft_weights = {});

  /**
   * @brief Exports a comparison overlay plot to file (SVG or PDF).
   * @param[in] filepath Target destination file path.
   * @param[in] datasets Comparison datasets to render.
   * @param[in] comparison_key Key for the compared curve.
   * @param[in] config Plot styling configuration.
   * @return Expected void on success or error description string on failure.
   */
  [[nodiscard]] static std::expected<void, std::string>
  exportComparison(const std::string &filepath,
                   std::span<const correlation::plotters::LabeledHistogram> datasets,
                   const std::string &comparison_key,
                   const correlation::plotters::PlotConfig &config);
};

} // namespace correlation::app
