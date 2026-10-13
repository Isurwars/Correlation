/**
 * @file PdfPlotter.hpp
 * @brief High-quality PDF export generator for distribution function histograms.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "plotters/PlotTypes.hpp"

#include <string>
#include <vector>

namespace correlation::plotters {

/**
 * @struct ComparisonQuery
 * @brief Parameters defining a multi-dataset comparison PDF query.
 */
struct ComparisonQuery {
  std::string key;
  std::string filepath;
};

/**
 * @brief Renders a single histogram dataset to a PDF document file.
 * @param hist Input histogram dataset.
 * @param filepath Destination filesystem path.
 * @param config Canvas layout, colors, and font configuration.
 */
void renderHistogramAsPdf(const correlation::analysis::Histogram &hist, const std::string &filepath,
                          const PlotConfig &config = PlotConfig());

/**
 * @brief Renders an overlaid multi-dataset comparison plot to a PDF document file.
 * @param datasets Collection of labeled histogram datasets.
 * @param query Comparison query specifying key and destination filepath.
 * @param config Canvas layout, colors, and font configuration.
 */
void renderComparisonPdf(const std::vector<LabeledHistogram> &datasets,
                         const ComparisonQuery &query, const PlotConfig &config = PlotConfig());

} // namespace correlation::plotters
