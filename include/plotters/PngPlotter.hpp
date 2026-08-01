/**
 * @file PngPlotter.hpp
 * @brief High-resolution PNG image plotter declaration.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "plotters/PlotTypes.hpp"

#include <map>
#include <string>
#include <vector>

namespace correlation::plotters {

/**
 * @brief Renders a single histogram dataset as a high-resolution PNG image.
 * @param hist Input histogram data.
 * @param filepath Output PNG file path.
 * @param config Visualization and canvas layout configuration.
 * @param ashcroft_weights Partial weighting map for Ashcroft-Langreth scaling.
 * @param dpi Target DPI resolution (e.g. 150, 300, 600).
 */
void renderHistogramAsPng(const correlation::analysis::Histogram &hist, const std::string &filepath,
                          const PlotConfig &config = {},
                          const std::map<std::string, real_t> &ashcroft_weights = {}, int dpi = 300);

/**
 * @brief Renders multiple comparison histograms as a high-resolution PNG image.
 * @param datasets Vector of labeled histogram comparison datasets.
 * @param partial_key Target partial histogram key (e.g. "Total", "Si-O").
 * @param filepath Output PNG file path.
 * @param config Visualization and canvas layout configuration.
 * @param dpi Target DPI resolution (e.g. 150, 300, 600).
 */
void renderComparisonPng(const std::vector<LabeledHistogram> &datasets, const std::string &partial_key,
                         const std::string &filepath, const PlotConfig &config = {}, int dpi = 300);

} // namespace correlation::plotters
