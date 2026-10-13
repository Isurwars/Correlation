/**
 * @file SvgComparisonRenderer.hpp
 * @brief Declaration of multi-histogram comparison rendering logic for the SVG plotter.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "plotters/PlotTypes.hpp"

#include <string>
#include <vector>

namespace correlation::plotters {

/**
 * @brief Renders multiple histograms overlaid on a single comparison plot as an SVG string.
 * @param datasets Collection of labeled histogram datasets to compare.
 * @param partial_key Partial key to visualize across datasets (default "Total").
 * @param config Canvas layout, themes, and configuration.
 * @param hover Real-time hover cursor coordinates and dimensions.
 * @return Formatted XML/SVG string representation of the overlaid comparison plot.
 */
[[nodiscard]] std::string renderComparisonSvg(const std::vector<LabeledHistogram> &datasets,
                                              const std::string &partial_key = "Total",
                                              const PlotConfig &config = PlotConfig(),
                                              const HoverInfo &hover = HoverInfo());

} // namespace correlation::plotters
