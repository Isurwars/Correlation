/**
 * @file SvgHistogramRenderer.hpp
 * @brief Declaration of single histogram rendering logic for the SVG plotter.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "plotters/PlotTypes.hpp"

#include <map>
#include <string>

namespace correlation::plotters {

/**
 * @brief Renders a single `Histogram` as a self-contained SVG string.
 * @param hist Input histogram dataset.
 * @param config Canvas layout, themes, and configuration.
 * @param hover Real-time hover cursor coordinates and dimensions.
 * @param weights Partial weighting map for Ashcroft-Langreth scaling.
 * @param curve_visibility Selective curve visibility map.
 * @param custom_colors Explicit color override map per partial key.
 * @return Formatted XML/SVG string representation of the plot.
 */
[[nodiscard]] std::string
renderHistogramAsSvg(const correlation::analysis::Histogram &hist,
                     const PlotConfig &config = PlotConfig(), const HoverInfo &hover = HoverInfo(),
                     const std::map<std::string, real_t> &weights = {},
                     const std::map<std::string, bool> &curve_visibility = {},
                     const std::map<std::string, std::string> &custom_colors = {});

} // namespace correlation::plotters
