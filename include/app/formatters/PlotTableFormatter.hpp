/**
 * @file PlotTableFormatter.hpp
 * @brief Formatter transforming Histogram analytical data into Slint table models.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "AppWindow.h"
#include "analysis/DistributionFunctions.hpp"
#include <slint.h>

#include <memory>
#include <string>
#include <vector>

namespace correlation::app {

/**
 * @struct PlotTableModel
 * @brief Container holding Slint VectorModels for table headers and rows.
 */
struct PlotTableModel {
  std::shared_ptr<slint::VectorModel<slint::SharedString>> headers;
  std::shared_ptr<slint::VectorModel<TableRow>> rows;
};

/**
 * @class PlotTableFormatter
 * @brief Stateless utility transforming Histogram data into Slint tabular models.
 */
class PlotTableFormatter {
public:
  /**
   * @brief Extracts sorted partial curve keys from a histogram, placing 'Total' first if present.
   * @param[in] hist Pointer to the analytical histogram.
   * @return Ordered list of column keys for partials.
   */
  [[nodiscard]] static std::vector<std::string>
  extractSortedPartialKeys(const correlation::analysis::Histogram *hist);

  /**
   * @brief Formats histogram bins and partial distributions into Slint table models.
   * @param[in] hist Pointer to the analytical histogram.
   * @return Formatted PlotTableModel containing headers and rows.
   */
  [[nodiscard]] static PlotTableModel formatTable(const correlation::analysis::Histogram *hist);
};

} // namespace correlation::app
