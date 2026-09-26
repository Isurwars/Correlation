/**
 * @file PlotTableFormatter.cpp
 * @brief Implementation of PlotTableFormatter.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/formatters/PlotTableFormatter.hpp"

#include <algorithm>
#include <format>

namespace correlation::app {

namespace {

[[nodiscard]] std::shared_ptr<slint::VectorModel<slint::SharedString>>
buildHeaders(const correlation::analysis::Histogram *hist, const std::vector<std::string> &keys) {
  auto headers = std::make_shared<slint::VectorModel<slint::SharedString>>();
  if (hist == nullptr) {
    return headers;
  }

  std::string x_header = hist->x_label;
  if (!hist->x_unit.empty()) {
    x_header += " (" + hist->x_unit + ")";
  }
  headers->push_back(slint::SharedString(x_header));
  for (const auto &key : keys) {
    headers->push_back(slint::SharedString(key));
  }
  return headers;
}

[[nodiscard]] std::shared_ptr<slint::VectorModel<TableRow>>
buildRows(const correlation::analysis::Histogram *hist,
          const std::map<std::string, std::vector<real_t>> &partials,
          const std::vector<std::string> &keys) {
  auto rows = std::make_shared<slint::VectorModel<TableRow>>();
  if (hist == nullptr || hist->bins.empty()) {
    return rows;
  }

  const size_t num_bins = hist->bins.size();
  for (size_t i = 0; i < num_bins; ++i) {
    TableRow row;
    auto row_values = std::make_shared<slint::VectorModel<slint::SharedString>>();
    row_values->push_back(slint::SharedString(std::format("{:.4f}", hist->bins[i])));

    for (const auto &key : keys) {
      real_t val = 0.0;
      auto iter = partials.find(key);
      if (iter != partials.end() && i < iter->second.size()) {
        val = iter->second[i];
      }
      row_values->push_back(slint::SharedString(std::format("{:.6g}", val)));
    }

    row.values = row_values;
    rows->push_back(row);
  }
  return rows;
}

} // namespace

std::vector<std::string>
PlotTableFormatter::extractSortedPartialKeys(const correlation::analysis::Histogram *hist) {
  if (hist == nullptr) {
    return {};
  }

  const auto &partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;

  std::vector<std::string> keys;
  for (const auto &pair : partials) {
    if (pair.first != "Total") {
      keys.push_back(pair.first);
    }
  }
  std::ranges::sort(keys);
  if (partials.contains("Total")) {
    keys.insert(keys.begin(), "Total");
  }
  return keys;
}

PlotTableModel PlotTableFormatter::formatTable(const correlation::analysis::Histogram *hist) {
  if (hist == nullptr) {
    return {
        .headers = std::make_shared<slint::VectorModel<slint::SharedString>>(),
        .rows = std::make_shared<slint::VectorModel<TableRow>>(),
    };
  }

  const auto &partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;
  const std::vector<std::string> keys = extractSortedPartialKeys(hist);

  return {
      .headers = buildHeaders(hist, keys),
      .rows = buildRows(hist, partials, keys),
  };
}

} // namespace correlation::app
