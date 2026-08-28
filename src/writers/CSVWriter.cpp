/**
 * @file CSVWriter.cpp
 * @brief Implementation of the CSV writer.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "writers/CSVWriter.hpp"
#include "writers/WriterFactory.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace correlation::writers {

// Automatic registration
const bool registered = WriterFactory::registerTypeSafe<CSVWriter>("CSVWriter");

void CSVWriter::writeAllCSVs(const std::string &base_path, const correlation::analysis::DistributionFunctions &dists,
                             bool /*write_smoothed*/) {
  const auto &all_histograms = dists.getAllHistograms();
  for (const auto &[name, hist] : all_histograms) {
    try {
      if (hist.partials.empty() || hist.bins.empty() || hist.file_suffix.empty()) {
        continue;
      }

      const std::string filename = base_path + hist.file_suffix + ".csv";

      // For normalized PAD/DAD, inject raw companion partials.
      const correlation::analysis::Histogram *raw_companion = nullptr;
      const std::string raw_key = name + "_raw";
      if ((name == "PAD" || name == "DAD") && all_histograms.contains(raw_key)) {
        raw_companion = &all_histograms.at(raw_key);
      }

      writeHistogramToCSV(filename, hist, raw_companion);

    } catch (const std::exception &e) {
      std::cerr << "Error writing file for '" << name << "': " << e.what() << "\n";
    }
  }
}

namespace {

struct ColumnDef {
  std::string name;
  std::string unit;
  std::string comment;
  const std::vector<correlation::real_t> *data{nullptr};
};

[[nodiscard]] std::vector<std::string>
getSortedKeys(const std::map<std::string, std::vector<correlation::real_t>> &map) {
  std::vector<std::string> keys;
  keys.reserve(map.size());
  for (const auto &elem : map) {
    keys.push_back(elem.first);
  }
  std::ranges::sort(keys);
  return keys;
}

[[nodiscard]] std::vector<ColumnDef> buildColumns(const correlation::analysis::Histogram &hist,
                                                  const correlation::analysis::Histogram *raw_companion) {
  std::vector<ColumnDef> cols;

  const std::string raw_unit =
      (raw_companion != nullptr && !raw_companion->y_unit.empty()) ? raw_companion->y_unit : "counts";
  const std::string data_unit = hist.y_unit.empty() ? "arbitrary units" : hist.y_unit;

  if (raw_companion != nullptr) {
    for (const auto &key : getSortedKeys(raw_companion->partials)) {
      cols.push_back(
          {.name = key + "_raw", .unit = raw_unit, .comment = key + "_raw", .data = &raw_companion->partials.at(key)});
    }
  }

  for (const auto &key : getSortedKeys(hist.partials)) {
    cols.push_back({.name = key, .unit = data_unit, .comment = key, .data = &hist.partials.at(key)});
  }

  for (const auto &key : getSortedKeys(hist.smoothed_partials)) {
    cols.push_back({.name = key + "_smoothed",
                    .unit = data_unit,
                    .comment = key + "_smoothed",
                    .data = &hist.smoothed_partials.at(key)});
  }

  return cols;
}

void writeHeader(std::ostream &file, const correlation::analysis::Histogram &hist, const std::vector<ColumnDef> &cols) {
  const std::string bin_unit = hist.x_unit.empty() ? "arbitrary units" : hist.x_unit;
  const std::string description = hist.description.empty() ? "Data export" : hist.description;
  const std::string dim_label = hist.x_label.empty() ? "x" : hist.x_label;

  // Line 1: Long Name
  file << dim_label;
  for (const auto &col : cols) {
    file << "," << col.name;
  }
  file << '\n';

  // Line 2: Units
  file << bin_unit;
  for (const auto &col : cols) {
    file << "," << col.unit;
  }
  file << '\n';

  // Line 3: Comments
  file << description;
  for (const auto &col : cols) {
    file << "," << col.comment;
  }
  file << '\n';
}

void writeDataRows(std::ostream &file, const std::vector<correlation::real_t> &bins,
                   const std::vector<ColumnDef> &cols) {
  const size_t num_rows = bins.size();
  for (size_t i = 0; i < num_rows; ++i) {
    file << std::fixed << std::setprecision(5) << bins[i];
    for (const auto &col : cols) {
      file << "," << (*col.data)[i];
    }
    file << '\n';
  }
}

} // namespace

void CSVWriter::writeHistogramToCSV(const std::string &filename, const correlation::analysis::Histogram &hist,
                                    const correlation::analysis::Histogram *raw_companion) {
  if (hist.partials.empty() || hist.bins.empty()) {
    return;
  }

  std::ofstream file(filename);
  if (!file) {
    throw std::runtime_error("Failed to open file for writing: " + filename);
  }

  const auto columns = buildColumns(hist, raw_companion);
  writeHeader(file, hist, columns);
  writeDataRows(file, hist.bins, columns);
}

} // namespace correlation::writers
