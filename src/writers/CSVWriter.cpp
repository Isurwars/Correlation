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
  for (const auto &[name, hist] : dists.getAllHistograms()) {
    try {
      if (hist.partials.empty() || hist.bins.empty() || hist.file_suffix.empty()) {
        continue;
      }

      std::string filename = base_path + hist.file_suffix + ".csv";
      writeHistogramToCSV(filename, hist);

    } catch (const std::exception &e) {
      std::cerr << "Error writing file for '" << name << "': " << e.what() << "\n";
    }
  }
}

void CSVWriter::writeHistogramToCSV(const std::string &filename, const correlation::analysis::Histogram &hist) {
  if (hist.partials.empty() || hist.bins.empty()) {
    return;
  }

  std::ofstream file(filename);
  if (!file) {
    throw std::runtime_error("Failed to open file for writing: " + filename);
  }

  // Get sorted keys for both raw and smoothed data
  std::vector<std::string> raw_keys;
  raw_keys.reserve(hist.partials.size());
  for (const auto &[key, value] : hist.partials) {
    raw_keys.push_back(key);
  }
  std::ranges::sort(raw_keys);

  std::vector<std::string> smoothed_keys;
  smoothed_keys.reserve(hist.smoothed_partials.size());
  for (const auto &[key, value] : hist.smoothed_partials) {
    smoothed_keys.push_back(key);
  }
  std::ranges::sort(smoothed_keys);

  // Get metadata if available
  std::string bin_unit = hist.x_unit.empty() ? "arbitrary units" : hist.x_unit;
  std::string data_unit = hist.y_unit.empty() ? "arbitrary units" : hist.y_unit;
  std::string description = hist.description.empty() ? "Data export" : hist.description;
  std::string dim_label = hist.x_label.empty() ? "x" : hist.x_label;

  // --- Write Header ---
  // Line 1: Long Name
  file << dim_label;
  for (const auto &key : raw_keys) {
    file << "," << key;
  }
  for (const auto &key : smoothed_keys) {
    file << "," << key << "_smoothed";
  }
  file << '\n';

  // Line 2: Units
  file << bin_unit;
  for (const auto &key : raw_keys) {
    file << "," << data_unit;
  }
  for (const auto &key : smoothed_keys) {
    file << "," << data_unit;
  }
  file << '\n';

  // Line 3: Comments
  file << description;
  for (const auto &key : raw_keys) {
    file << "," << key;
  }
  for (const auto &key : smoothed_keys) {
    file << "," << key << "_smoothed";
  }
  file << '\n';

  // --- Write Data Rows ---
  const size_t num_rows = hist.bins.size();
  for (size_t i = 0; i < num_rows; ++i) {
    file << std::fixed << std::setprecision(5) << hist.bins[i];

    // Write raw data
    for (const auto &key : raw_keys) {
      file << "," << hist.partials.at(key)[i];
    }

    // Write smoothed data
    for (const auto &key : smoothed_keys) {
      file << "," << hist.smoothed_partials.at(key)[i];
    }
    file << '\n';
  }
}

} // namespace correlation::writers
