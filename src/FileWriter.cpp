// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/FileWriter.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

FileWriter::FileWriter(const DistributionFunctions &df) : df_(df) {}

void FileWriter::writeHistogramToCSV(const std::string &filename,
                                     const Histogram &hist,
                                     bool use_smoothed_data) const {
  const auto &partials_to_write =
      use_smoothed_data ? hist.smoothed_partials : hist.partials;

  if (partials_to_write.empty() || hist.bins.empty()) {
    // Don't write files for empty or uncalculated histograms
    return;
  }

  std::ofstream file(filename);
  if (!file) {
    throw std::runtime_error("Failed to open file for writing: " + filename);
  }

  // Get a sorted list of keys to ensure a consistent column order in the CSV
  std::vector<std::string> keys;
  for (const auto &[key, val] : partials_to_write) {
    keys.push_back(key);
  }
  std::sort(keys.begin(), keys.end());

  // --- Write Header ---
  // The header is now generated dynamically from the keys of the data map.
  file << std::setw(12) << std::right << hist.bin_label << ",";
  for (size_t i = 0; i < keys.size(); ++i) {
    file << std::setw(12) << std::right << keys[i]
         << (i == keys.size() - 1 ? "" : ",");
  }
  file << '\n';

  // --- Write Data Rows ---
  // This approach iterates row-by-row, which is more memory-efficient.
  const size_t num_rows = hist.bins.size();
  for (size_t i = 0; i < num_rows; ++i) {
    file << std::fixed << std::setprecision(5) << std::setw(12) << std::right
         << hist.bins[i] << ",";
    for (size_t j = 0; j < keys.size(); ++j) {
      // Use .at() to access map elements by key
      file << std::setw(12) << std::right << partials_to_write.at(keys[j])[i]
           << (j == keys.size() - 1 ? "" : ",");
    }
    file << '\n';
  }
  std::cout << "Successfully wrote to file: " << filename << std::endl;
}

void FileWriter::writeAllCSVs(const std::string &base_path,
                              bool write_smoothed) const {
  // A simple map to associate internal histogram names with desired file
  // suffixes.
  const std::map<std::string, std::string> file_map = {
      {"g(r)", "_g.csv"},       {"J(r)", "_J.csv"}, {"G(r)", "__G.csv"},
      {"f(theta)", "_PAD.csv"}, {"S(Q)", "_S.csv"}, {"XRD", "_XRD.csv"},
      {"CN", "_CN.csv"}};

  for (const auto &[name, suffix] : file_map) {
    try {
      const auto &hist = df_.getHistogram(name);

      // Write the raw (non-smoothed) data
      std::string filename = base_path + suffix;
      writeHistogramToCSV(filename, hist, false);

      // If requested, write the smoothed data to a separate file
      if (write_smoothed && !hist.smoothed_partials.empty()) {
        std::string smoothed_suffix = suffix;
        // Insert "_smoothed" before the file extension (e.g., _g.csv ->
        // _g_smoothed.csv)
        smoothed_suffix.insert(suffix.rfind('.'), "_smoothed");
        std::string smoothed_filename = base_path + smoothed_suffix;
        writeHistogramToCSV(smoothed_filename, hist, true);
      }

    } catch (const std::out_of_range &) {
      // This is not an error; it just means the histogram wasn't calculated.
      // We can silently skip it.
    } catch (const std::exception &e) {
      std::cerr << "Error writing file for '" << name << "': " << e.what()
                << std::endl;
    }
  }
}
