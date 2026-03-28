// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "writers/CSVWriter.hpp"
#include "writers/WriterFactory.hpp"
#include "writers/WriterUtils.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Writer {

// Automatic registration
static bool registered =
    WriterFactory::instance().registerWriter(std::make_unique<CSVWriter>());

void CSVWriter::writeAllCSVs(const std::string &base_path,
                             const DistributionFunctions &df,
                             bool /*write_smoothed*/) const {
  const std::map<std::string, std::string> file_map = {
      {"g_r", "_g.csv"},
      {"J_r", "_J.csv"},
      {"G_r", "_G_reduced.csv"},
      {"BAD", "_PAD.csv"},
      {"DAD", "_DAD.csv"},
      {"RD", "_RD.csv"},
      {"fft_S_q", "_fft_S.csv"},
      {"XRD", "_XRD.csv"},
      {"CN", "_CN.csv"},
      {"VACF", "_VACF.csv"},
      {"Normalized VACF", "_VACF_norm.csv"},
      {"VDOS", "_VDOS.csv"},
      {"debye_S_q", "_Debye_S.csv"},
      {"S_q", "_S.csv"},
      {"MSD", "_MSD.csv"},
      {"D_eff", "_D_eff.csv"}};

  for (const auto &[name, suffix] : file_map) {
    try {
      const auto &hist = df.getHistogram(name);

      std::string filename = base_path + suffix;
      writeHistogramToCSV(filename, name, hist);

    } catch (const std::out_of_range &) {
      // This is not an error; it just means the histogram wasn't calculated.
      // We can silently skip it.
    } catch (const std::exception &e) {
      std::cerr << "Error writing file for '" << name << "': " << e.what()
                << std::endl;
    }
  }
}

void CSVWriter::writeHistogramToCSV(const std::string &filename,
                                    const std::string &name,
                                    const Histogram &hist) const {
  if (hist.partials.empty() || hist.bins.empty()) {
    return;
  }

  std::ofstream file(filename);
  if (!file) {
    throw std::runtime_error("Failed to open file for writing: " + filename);
  }

  // Get sorted keys for both raw and smoothed data
  std::vector<std::string> raw_keys;
  for (const auto &[k, v] : hist.partials) {
    raw_keys.push_back(k);
  }
  std::sort(raw_keys.begin(), raw_keys.end());

  std::vector<std::string> smoothed_keys;
  for (const auto &[k, v] : hist.smoothed_partials) {
    smoothed_keys.push_back(k);
  }
  std::sort(smoothed_keys.begin(), smoothed_keys.end());

  // Get metadata if available
  std::string bin_unit = "arbitrary units";
  std::string data_unit = "arbitrary units";
  std::string description = "Data export";
  std::string dim_label = hist.bin_label.empty() ? "x" : hist.bin_label;

  const auto &metadata_map = Correlation::WriterUtils::metadata_map;

  if (metadata_map.count(name)) {
    const auto &meta = metadata_map.at(name);
    bin_unit = meta.bin_unit;
    data_unit = meta.data_unit;
    description = meta.description;
    if (!meta.bin_label.empty()) {
      dim_label = meta.bin_label;
    }
  }

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
    std::string current_data_unit = data_unit;
    if (metadata_map.count(key)) {
      current_data_unit = metadata_map.at(key).bin_unit;
    }
    file << "," << current_data_unit;
  }
  for (const auto &key : smoothed_keys) {
    std::string current_data_unit = data_unit;
    if (metadata_map.count(key)) {
      current_data_unit = metadata_map.at(key).bin_unit;
    }
    file << "," << current_data_unit;
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

} // namespace Writer
