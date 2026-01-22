// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
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

#include <highfive/highfive.hpp>

FileWriter::FileWriter(const DistributionFunctions &df) : df_(df) {}

void FileWriter::writeHistogramToCSV(const std::string &filename,
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
  for (const auto &[key, val] : hist.partials) {
    raw_keys.push_back(key);
  }
  std::sort(raw_keys.begin(), raw_keys.end());

  std::vector<std::string> smoothed_keys;
  for (const auto &[key, val] : hist.smoothed_partials) {
    smoothed_keys.push_back(key);
  }
  std::sort(smoothed_keys.begin(), smoothed_keys.end());

  // --- Write Header ---
  file << hist.bin_label;
  // Write headers for raw data
  for (const auto &key : raw_keys) {
    file << "," << key;
  }
  // Write headers for smoothed data
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

void FileWriter::writeAllCSVs(const std::string &base_path,
                              bool write_smoothed) const {
  const std::map<std::string, std::string> file_map = {
      {"g(r)", "_g.csv"},       {"J(r)", "_J.csv"}, {"G(r)", "__G.csv"},
      {"f(theta)", "_PAD.csv"}, {"S(Q)", "_S.csv"}, {"XRD", "_XRD.csv"},
      {"CN", "_CN.csv"}};

  for (const auto &[name, suffix] : file_map) {
    try {
      const auto &hist = df_.getHistogram(name);

      std::string filename = base_path + suffix;
      writeHistogramToCSV(filename, hist);

    } catch (const std::out_of_range &) {
      // This is not an error; it just means the histogram wasn't calculated.
      // We can silently skip it.
    } catch (const std::exception &e) {
      std::cerr << "Error writing file for '" << name << "': " << e.what()
                << std::endl;
    }
  }
}

void FileWriter::writeHDF(const std::string &filename) const {
  try {
    HighFive::File file(filename, HighFive::File::ReadWrite |
                                      HighFive::File::Create |
                                      HighFive::File::Truncate);

    for (const auto &[name, hist] : df_.getAllHistograms()) {
      if (hist.partials.empty() || hist.bins.empty()) {
        continue;
      }

      // Create a group for each histogram (e.g., "g(r)")
      // Replace '/' with '_' in name to avoid HDF5 path issues
      std::string group_name = name;
      std::replace(group_name.begin(), group_name.end(), '/', '_');
      HighFive::Group group = file.createGroup(group_name);

      // Store bins as a dataset
      group.createDataSet("bins", hist.bins);
      // Add attribute for bin label
      group.createAttribute<std::string>("bin_label",
                                         HighFive::DataSpace::From(hist.bin_label))
          .write(hist.bin_label);

      // Store raw partials
      HighFive::Group raw_group = group.createGroup("raw");
      for (const auto &[key, data] : hist.partials) {
        raw_group.createDataSet(key, data);
      }

      // Store smoothed partials if any
      if (!hist.smoothed_partials.empty()) {
        HighFive::Group smoothed_group = group.createGroup("smoothed");
        for (const auto &[key, data] : hist.smoothed_partials) {
          smoothed_group.createDataSet(key, data);
        }
      }
    }
  } catch (const HighFive::Exception &err) {
    throw std::runtime_error("HDF5 Error: " + std::string(err.what()));
  }
}
