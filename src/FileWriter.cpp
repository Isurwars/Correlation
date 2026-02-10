// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "FileWriter.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <highfive/highfive.hpp>
#include <H5DSpublic.h>

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
      {"CN", "_CN.csv"}, {"VACF", "_VACF.csv"},
      {"Normalized VACF", "_VACF_norm.csv"}};

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
  // Define metadata structure
  struct FunctionMetadata {
    std::string bin_label;
    std::string bin_unit;
    std::string data_unit;
    std::string description;
  };

  // Metadata mapping
  const std::map<std::string, FunctionMetadata> metadata_map = {
      {"g(r)", {"r (Å)", "Å", "Å^-1", "Radial Distribution Function"}},
      {"J(r)",
       {"r (Å)", "Å", "Å^-1", "Radial Distribution of Electron Density"}},
      {"G(r)",
       {"r (Å)", "Å", "Å^-1", "Reduced Radial Distribution Function"}},
      {"f(theta)", {"theta (angle)", "Degrees", "degree^-1", "Bond Angle Distribution"}},
      {"S(Q)", {"q (Å^-1)", "Å^-1", "arbitrary units", "Structure Factor"}},
      {"XRD", {"2theta", "Degrees (2theta)", "Intensity", "X-Ray Diffraction Pattern"}},
      {"CN", {"counts", "neighbors", "Count", "Coordination Number"}},
      {"VACF", {"Time (fs)", "fs", "Å^2/fs^2", "Velocity Autocorrelation Function"}},
      {"Normalized VACF",
       {"Time (fs)", "fs", "normalized", "Normalized Velocity Autocorrelation Function"}}};

  try {
    HighFive::File file(filename, HighFive::File::ReadWrite |
                                      HighFive::File::Create |
                                      HighFive::File::Truncate);

    for (const auto &[name, hist] : df_.getAllHistograms()) {
      if (hist.partials.empty() || hist.bins.empty()) {
        continue;
      }

      // Sanitize group name: replace '(' with '_' and remove ')'
      std::string group_name = name;
      std::replace(group_name.begin(), group_name.end(), '(', '_');
      group_name.erase(std::remove(group_name.begin(), group_name.end(), ')'),
                       group_name.end());
      
      // Also replace '/' with '_' just in case
      std::replace(group_name.begin(), group_name.end(), '/', '_');
      // Replace spaces with underscores
      std::replace(group_name.begin(), group_name.end(), ' ', '_');

      // Create a property list for group creation
      HighFive::GroupCreateProps props;
      props.add(HighFive::LinkCreationOrder(HighFive::CreationOrder::Tracked | HighFive::CreationOrder::Indexed));

      HighFive::Group group = file.createGroup(group_name, props);

      // Get metadata if available
      std::string bin_unit = "arbitrary units";
      std::string data_unit = "arbitrary units";
      std::string description = "";
      std::string dim_label = hist.bin_label.empty() ? "x" : hist.bin_label;

      if (metadata_map.count(name)) {
        const auto &meta = metadata_map.at(name);
        bin_unit = meta.bin_unit;
        data_unit = meta.data_unit;
        description = meta.description;
        if (!meta.bin_label.empty()) {
            dim_label = meta.bin_label;
        }
      }

      // Add description attribute to the group
      if (!description.empty()) {
        group.createAttribute<std::string>(
                 "description", HighFive::DataSpace::From(description))
            .write(description);
      }

      // Prepare headers
      std::vector<std::string> headers;
      // Start with bins
      headers.push_back(hist.bin_label);
      // Then raw data keys
      std::vector<std::string> raw_keys;
      for (const auto &[key, _] : hist.partials) {
        raw_keys.push_back(key);
      }
      std::sort(raw_keys.begin(), raw_keys.end());
      headers.insert(headers.end(), raw_keys.begin(), raw_keys.end());

      // Then smoothed data keys
      std::vector<std::string> smoothed_keys;
      if (!hist.smoothed_partials.empty()) {
        for (const auto &[key, _] : hist.smoothed_partials) {
          smoothed_keys.push_back(key + "_smoothed");
        }
        std::sort(smoothed_keys.begin(), smoothed_keys.end());
        headers.insert(headers.end(), smoothed_keys.begin(), smoothed_keys.end());
      }
      
      // Determine dimensions
      size_t n_rows = hist.bins.size();
      size_t n_cols = headers.size();
      
      // Create and fill 2D data structure
      std::vector<std::vector<double>> data(n_rows, std::vector<double>(n_cols));
      
      for (size_t i = 0; i < n_rows; ++i) {
          size_t col = 0;
          // Col 0: bins
          data[i][col++] = hist.bins[i];
          
          // Next cols: raw data
          for (const auto &key : raw_keys) {
              data[i][col++] = hist.partials.at(key)[i];
          }
          
          // Next cols: smoothed data
          if (!hist.smoothed_partials.empty()) {
              // Extract smoothed keys without suffix for lookup
              for (const auto &s_key : smoothed_keys) {
                  // remove "_smoothed" suffix
                  std::string original_key = s_key.substr(0, s_key.size() - 9); 
                  data[i][col++] = hist.smoothed_partials.at(original_key)[i];
              }
          }
      }
      
      // Create Dataset
      HighFive::DataSet ds = group.createDataSet("data", data);
      
      // Add attributes
      ds.createAttribute<std::string>("units", HighFive::DataSpace::From(data_unit)).write(data_unit);
      
      // Try to write headers as attribute - vector of strings support varies, let's try standard way
      try {
           ds.createAttribute<std::string>("column_names", HighFive::DataSpace::From(headers)).write(headers);
      } catch (const HighFive::Exception& e) {
           std::cerr << "Warning: Could not write column_names attribute for " << group_name << ": " << e.what() << std::endl;
      }
      
      // Add bin label attribute
      ds.createAttribute<std::string>("bin_label", HighFive::DataSpace::From(dim_label)).write(dim_label);
      ds.createAttribute<std::string>("bin_units", HighFive::DataSpace::From(bin_unit)).write(bin_unit);
    }
  } catch (const HighFive::Exception &err) {
    throw std::runtime_error("HDF5 Error: " + std::string(err.what()));
  }
}
