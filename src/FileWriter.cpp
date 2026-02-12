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

#include <hdf5_hl.h>
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
      {"g(r)", "_g.csv"},
      {"J(r)", "_J.csv"},
      {"G(r)", "__G.csv"},
      {"f(theta)", "_PAD.csv"},
      {"S(Q)", "_S.csv"},
      {"XRD", "_XRD.csv"},
      {"CN", "_CN.csv"},
      {"VACF", "_VACF.csv"},
      {"Normalized VACF", "_VACF_norm.csv"},
      {"VDOS", "_VDOS.csv"}};

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
      {"G(r)", {"r (Å)", "Å", "Å^-1", "Reduced Radial Distribution Function"}},
      {"f(theta)",
       {"theta (angle)", "Degrees", "degree^-1", "Bond Angle Distribution"}},
      {"S(Q)", {"q (Å^-1)", "Å^-1", "arbitrary units", "Structure Factor"}},
      {"XRD",
       {"2theta", "Degrees (2theta)", "Intensity",
        "X-Ray Diffraction Pattern"}},
      {"CN", {"counts", "neighbors", "Count", "Coordination Number"}},
      {"VACF",
       {"Time (fs)", "fs", "Å^2/fs^2", "Velocity Autocorrelation Function"}},
      {"Normalized VACF",
       {"Time (fs)", "fs", "normalized",
        "Normalized Velocity Autocorrelation Function"}},
      {"VDOS",
       {"Frequency (THz)", "THz", "arbitrary units",
        "Vibrational Density of States"}}};

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
      props.add(HighFive::LinkCreationOrder(HighFive::CreationOrder::Tracked |
                                            HighFive::CreationOrder::Indexed));

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
        group
            .createAttribute<std::string>(
                "description", HighFive::DataSpace::From(description))
            .write(description);
      }

      // Prepare headers and data keys
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
        headers.insert(headers.end(), smoothed_keys.begin(),
                       smoothed_keys.end());
      }

      // Determine dimensions
      size_t n_rows = hist.bins.size();
      size_t n_cols = headers.size();

      for (size_t col = 0; col < n_cols; ++col) {
        std::string col_name = headers[col];
        // Sanitize column name
        std::replace(col_name.begin(), col_name.end(), '(', '_');
        std::replace(col_name.begin(), col_name.end(), ')', '_');
        std::replace(col_name.begin(), col_name.end(), '/', '_');
        std::replace(col_name.begin(), col_name.end(), ' ', '_');

        // Add prefix
        std::stringstream ss;
        ss << std::setw(2) << std::setfill('0') << col << "_" << col_name;
        std::string dataset_name = ss.str();

        std::vector<double> col_data(n_rows);
        if (col == 0) {
          col_data = hist.bins;
        } else if (col <= raw_keys.size()) {
          col_data = hist.partials.at(raw_keys[col - 1]);
        } else {
          std::string original_key = smoothed_keys[col - 1 - raw_keys.size()];
          // remove "_smoothed" suffix for lookup
          original_key = original_key.substr(0, original_key.size() - 9);
          col_data = hist.smoothed_partials.at(original_key);
        }

        HighFive::DataSet ds = group.createDataSet(dataset_name, col_data);

        // Add attributes
        if (col == 0) {
          ds.createAttribute<std::string>("Long Name",
                                          HighFive::DataSpace::From(dim_label))
              .write(dim_label);
          ds.createAttribute<std::string>("Units",
                                          HighFive::DataSpace::From(bin_unit))
              .write(bin_unit);
          ds.createAttribute<std::string>("Comments",
                                          HighFive::DataSpace::From(dim_label))
              .write(dim_label);

          // Make this a dimension scale
          herr_t status = H5DSset_scale(ds.getId(), dim_label.c_str());
          if (status < 0) {
            std::cerr << "Warning: Failed to set dimension scale for "
                      << dataset_name << std::endl;
          }
        } else {
          ds.createAttribute<std::string>(
                "Long Name", HighFive::DataSpace::From(headers[col]))
              .write(headers[col]);
          ds.createAttribute<std::string>("Units",
                                          HighFive::DataSpace::From(data_unit))
              .write(data_unit);
          ds.createAttribute<std::string>(
                "Comments", HighFive::DataSpace::From(headers[col]))
              .write(headers[col]);

          // Attach the bin dimension scale (col 0) to this dataset
          // Re-open the bin dataset to get its ID.
          // Since we are in a loop, we can construct the name of the 0th
          // dataset.
          std::stringstream ss_bin;
          ss_bin << std::setw(2) << std::setfill('0') << 0 << "_" << headers[0];
          // Actually, we can just keep the hid_t of the bin dataset if we are
          // careful, but HighFive::DataSet might close it when it goes out of
          // scope if we don't keep it alive. Easier to just open it again or
          // store the ID if we can ensure it's valid. Given the scope, ds for
          // col=0 is gone. Let's retrieve it from the group.

          // Re-construct the bin dataset name (it's always index 0)
          // We need to sanitize the bin label same as we did before
          std::string bin_col_name = headers[0];
          std::replace(bin_col_name.begin(), bin_col_name.end(), '(', '_');
          std::replace(bin_col_name.begin(), bin_col_name.end(), ')', '_');
          std::replace(bin_col_name.begin(), bin_col_name.end(), '/', '_');
          std::replace(bin_col_name.begin(), bin_col_name.end(), ' ', '_');

          std::stringstream ss0;
          ss0 << "00_" << bin_col_name;
          std::string bin_ds_name = ss0.str();

          if (group.exist(bin_ds_name)) {
            HighFive::DataSet bin_ds = group.getDataSet(bin_ds_name);
            herr_t status = H5DSattach_scale(ds.getId(), bin_ds.getId(), 0);
            if (status < 0) {
              std::cerr << "Warning: Failed to attach dimension scale for "
                        << dataset_name << std::endl;
            }
          }
        }
      }
    }
  } catch (const HighFive::Exception &err) {
    throw std::runtime_error("HDF5 Error: " + std::string(err.what()));
  }
}
