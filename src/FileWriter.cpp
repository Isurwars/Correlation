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
    std::string bin_label; // Added bin_label
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

      // Store raw partials
      for (const auto &[key, data] : hist.partials) {
        // Create (N, 1) datasource for column vector
        HighFive::DataSpace dataspace({data.size(), 1});
        HighFive::DataSet ds = group.createDataSet(key, dataspace, HighFive::create_datatype<double>());
        ds.write_raw(data.data()); // Use write_raw to bypass dimension check for (N) vs (N, 1)
        
        ds.createAttribute<std::string>("units",
                                        HighFive::DataSpace::From(data_unit))
            .write(data_unit);
      }

      // Store smoothed partials if any
      if (!hist.smoothed_partials.empty()) {
        for (const auto &[key, data] : hist.smoothed_partials) {
          std::string smoothed_key = key + "_smoothed";
          // Create (N, 1) datasource for column vector
          HighFive::DataSpace dataspace({data.size(), 1});
          HighFive::DataSet ds = group.createDataSet(smoothed_key, dataspace, HighFive::create_datatype<double>());
          ds.write_raw(data.data()); // Use write_raw

          ds.createAttribute<std::string>("units",
                                          HighFive::DataSpace::From(data_unit))
              .write(data_unit);
        }
      }

      // Store bins as a dataset - WRITTEN LAST to hopefully appear first in LIFO readers
      // Create (N, 1) datasource for column vector
      HighFive::DataSpace user_bins_dataspace({hist.bins.size(), 1});
      HighFive::DataSet bins_ds = group.createDataSet("bins", user_bins_dataspace, HighFive::create_datatype<double>());
      bins_ds.write_raw(hist.bins.data()); // Use write_raw
      
      // Make it a dimension scale using C API
      if (H5DSset_scale(bins_ds.getId(), dim_label.c_str()) < 0) {
          std::cerr << "Warning: Failed to set dimension scale for " << group_name << std::endl;
      }

      // Add units to bins
      bins_ds.createAttribute<std::string>("units",
                                           HighFive::DataSpace::From(bin_unit))
          .write(bin_unit);
      
      // Add attribute for bin label (legacy but useful)
      group.createAttribute<std::string>(
               "bin_label", HighFive::DataSpace::From(dim_label))
          .write(dim_label);

      // Attach scales to all previously created datasets
      // We have to iterate again or keep track of them. 
      // Since we just created them, let's iterate over the group's children to attach scales
      for (const auto& ds_name : group.listObjectNames()) {
          if (ds_name == "bins") continue;
          
          HighFive::DataSet ds = group.getDataSet(ds_name);
          if (H5DSattach_scale(ds.getId(), bins_ds.getId(), 0) < 0) {
               std::cerr << "Warning: Failed to attach dimension scale for " << ds_name << std::endl;
          }
      }
    }
  } catch (const HighFive::Exception &err) {
    throw std::runtime_error("HDF5 Error: " + std::string(err.what()));
  }
}
