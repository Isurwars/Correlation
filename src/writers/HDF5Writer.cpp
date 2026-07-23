/**
 * @file HDF5Writer.cpp
 * @brief Implementation of the HDF5 writer.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "writers/HDF5Writer.hpp"
#include "writers/WriterFactory.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#include <hdf5_hl.h>
#include <highfive/highfive.hpp>

namespace correlation::writers {

namespace {

// Automatic registration
const bool registered = WriterFactory::registerTypeSafe<HDF5Writer>("HDF5Writer");

struct DatasetWriteQuery {
  std::string dim_label;
  std::string bin_unit;
  std::string data_unit;
  std::string description;
  std::string bin_ds_name;
};

void writeColumnDataset(HighFive::Group &group, size_t col, const std::string &col_name,
                        const std::vector<float> &col_data, const DatasetWriteQuery &query) {
  // Add prefix
  std::stringstream str_stream;
  str_stream << std::setw(2) << std::setfill('0') << col << "_" << col_name;
  std::string dataset_name = str_stream.str();

  HighFive::DataSet dataset = group.createDataSet(dataset_name, col_data);

  // Add attributes
  if (col == 0) {
    dataset.createAttribute<std::string>("Long Name", HighFive::DataSpace::From(query.dim_label))
        .write(query.dim_label);
    dataset.createAttribute<std::string>("Units", HighFive::DataSpace::From(query.bin_unit)).write(query.bin_unit);
    dataset.createAttribute<std::string>("Comments", HighFive::DataSpace::From(query.description))
        .write(query.description);

    // Make this a dimension scale
    herr_t status = H5DSset_scale(dataset.getId(), query.dim_label.c_str());
    if (status < 0) {
      std::cerr << "Warning: Failed to set dimension scale for " << dataset_name << '\n';
    }
  } else {
    dataset.createAttribute<std::string>("Long Name", HighFive::DataSpace::From(col_name)).write(col_name);
    dataset.createAttribute<std::string>("Units", HighFive::DataSpace::From(query.data_unit)).write(query.data_unit);
    dataset.createAttribute<std::string>("Comments", HighFive::DataSpace::From(col_name)).write(col_name);

    if (group.exist(query.bin_ds_name)) {
      HighFive::DataSet bin_ds = group.getDataSet(query.bin_ds_name);
      herr_t status = H5DSattach_scale(dataset.getId(), bin_ds.getId(), 0);
      if (status < 0) {
        std::cerr << "Warning: Failed to attach dimension scale for " << dataset_name << '\n';
      }
    }
  }
}

void writeHistogramToGroup(HighFive::File &file, const std::string &name,
                           const correlation::analysis::Histogram &hist) {
  // Sanitize group name: replace '(' with '_' and remove ')'
  std::string group_name = name;
  std::replace(group_name.begin(), group_name.end(), '(', '_');
  group_name.erase(std::remove(group_name.begin(), group_name.end(), ')'), group_name.end());

  // Also replace '/' with '_' just in case
  std::replace(group_name.begin(), group_name.end(), '/', '_');
  // Replace spaces with underscores
  std::replace(group_name.begin(), group_name.end(), ' ', '_');

  // Create a property list for group creation
  HighFive::GroupCreateProps props;
  props.add(HighFive::LinkCreationOrder(HighFive::CreationOrder::Tracked | HighFive::CreationOrder::Indexed));

  HighFive::Group group = file.createGroup(group_name, props);

  // Get metadata if available
  std::string bin_unit = hist.x_unit.empty() ? "arbitrary units" : hist.x_unit;
  std::string data_unit = hist.y_unit.empty() ? "arbitrary units" : hist.y_unit;
  std::string description = hist.description.empty() ? "" : hist.description;
  std::string dim_label = hist.x_label.empty() ? "x" : hist.x_label;

  // Add description attribute to the group
  if (!description.empty()) {
    group.createAttribute<std::string>("description", HighFive::DataSpace::From(description)).write(description);
  }

  // Prepare headers and data keys
  std::vector<std::string> headers;
  // Start with bins
  headers.push_back(dim_label);
  // Then raw data keys
  std::vector<std::string> raw_keys;
  raw_keys.reserve(hist.partials.size());
  for (const auto &[key, unsued] : hist.partials) {
    raw_keys.push_back(key);
  }
  std::ranges::sort(raw_keys);
  headers.insert(headers.end(), raw_keys.begin(), raw_keys.end());

  // Then smoothed data keys
  std::vector<std::string> smoothed_keys;
  if (!hist.smoothed_partials.empty()) {
    smoothed_keys.reserve(hist.smoothed_partials.size());
    for (const auto &[key, unsued] : hist.smoothed_partials) {
      smoothed_keys.push_back(key + "_smoothed");
    }
    std::ranges::sort(smoothed_keys);
    headers.insert(headers.end(), smoothed_keys.begin(), smoothed_keys.end());
  }

  // Determine dimensions
  size_t n_rows = hist.bins.size();
  size_t n_cols = headers.size();

  // Prepare the dimension scale name (column 0)
  std::string bin_col_name = headers[0];
  std::replace(bin_col_name.begin(), bin_col_name.end(), '(', '_');
  std::replace(bin_col_name.begin(), bin_col_name.end(), ')', '_');
  std::replace(bin_col_name.begin(), bin_col_name.end(), '/', '_');
  std::replace(bin_col_name.begin(), bin_col_name.end(), ' ', '_');

  std::stringstream ss0;
  ss0 << "00_" << bin_col_name;
  std::string bin_ds_name = ss0.str();

  DatasetWriteQuery query{dim_label, bin_unit, data_unit, description, bin_ds_name};

  for (size_t col = 0; col < n_cols; ++col) {
    std::string col_name = headers[col];
    // Sanitize column name
    std::replace(col_name.begin(), col_name.end(), '(', '_');
    std::replace(col_name.begin(), col_name.end(), ')', '_');
    std::replace(col_name.begin(), col_name.end(), '/', '_');
    std::replace(col_name.begin(), col_name.end(), ' ', '_');

    std::vector<float> col_data(n_rows);
    if (col == 0) {
      col_data.assign(hist.bins.begin(), hist.bins.end());
    } else if (col <= raw_keys.size()) {
      const auto &src = hist.partials.at(raw_keys[col - 1]);
      col_data.assign(src.begin(), src.end());
    } else {
      std::string original_key = smoothed_keys[col - 1 - raw_keys.size()];
      // remove "_smoothed" suffix for lookup
      original_key = original_key.substr(0, original_key.size() - 9);
      const auto &src = hist.smoothed_partials.at(original_key);
      col_data.assign(src.begin(), src.end());
    }

    writeColumnDataset(group, col, col_name, col_data, query);
  }
}

} // namespace

void HDF5Writer::writeHDF(const std::string &filename, const correlation::analysis::DistributionFunctions &dists) {
  try {
    HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

    // Iterate over all calculated histograms and write them as HDF5 groups
    for (const auto &[name, hist] : dists.getAllHistograms()) {
      if (hist.partials.empty() || hist.bins.empty()) {
        continue;
      }
      writeHistogramToGroup(file, name, hist);
    }
  } catch (const HighFive::Exception &err) {
    throw std::runtime_error("HDF5 Error: " + std::string(err.what()));
  }
}

} // namespace correlation::writers
