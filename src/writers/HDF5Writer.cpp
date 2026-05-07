/**
 * @file HDF5Writer.cpp
 * @brief Implementation of the HDF5 writer.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "writers/HDF5Writer.hpp"
#include "writers/WriterFactory.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <hdf5_hl.h>
#include <highfive/highfive.hpp>

namespace correlation::writers {

// Automatic registration
static bool registered =
    WriterFactory::instance().registerWriter(std::make_unique<HDF5Writer>());

void HDF5Writer::writeHDF(
    const std::string &filename,
    const correlation::analysis::DistributionFunctions &df) const {
  try {
    HighFive::File file(filename, HighFive::File::ReadWrite |
                                      HighFive::File::Create |
                                      HighFive::File::Truncate);

    // Iterate over all calculated histograms and write them as HDF5 groups
    for (const auto &[name, hist] : df.getAllHistograms()) {
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
      std::string bin_unit =
          hist.x_unit.empty() ? "arbitrary units" : hist.x_unit;
      std::string data_unit =
          hist.y_unit.empty() ? "arbitrary units" : hist.y_unit;
      std::string description =
          hist.description.empty() ? "" : hist.description;
      std::string dim_label = hist.x_label.empty() ? "x" : hist.x_label;

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
      headers.push_back(dim_label);
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
          ds.createAttribute<std::string>(
                "Comments", HighFive::DataSpace::From(description))
              .write(description);

          // Make this a dimension scale
          herr_t status = H5DSset_scale(ds.getId(), dim_label.c_str());
          if (status < 0) {
            std::cerr << "Warning: Failed to set dimension scale for "
                      << dataset_name << std::endl;
          }
        } else {
          std::string current_long_name = headers[col];
          std::string current_data_unit = data_unit;
          std::string current_comment = headers[col];

          ds.createAttribute<std::string>(
                "Long Name", HighFive::DataSpace::From(current_long_name))
              .write(current_long_name);
          ds.createAttribute<std::string>(
                "Units", HighFive::DataSpace::From(current_data_unit))
              .write(current_data_unit);
          ds.createAttribute<std::string>(
                "Comments", HighFive::DataSpace::From(current_comment))
              .write(current_comment);
          std::stringstream ss_bin;
          ss_bin << std::setw(2) << std::setfill('0') << 0 << "_" << headers[0];

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

} // namespace correlation::writers
