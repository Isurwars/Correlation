/**
 * @file HDF5Writer.cpp
 * @brief Implementation of the simplified HDF5 writer.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "writers/HDF5Writer.hpp"
#include "writers/WriterFactory.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <highfive/highfive.hpp>

namespace correlation::writers {

namespace {

// Automatic registration
const bool registered = WriterFactory::registerTypeSafe<HDF5Writer>("HDF5Writer");

void writeHistogramToGroup(HighFive::File &file, const std::string &name,
                           const correlation::analysis::Histogram &hist) {
  // Sanitize group name
  std::string group_name = name;
  std::ranges::replace(group_name, '(', '_');
  std::ranges::replace(group_name, ')', '_');
  std::ranges::replace(group_name, '/', '_');
  std::ranges::replace(group_name, ' ', '_');

  HighFive::Group group = file.createGroup(group_name);

  // Metadata matching CSV header standard
  const std::string bin_unit = hist.x_unit.empty() ? "arbitrary units" : hist.x_unit;
  const std::string data_unit = hist.y_unit.empty() ? "arbitrary units" : hist.y_unit;
  const std::string description = hist.description.empty() ? "Data export" : hist.description;
  const std::string dim_label = hist.x_label.empty() ? "x" : hist.x_label;

  group.createAttribute<std::string>("description", HighFive::DataSpace::From(description)).write(description);

  // 1. Write Coordinate / Bin Dataset
  const std::vector<float> bin_data(hist.bins.begin(), hist.bins.end());
  HighFive::DataSet bin_dataset = group.createDataSet(dim_label, bin_data);
  bin_dataset.createAttribute<std::string>("units", HighFive::DataSpace::From(bin_unit)).write(bin_unit);
  bin_dataset.createAttribute<std::string>("label", HighFive::DataSpace::From(dim_label)).write(dim_label);
  bin_dataset.createAttribute<std::string>("description", HighFive::DataSpace::From(description)).write(description);

  // 2. Write Raw Partial Datasets
  std::vector<std::string> raw_keys;
  raw_keys.reserve(hist.partials.size());
  for (const auto &[key, value] : hist.partials) {
    raw_keys.push_back(key);
  }
  std::ranges::sort(raw_keys);

  for (const auto &key : raw_keys) {
    const auto &values = hist.partials.at(key);
    const std::vector<float> float_values(values.begin(), values.end());

    HighFive::DataSet dataset = group.createDataSet(key, float_values);
    dataset.createAttribute<std::string>("units", HighFive::DataSpace::From(data_unit)).write(data_unit);
    dataset.createAttribute<std::string>("label", HighFive::DataSpace::From(key)).write(key);
  }

  // 3. Write Smoothed Partial Datasets
  std::vector<std::string> smoothed_keys;
  smoothed_keys.reserve(hist.smoothed_partials.size());
  for (const auto &[key, value] : hist.smoothed_partials) {
    smoothed_keys.push_back(key);
  }
  std::ranges::sort(smoothed_keys);

  for (const auto &key : smoothed_keys) {
    const auto &values = hist.smoothed_partials.at(key);
    const std::vector<float> float_values(values.begin(), values.end());
    const std::string smoothed_name = key + "_smoothed";

    HighFive::DataSet dataset = group.createDataSet(smoothed_name, float_values);
    dataset.createAttribute<std::string>("units", HighFive::DataSpace::From(data_unit)).write(data_unit);
    dataset.createAttribute<std::string>("label", HighFive::DataSpace::From(smoothed_name)).write(smoothed_name);
  }
}

} // namespace

void HDF5Writer::writeHDF(const std::string &filename, const correlation::analysis::DistributionFunctions &dists) {
  try {
    HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

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
