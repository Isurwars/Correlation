/**
 * @file ArrowWriter.cpp
 * @brief Implementation of the Apache Arrow/Parquet writer.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "writers/ArrowWriter.hpp"
#include "writers/WriterFactory.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#include <parquet/exception.h>

namespace correlation::writers {

namespace {

// Automatic registration
const bool registered = WriterFactory::registerTypeSafe<ArrowWriter>("ArrowWriter");

std::vector<std::string> getSortedKeys(const std::map<std::string, std::vector<real_t>> &map) {
  std::vector<std::string> keys;
  keys.reserve(map.size());
  for (const auto &elem : map) {
    keys.push_back(elem.first);
  }
  std::ranges::sort(keys);
  return keys;
}

void addFloatColumn(const std::string &name, const std::vector<real_t> &values, arrow::FieldVector &fields,
                    std::vector<std::shared_ptr<arrow::Array>> &arrays) {
  fields.push_back(arrow::field(name, arrow::float32()));
  arrow::FloatBuilder builder;
  std::vector<float> float_values(values.begin(), values.end());
  PARQUET_THROW_NOT_OK(builder.AppendValues(float_values.data(), static_cast<int64_t>(float_values.size())));
  std::shared_ptr<arrow::Array> array;
  PARQUET_THROW_NOT_OK(builder.Finish(&array));
  arrays.push_back(array);
}

void writeTableToParquet(const std::string &filename, const std::shared_ptr<arrow::Table> &table) {
  std::shared_ptr<arrow::io::FileOutputStream> outfile;
  PARQUET_ASSIGN_OR_THROW(outfile, arrow::io::FileOutputStream::Open(filename));

  std::shared_ptr<parquet::WriterProperties> props = parquet::WriterProperties::Builder().build();

  PARQUET_THROW_NOT_OK(
      parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 1024LL * 1024LL, props));
  PARQUET_THROW_NOT_OK(outfile->Close());
}

} // namespace

void ArrowWriter::writeAllParquet(const std::string &base_path,
                                  const correlation::analysis::DistributionFunctions &dists, bool /*write_smoothed*/) {
  const auto &all_histograms = dists.getAllHistograms();
  for (const auto &[name, hist] : all_histograms) {
    try {
      if (hist.partials.empty() || hist.bins.empty() || hist.file_suffix.empty()) {
        continue;
      }

      std::string filename = base_path + hist.file_suffix + ".parquet";

      // For normalized PAD/DAD, inject raw companion partials.
      const correlation::analysis::Histogram *raw_companion = nullptr;
      const std::string raw_key = name + "_raw";
      if ((name == "PAD" || name == "DAD") && all_histograms.contains(raw_key)) {
        raw_companion = &all_histograms.at(raw_key);
      }

      writeHistogramToParquet(filename, hist, raw_companion);

    } catch (const std::exception &e) {
      std::cerr << "Error writing Parquet file for '" << name << "': " << e.what() << "\n";
    }
  }
}

void ArrowWriter::writeHistogramToParquet(const std::string &filename, const correlation::analysis::Histogram &hist,
                                          const correlation::analysis::Histogram *raw_companion) {
  if (hist.partials.empty() || hist.bins.empty()) {
    return;
  }

  // Extract metadata matching CSV header standard
  std::string bin_unit = hist.x_unit.empty() ? "arbitrary units" : hist.x_unit;
  std::string data_unit = hist.y_unit.empty() ? "arbitrary units" : hist.y_unit;
  std::string description = hist.description.empty() ? "Data export" : hist.description;
  std::string dim_label = hist.x_label.empty() ? "x" : hist.x_label;

  // Get sorted keys for raw companion, normalized, and smoothed data
  std::vector<std::string> companion_keys;
  if (raw_companion != nullptr) {
    companion_keys.reserve(raw_companion->partials.size());
    for (const auto &elem : raw_companion->partials) {
      companion_keys.push_back(elem.first);
    }
    std::ranges::sort(companion_keys);
  }

  std::vector<std::string> raw_keys = getSortedKeys(hist.partials);
  std::vector<std::string> smoothed_keys = getSortedKeys(hist.smoothed_partials);

  // Build Arrow Schema and Data Columns
  arrow::FieldVector fields;
  std::vector<std::shared_ptr<arrow::Array>> arrays;

  // 1. Coordinate / Bin Column
  addFloatColumn(dim_label, hist.bins, fields, arrays);

  // 2. Raw Companion Columns (non-normalized counts)
  for (const auto &key : companion_keys) {
    addFloatColumn(key + "_raw", raw_companion->partials.at(key), fields, arrays);
  }

  // 3. Normalized Data Columns
  for (const auto &key : raw_keys) {
    addFloatColumn(key, hist.partials.at(key), fields, arrays);
  }

  // 4. Smoothed Data Columns
  for (const auto &key : smoothed_keys) {
    addFloatColumn(key + "_smoothed", hist.smoothed_partials.at(key), fields, arrays);
  }

  // Attach CSV-aligned header metadata to schema
  auto metadata = arrow::key_value_metadata(
      {"dim_label", "bin_unit", "data_unit", "description"},
      {dim_label, bin_unit, data_unit, description});
  auto schema = arrow::schema(fields, metadata);

  int64_t num_rows = static_cast<int64_t>(hist.bins.size());
  std::vector<std::shared_ptr<arrow::ChunkedArray>> chunked_arrays;
  chunked_arrays.reserve(arrays.size());
  for (const auto &arr : arrays) {
    chunked_arrays.push_back(std::make_shared<arrow::ChunkedArray>(arr));
  }
  auto table = arrow::Table::Make(schema, chunked_arrays, num_rows);

  // Write to Parquet File using default properties
  writeTableToParquet(filename, table);
}

} // namespace correlation::writers
