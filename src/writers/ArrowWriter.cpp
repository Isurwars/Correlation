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
// NOLINTNEXTLINE(cert-err58-cpp, bugprone-throwing-static-initialization)
const bool registered = WriterFactory::instance().registerWriter(std::make_unique<ArrowWriter>());

std::vector<std::string> getSortedKeys(const std::map<std::string, std::vector<double>> &map) {
  std::vector<std::string> keys;
  keys.reserve(map.size());
  for (const auto &elem : map) {
    keys.push_back(elem.first);
  }
  std::ranges::sort(keys);
  return keys;
}

void addDoubleColumn(const std::string &name, const std::vector<double> &values, arrow::FieldVector &fields,
                     std::vector<std::shared_ptr<arrow::Array>> &arrays) {
  fields.push_back(arrow::field(name, arrow::float64()));
  arrow::DoubleBuilder builder;
  PARQUET_THROW_NOT_OK(builder.AppendValues(values));
  std::shared_ptr<arrow::Array> array;
  PARQUET_THROW_NOT_OK(builder.Finish(&array));
  arrays.push_back(array);
}

void writeTableToParquet(const std::string &filename, const std::shared_ptr<arrow::Table> &table) {
  std::shared_ptr<arrow::io::FileOutputStream> outfile;
  PARQUET_ASSIGN_OR_THROW(outfile, arrow::io::FileOutputStream::Open(filename));

  // Writer properties
  std::shared_ptr<parquet::WriterProperties> props =
      parquet::WriterProperties::Builder().compression(parquet::Compression::UNCOMPRESSED)->build();

  PARQUET_THROW_NOT_OK(
      parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 1024LL * 1024LL, props));
}

} // namespace

void ArrowWriter::writeAllParquet(const std::string &base_path,
                                  const correlation::analysis::DistributionFunctions &dists, bool /*write_smoothed*/) {
  for (const auto &[name, hist] : dists.getAllHistograms()) {
    try {
      if (hist.partials.empty() || hist.bins.empty() || hist.file_suffix.empty()) {
        continue;
      }

      std::string filename = base_path + hist.file_suffix + ".parquet";
      writeHistogramToParquet(filename, hist);

    } catch (const std::exception &e) {
      std::cerr << "Error writing Parquet file for '" << name << "': " << e.what() << "\n";
    }
  }
}

void ArrowWriter::writeHistogramToParquet(const std::string &filename, const correlation::analysis::Histogram &hist) {
  if (hist.partials.empty() || hist.bins.empty()) {
    return;
  }

  // Get sorted keys for both raw and smoothed data
  std::vector<std::string> raw_keys = getSortedKeys(hist.partials);
  std::vector<std::string> smoothed_keys = getSortedKeys(hist.smoothed_partials);

  // Get metadata if available
  std::string dim_label = hist.x_label.empty() ? "x" : hist.x_label;

  // --- Build Arrow Schema and Data Arrays ---
  arrow::FieldVector fields;
  std::vector<std::shared_ptr<arrow::Array>> arrays;

  // 1. Bin Column
  addDoubleColumn(dim_label, hist.bins, fields, arrays);

  // 2. Raw Data Columns
  for (const auto &key : raw_keys) {
    addDoubleColumn(key, hist.partials.at(key), fields, arrays);
  }

  // 3. Smoothed Data Columns
  for (const auto &key : smoothed_keys) {
    addDoubleColumn(key + "_smoothed", hist.smoothed_partials.at(key), fields, arrays);
  }

  auto schema = arrow::schema(fields);
  int64_t num_rows = static_cast<int64_t>(hist.bins.size());
  auto table = arrow::Table::Make(schema, arrays, num_rows);

  // --- Write to Parquet File ---
  writeTableToParquet(filename, table);
}

} // namespace correlation::writers
