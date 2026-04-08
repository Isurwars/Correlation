/**
 * @file ArrowWriter.cpp
 * @brief Implementation of the Apache Arrow/Parquet writer.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "writers/ArrowWriter.hpp"
#include "writers/WriterFactory.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#include <parquet/exception.h>

namespace Writer {

// Automatic registration
static bool registered =
    WriterFactory::instance().registerWriter(std::make_unique<ArrowWriter>());

void ArrowWriter::writeAllParquet(const std::string &base_path,
                                  const DistributionFunctions &df,
                                  bool /*write_smoothed*/) const {
  const std::map<std::string, std::string> file_map = {
      {"g_r", "_g.parquet"},
      {"J_r", "_J.parquet"},
      {"G_r", "_G_reduced.parquet"},
      {"BAD", "_PAD.parquet"},
      {"DAD", "_DAD.parquet"},
      {"RD", "_RD.parquet"},
      {"S_q", "_S.parquet"},
      {"XRD", "_XRD.parquet"},
      {"CN", "_CN.parquet"},
      {"VACF", "_VACF.parquet"},
      {"Normalized VACF", "_VACF_norm.parquet"},
      {"VDOS", "_VDOS.parquet"},
      {"MSD", "_MSD.parquet"},
      {"D_eff", "_D_eff.parquet"}};

  for (const auto &[name, suffix] : file_map) {
    try {
      const auto &hist = df.getHistogram(name);

      if (hist.partials.empty() || hist.bins.empty() || hist.bins.size() == 0) {
        continue;
      }

      std::string filename = base_path + suffix;
      writeHistogramToParquet(filename, name, hist);

    } catch (const std::out_of_range &) {
      // This is not an error; it just means the histogram wasn't calculated.
      // We can silently skip it.
    } catch (const std::exception &e) {
      std::cerr << "Error writing Parquet file for '" << name
                << "': " << e.what() << std::endl;
    }
  }
}

void ArrowWriter::writeHistogramToParquet(const std::string &filename,
                                          const std::string &name,
                                          const Histogram &hist) const {
  if (hist.partials.empty() || hist.bins.empty()) {
    return;
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
  std::string dim_label = hist.x_label.empty() ? "x" : hist.x_label;

  // --- Build Arrow Schema and Data Arrays ---
  arrow::FieldVector fields;
  std::vector<std::shared_ptr<arrow::Array>> arrays;

  // 1. Bin Column
  fields.push_back(arrow::field(dim_label, arrow::float64()));
  arrow::DoubleBuilder bin_builder;
  PARQUET_THROW_NOT_OK(bin_builder.AppendValues(hist.bins));
  std::shared_ptr<arrow::Array> bin_array;
  PARQUET_THROW_NOT_OK(bin_builder.Finish(&bin_array));
  arrays.push_back(bin_array);

  // 2. Raw Data Columns
  for (const auto &key : raw_keys) {
    fields.push_back(arrow::field(key, arrow::float64()));
    arrow::DoubleBuilder data_builder;
    PARQUET_THROW_NOT_OK(data_builder.AppendValues(hist.partials.at(key)));
    std::shared_ptr<arrow::Array> data_array;
    PARQUET_THROW_NOT_OK(data_builder.Finish(&data_array));
    arrays.push_back(data_array);
  }

  // 3. Smoothed Data Columns
  for (const auto &key : smoothed_keys) {
    std::string col_name = key + "_smoothed";
    fields.push_back(arrow::field(col_name, arrow::float64()));
    arrow::DoubleBuilder smoothed_builder;
    PARQUET_THROW_NOT_OK(
        smoothed_builder.AppendValues(hist.smoothed_partials.at(key)));
    std::shared_ptr<arrow::Array> smoothed_array;
    PARQUET_THROW_NOT_OK(smoothed_builder.Finish(&smoothed_array));
    arrays.push_back(smoothed_array);
  }

  auto schema = arrow::schema(fields);
  int64_t num_rows = hist.bins.size();
  auto table = arrow::Table::Make(schema, arrays, num_rows);

  // --- Write to Parquet File ---
  std::shared_ptr<arrow::io::FileOutputStream> outfile;
  PARQUET_ASSIGN_OR_THROW(outfile, arrow::io::FileOutputStream::Open(filename));

  // Writer properties
  std::shared_ptr<parquet::WriterProperties> props =
      parquet::WriterProperties::Builder()
          .compression(parquet::Compression::UNCOMPRESSED)
          ->build();

  PARQUET_THROW_NOT_OK(parquet::arrow::WriteTable(
      *table, arrow::default_memory_pool(), outfile, 1024 * 1024, props));
}

} // namespace Writer
