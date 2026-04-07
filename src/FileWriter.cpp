/**
 * @file FileWriter.cpp
 * @brief Implementation of the unified file writer.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "FileWriter.hpp"
#include "writers/WriterFactory.hpp"

FileWriter::FileWriter(const DistributionFunctions &df) : df_(df) {}

void FileWriter::write(const std::string &base_path, bool use_csv,
                       bool use_hdf5, bool use_parquet, bool smoothing) const {
  auto &factory = Writer::WriterFactory::instance();

  if (use_csv) {
    if (auto writer = factory.getWriter("CSV")) {
      writer->write(base_path, df_, smoothing);
    }
  }

  if (use_hdf5) {
    if (auto writer = factory.getWriter("HDF5")) {
      writer->write(base_path, df_, smoothing);
    }
  }

  if (use_parquet) {
    if (auto writer = factory.getWriter("Parquet")) {
      writer->write(base_path, df_, smoothing);
    }
  }
}
