/**
 * @file FileWriter.cpp
 * @brief Implementation of the unified file writer.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "writers/FileWriter.hpp"
#include "writers/WriterFactory.hpp"

namespace correlation::writers {

FileWriter::FileWriter(const correlation::analysis::DistributionFunctions &df)
    : df_(df) {}

void FileWriter::write(const std::string &base_path, bool use_csv,
                       bool use_hdf5, bool use_parquet, bool smoothing) const {
  auto &factory = WriterFactory::instance();

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

} // namespace correlation::writers
