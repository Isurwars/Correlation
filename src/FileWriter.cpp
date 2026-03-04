// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "FileWriter.hpp"
#include "writers/CSVWriter.hpp"
#include "writers/HDF5Writer.hpp"
// #include "writers/ParquetWriter.hpp"

FileWriter::FileWriter(const DistributionFunctions &df) : df_(df) {}

void FileWriter::write(const std::string &base_path, bool use_csv,
                       bool use_hdf5, bool use_parquet, bool smoothing) const {
  if (use_csv) {
    CSVWriter csv_writer(df_);
    csv_writer.writeAllCSVs(base_path, smoothing);
  }

  if (use_hdf5) {
    HDF5Writer hdf5_writer(df_);
    hdf5_writer.writeHDF(base_path + ".h5");
  }

  if (use_parquet) {
    // ParquetWriter parquet_writer(df_);
    // parquet_writer.writeAll(base_path, smoothing);
  }
}
