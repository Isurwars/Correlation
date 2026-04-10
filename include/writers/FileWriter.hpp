/**
 * @file FileWriter.hpp
 * @brief Unified file writing interface for analysis output.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"

#include <string>

namespace correlation::writers {

/**
 * @class FileWriter
 * @brief Facade class that manages writing data to various file formats.
 *
 * This class orchestrates the usage of specific writers (CSV, HDF5, Parquet)
 * based on provided options.
 */
class FileWriter {
public:
  /**
   * @brief Constructs a FileWriter linked to a DistributionFunctions object.
   * @param df The DistributionFunctions object containing the data to be
   * written.
   */
  explicit FileWriter(const correlation::analysis::DistributionFunctions &df);

  /**
   * @brief Writes the available histograms using specified formats.
   *
   * @param base_path The base name for the output files (e.g.,
   * "output/my_sample").
   * @param use_csv Whether to write CSV files.
   * @param use_hdf5 Whether to write an HDF5 file.
   * @param use_parquet Whether to write Parquet files.
   * @param smoothing Whether to include smoothed data.
   */
  void write(const std::string &base_path, bool use_csv, bool use_hdf5,
             bool use_parquet, bool smoothing) const;

private:
  /** @brief Reference to the analysis data structure. */
  const correlation::analysis::DistributionFunctions &df_;
};

} // namespace correlation::writers
