/**
 * @file FileWriter.hpp
 * @brief Unified file writing interface for analysis output.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"

#include <expected>
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
   * @param dists The DistributionFunctions object containing the data to be
   * written.
   */
  explicit FileWriter(const correlation::analysis::DistributionFunctions &dists);

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
  void write(const std::string &base_path, bool use_csv, bool use_hdf5, bool use_parquet,
             bool smoothing) const;

  /**
   * @brief Writes a consolidated ZIP bundle containing plots, tables, and summary.
   *
   * @param zip_path The target path for the zip file (e.g., "output/bundle.zip").
   * @param use_csv Whether to include CSV tables.
   * @param use_hdf5 Whether to include HDF5 data.
   * @param use_parquet Whether to include Parquet data.
   * @param include_svg Whether to include SVG plots.
   * @param smoothing Whether to include smoothed curves.
   * @return std::expected<void, std::string> Success or error message.
   */
  std::expected<void, std::string>
  writeBundle(const std::string &zip_path, bool use_csv = true, bool use_hdf5 = false,
              bool use_parquet = false, bool include_svg = true, bool smoothing = true,
              const std::vector<std::string> &selected_algorithms = {}) const;

  /**
   * @brief Writes analysis outputs into categorized subfolders within a destination folder.
   *
   * @param dest_dir Destination directory path.
   * @param use_csv Whether to include CSV tables in `csv/`.
   * @param use_hdf5 Whether to include HDF5 data in `hdf5/`.
   * @param use_parquet Whether to include Parquet data in `arrow/`.
   * @param include_svg Whether to include SVG plots in `plots/`.
   * @param smoothing Whether to include smoothed curves.
   * @param selected_algorithms Optional filter list of algorithms to export (empty for all).
   * @return std::expected<void, std::string> Success or error message.
   */
  std::expected<void, std::string>
  writeFolder(const std::string &dest_dir, bool use_csv = true, bool use_hdf5 = false,
              bool use_parquet = false, bool include_svg = true, bool smoothing = true,
              const std::vector<std::string> &selected_algorithms = {}) const;

private:
  void writeSummaryFile(const std::string &base_path) const;

  /** @brief Pointer to the analysis data structure. */
  const correlation::analysis::DistributionFunctions *df_ = nullptr;
};

} // namespace correlation::writers
