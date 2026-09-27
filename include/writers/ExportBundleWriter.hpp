/**
 * @file ExportBundleWriter.hpp
 * @brief Consolidated ZIP export bundle writer packaging SVGs, CSVs, HDF5, Parquet, and summary.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "plotters/PlotTypes.hpp"

#include <expected>
#include <filesystem>
#include <string>

namespace correlation::writers {

/**
 * @struct ExportBundleOptions
 * @brief Configuration parameters for consolidated export bundle generation.
 */
struct ExportBundleOptions {
  bool include_svg{true};                          ///< Package SVG plots into `plots/`
  bool include_csv{true};                          ///< Package CSV tables into `csv/`
  bool include_hdf5{false};                        ///< Package HDF5 archive into `hdf5/`
  bool include_parquet{false};                     ///< Package Parquet tables into `parquet/`
  bool include_summary{true};                      ///< Package summary text file
  bool smoothing{true};                            ///< Include smoothed data in CSV / Parquet
  correlation::plotters::PlotConfig plot_config{}; ///< Configuration for rendered SVGs
};

/**
 * @class ExportBundleWriter
 * @brief Packages scientific analysis outputs into a structured ZIP archive.
 */
class ExportBundleWriter {
public:
  /**
   * @brief Writes a consolidated ZIP bundle from DistributionFunctions.
   * @param zip_path Destination path for the `.zip` archive.
   * @param dists The DistributionFunctions containing computed analysis data.
   * @param options Bundle configuration options.
   * @return std::expected<void, std::string> Success or descriptive error message.
   */
  static std::expected<void, std::string>
  writeBundle(const std::filesystem::path &zip_path,
              const correlation::analysis::DistributionFunctions &dists,
              const ExportBundleOptions &options = {});
};

} // namespace correlation::writers
