/**
 * @file ArrowWriter.hpp
 * @brief Writer for Apache Arrow/Parquet columnar output format.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseWriter.hpp"
#include "analysis/DistributionFunctions.hpp"

#include <string>
#include <vector>

namespace correlation::writers {

/**
 * @class ArrowWriter
 * @brief Manages the writing of distribution function data to Parquet files via
 * Apache Arrow.
 *
 * This class takes a DistributionFunctions object and provides methods to
 * write the calculated histograms to Parquet files.
 */
class ArrowWriter : public BaseWriter {
public:
  ArrowWriter() = default;

  std::string getName() const override { return "Parquet"; }
  std::vector<std::string> getExtensions() const override { return {".parquet"}; }

  void write(const std::string &base_path, const correlation::analysis::DistributionFunctions &dists,
             bool smoothing) const override {
    writeAllParquet(base_path, dists, smoothing);
  }

  /**
   * @brief Writes all available histograms to Parquet files.
   *
   * For each histogram in the DistributionFunctions object (e.g., "g(r)"),
   * this method will create a corresponding Parquet file (e.g.,
   * "base_path_g.parquet").
   *
   * @param base_path The base name for the output files (e.g.,
   * "output/my_sample").
   * @param dists The DistributionFunctions object containing the data.
   * @param write_smoothed If true, also writes smoothed data columns.
   */
  static void writeAllParquet(const std::string &base_path, const correlation::analysis::DistributionFunctions &dists,
                              bool write_smoothed = false);

private:
  /**
   * @brief The core implementation for writing a single histogram
   * and its smoothed histograms to a single Parquet file.
   * @param filename The full path of the file to write.
   * @param hist The Histogram data structure to write.
   */
  static void writeHistogramToParquet(const std::string &filename, const correlation::analysis::Histogram &hist);
};

} // namespace correlation::writers
