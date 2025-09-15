// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <string>

#include "DistributionFunctions.hpp"

/**
 * @class FileWriter
 * @brief Manages the writing of distribution function data to CSV files.
 *
 * This class takes a DistributionFunctions object and provides methods to
 * write the calculated histograms (both raw and smoothed) to disk in a
 * well-formatted, data-driven way.
 */
class FileWriter {
public:
  /**
   * @brief Constructs a FileWriter linked to a DistributionFunctions object.
   * @param df The DistributionFunctions object containing the data to be
   * written.
   */
  explicit FileWriter(const DistributionFunctions &df);

  /**
   * @brief Writes all available histograms to appropriately named CSV files.
   *
   * For each histogram in the DistributionFunctions object (e.g., "g(r)"),
   * this method will create a corresponding CSV file (e.g., "base_path_g.csv").
   *
   * @param base_path The base name for the output files (e.g.,
   * "output/my_sample").
   * @param write_smoothed If true, also writes smoothed data to separate files
   * (e.g., "base_path_g_smoothed.csv").
   */
  void writeAllCSVs(const std::string &base_path,
                    bool write_smoothed = false) const;

private:
  /**
   * @brief The core implementation for writing a single histogram to a CSV
   * file.
   * @param filename The full path of the file to write.
   * @param hist The Histogram data structure to write.
   * @param use_smoothed_data If true, writes the smoothed partials; otherwise,
   * writes the raw partials.
   */
  void writeHistogramToCSV(const std::string &filename, const Histogram &hist,
                           bool use_smoothed_data) const;

  const DistributionFunctions &df_;
};
