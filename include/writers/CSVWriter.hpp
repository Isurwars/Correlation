/**
 * @file CSVWriter.hpp
 * @brief Writer for CSV tabular output format.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <string>
#include <vector>

#include "BaseWriter.hpp"
#include "DistributionFunctions.hpp"

namespace correlation::writers {

/**
 * @class CSVWriter
 * @brief Manages the writing of distribution function data to CSV files.
 *
 * This class takes a DistributionFunctions object and provides methods to
 * write the calculated histograms (both raw and smoothed) to disk in a
 * well-formatted, data-driven way.
 */
class CSVWriter : public BaseWriter {
public:
  CSVWriter() = default;

  std::string getName() const override { return "CSV"; }
  std::vector<std::string> getExtensions() const override { return {".csv"}; }

  void write(const std::string &base_path, const DistributionFunctions &df,
             bool smoothing) const override {
    writeAllCSVs(base_path, df, smoothing);
  }

  /**
   * @brief Writes all available histograms to appropriately named CSV files.
   *
   * For each histogram in the DistributionFunctions object (e.g., "g(r)"),
   * this method will create a corresponding CSV file (e.g., "base_path_g.csv").
   *
   * @param base_path The base name for the output files (e.g.,
   * "output/my_sample").
   * @param df The DistributionFunctions object containing the data.
   * @param write_smoothed If true, also writes smoothed data to separate files
   * (e.g., "base_path_g_smoothed.csv").
   */
  void writeAllCSVs(const std::string &base_path, const DistributionFunctions &df,
                    bool write_smoothed = false) const;

private:
  /**
   * @brief The core implementation for writing a single histogram
   * and it's smoothed histograms to a single CSV file.
   * @param filename The full path of the file to write.
   * @param name The name of the histogram.
   * @param hist The Histogram data structure to write.
   */
  void writeHistogramToCSV(const std::string &filename, const std::string &name,
                           const Histogram &hist) const;
};

} // namespace correlation::writers
