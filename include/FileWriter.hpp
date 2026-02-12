// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
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
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//
  explicit FileWriter(const DistributionFunctions &df);

  //-------------------------------------------------------------------------//
  //-------------------------------- Methods --------------------------------//
  //-------------------------------------------------------------------------//

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

  /**
   * @brief Writes all available histograms to a single HDF5 file.
   * @param filename The full path of the HDF5 file to write.
   */
  void writeHDF(const std::string &filename) const;

private:
  /**
   * @brief The core implementation for writing a single histogram
   * and it's smoothed histograms to a single CSV file.
   * @param filename The full path of the file to write.
   * @param hist The Histogram data structure to write.
   */
  void writeHistogramToCSV(const std::string &filename,
                           const Histogram &hist) const;

  const DistributionFunctions &df_;
};
