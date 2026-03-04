// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <string>

#include "DistributionFunctions.hpp"

/**
 * @class ArrowWriter
 * @brief Manages the writing of distribution function data to Parquet files via Apache Arrow.
 *
 * This class takes a DistributionFunctions object and provides methods to
 * write the calculated histograms to Parquet files.
 */
class ArrowWriter {
public:
  /**
   * @brief Constructs a ArrowWriter linked to a DistributionFunctions object.
   * @param df The DistributionFunctions object containing the data to be
   * written.
   */
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//
  explicit ArrowWriter(const DistributionFunctions &df);

  //-------------------------------------------------------------------------//
  //-------------------------------- Methods --------------------------------//
  //-------------------------------------------------------------------------//

  /**
   * @brief Writes all available histograms to Parquet files.
   *
   * For each histogram in the DistributionFunctions object (e.g., "g(r)"),
   * this method will create a corresponding Parquet file (e.g., "base_path_g.parquet").
   *
   * @param base_path The base name for the output files (e.g.,
   * "output/my_sample").
   * @param write_smoothed If true, also writes smoothed data columns.
   */
  void writeAllParquet(const std::string &base_path,
                       bool write_smoothed = false) const;

private:
  //-------------------------------------------------------------------------//
  //------------------------------ Helpers ----------------------------------//
  //-------------------------------------------------------------------------//

  /**
   * @brief The core implementation for writing a single histogram
   * and its smoothed histograms to a single Parquet file.
   * @param filename The full path of the file to write.
   * @param name The name of the histogram.
   * @param hist The Histogram data structure to write.
   */
  void writeHistogramToParquet(const std::string &filename, const std::string &name,
                               const Histogram &hist) const;

  const DistributionFunctions &df_;
};
