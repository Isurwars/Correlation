// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <string>

#include "DistributionFunctions.hpp"

/**
 * @class HDF5Writer
 * @brief Manages the writing of distribution function data to HDF5 files.
 *
 * This class takes a DistributionFunctions object and provides methods to
 * write the calculated histograms to an HDF5 file.
 */
class HDF5Writer {
public:
  /**
   * @brief Constructs a HDF5Writer linked to a DistributionFunctions object.
   * @param df The DistributionFunctions object containing the data to be
   * written.
   */
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//
  explicit HDF5Writer(const DistributionFunctions &df);

  //-------------------------------------------------------------------------//
  //-------------------------------- Methods --------------------------------//
  //-------------------------------------------------------------------------//

  /**
   * @brief Writes all available histograms to a single HDF5 file.
   * @param filename The full path of the HDF5 file to write.
   */
  void writeHDF(const std::string &filename) const;

private:
  const DistributionFunctions &df_;
};
