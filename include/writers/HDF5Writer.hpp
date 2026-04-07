/**
 * @file HDF5Writer.hpp
 * @brief Writer for HDF5 hierarchical data format.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <string>
#include <vector>

#include "BaseWriter.hpp"
#include "DistributionFunctions.hpp"

namespace Writer {

/**
 * @class HDF5Writer
 * @brief Manages the writing of distribution function data to HDF5 files.
 *
 * This class takes a DistributionFunctions object and provides methods to
 * write the calculated histograms to an HDF5 file.
 */
class HDF5Writer : public BaseWriter {
public:
  HDF5Writer() = default;

  std::string getName() const override { return "HDF5"; }
  std::vector<std::string> getExtensions() const override { return {".h5", ".hdf5"}; }

  void write(const std::string &base_path, const DistributionFunctions &df,
             bool smoothing) const override {
    writeHDF(base_path + ".h5", df);
  }

  /**
   * @brief Writes all available histograms to a single HDF5 file.
   * @param filename The full path of the HDF5 file to write.
   * @param df The DistributionFunctions object containing the data.
   */
  void writeHDF(const std::string &filename, const DistributionFunctions &df) const;
};

} // namespace Writer
