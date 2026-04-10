/**
 * @file BaseWriter.hpp
 * @brief Abstract base class for all file format writers.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"

#include <string>
#include <vector>

namespace correlation::writers {

/**
 * @brief Base class for all file writers.
 */
class BaseWriter {
public:
  virtual ~BaseWriter() = default;

  /**
   * @brief Returns the display name of the format.
   */
  virtual std::string getName() const = 0;

  /**
   * @brief Returns a list of supported file extensions.
   */
  virtual std::vector<std::string> getExtensions() const = 0;

  /**
   * @brief Writes the distribution function data to file(s).
   * @param base_path The base name for the output files.
   * @param df The DistributionFunctions object containing the data.
   * @param smoothing Whether to include smoothed data.
   */
  virtual void write(const std::string &base_path,
                     const correlation::analysis::DistributionFunctions &df,
                     bool smoothing) const = 0;
};

} // namespace correlation::writers
