/**
 * @file BaseReader.hpp
 * @brief Abstract base class for all file format readers.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <string>
#include <vector>
#include <functional>
#include "Cell.hpp"
#include "Trajectory.hpp"

namespace FileReader {

/**
 * @brief Base class for all file readers.
 */
class BaseReader {
public:
  virtual ~BaseReader() = default;

  /**
   * @brief Returns the display name of the format.
   */
  virtual std::string getName() const = 0;

  /**
   * @brief Returns a list of supported file extensions.
   */
  virtual std::vector<std::string> getExtensions() const = 0;

  /**
   * @brief Returns true if this reader supports trajectory files.
   */
  virtual bool isTrajectory() const = 0;

  /**
   * @brief Reads an atomic structure from a file.
   */
  virtual Cell readStructure(const std::string &filename,
                             std::function<void(float, const std::string &)>
                                 progress_callback = nullptr) = 0;

  /**
   * @brief Reads a trajectory from a file.
   */
  virtual Trajectory readTrajectory(const std::string &filename,
                                    std::function<void(float, const std::string &)>
                                        progress_callback = nullptr) = 0;
};

} // namespace FileReader
