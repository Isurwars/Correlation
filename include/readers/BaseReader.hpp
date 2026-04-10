/**
 * @file BaseReader.hpp
 * @brief Abstract base class for all file format readers.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <functional>
#include <string>
#include <vector>

namespace correlation::readers {

/**
 * @brief Base class for all file readers.
 */
class BaseReader {
public:
  virtual ~BaseReader() = default;

  /**
   * @brief Returns the display name of the format.
   * @return A string representing the format name (e.g. "LAMMPS Dump", "CIF").
   */
  virtual std::string getName() const = 0;

  /**
   * @brief Returns a list of supported file extensions.
   * @return A vector of strings containing extensions (e.g. {".cif", ".ext"}).
   */
  virtual std::vector<std::string> getExtensions() const = 0;

  /**
   * @brief Returns true if this reader supports trajectory files.
   * @return True if the format can store multiple time steps.
   */
  virtual bool isTrajectory() const = 0;

  /**
   * @brief Reads an atomic structure from a file.
   * @param filename Path to the file to be read.
   * @param progress_callback Optional callback to report progress [0.0, 1.0].
   * @return A Cell object containing the loaded structure.
   */
  virtual correlation::core::Cell readStructure(
      const std::string &filename,
      std::function<void(float, const std::string &)> progress_callback =
          nullptr) = 0;

  /**
   * @brief Reads a trajectory from a file.
   * @param filename Path to the file to be read.
   * @param progress_callback Optional callback to report progress [0.0, 1.0].
   * @return A Trajectory object containing all loaded frames.
   */
  virtual correlation::core::Trajectory readTrajectory(
      const std::string &filename,
      std::function<void(float, const std::string &)> progress_callback =
          nullptr) = 0;
};

} // namespace correlation::readers
