/**
 * @file ITrajectoryReader.hpp
 * @brief Abstract interface for parsing multi-frame time-series trajectories.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Trajectory.hpp"

#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace correlation::readers {

/**
 * @brief Progress report callback function type.
 * Invoked with normalized progress in [0.0, 1.0] and a status message.
 */
using ProgressCallback = std::function<void(float, const std::string &)>;

/**
 * @brief Segregated contract for file format parsers that read multi-frame simulation trajectories.
 */
class ITrajectoryReader {
public:
  ITrajectoryReader() = default;
  virtual ~ITrajectoryReader() = default;

  ITrajectoryReader(const ITrajectoryReader &) = delete;
  ITrajectoryReader &operator=(const ITrajectoryReader &) = delete;
  ITrajectoryReader(ITrajectoryReader &&) = delete;
  ITrajectoryReader &operator=(ITrajectoryReader &&) = delete;

  /**
   * @brief Reads a trajectory from a file.
   * @param[in] filename Path to the file to be read.
   * @param[in] progress_callback Optional callback to report progress [0.0, 1.0].
   * @return A Trajectory object containing all loaded frames.
   */
  [[nodiscard]] virtual correlation::core::Trajectory
  readTrajectory(const std::string &filename, ProgressCallback progress_callback = nullptr) = 0;

  /**
   * @brief Reads a trajectory concatenated across multiple files.
   * @param[in] filenames List of paths to files to be read.
   * @param[in] progress_callback Optional callback to report progress [0.0, 1.0].
   * @return A Trajectory object containing all loaded frames.
   */
  [[nodiscard]] virtual correlation::core::Trajectory
  readTrajectory(const std::vector<std::string> &filenames,
                 ProgressCallback progress_callback = nullptr) {
    if (filenames.empty()) {
      throw std::runtime_error("No files provided for trajectory reading.");
    }
    if (filenames.size() == 1) {
      return readTrajectory(filenames.front(), progress_callback);
    }
    throw std::runtime_error("Multi-file trajectory reading is not supported for this format.");
  }
};

} // namespace correlation::readers
