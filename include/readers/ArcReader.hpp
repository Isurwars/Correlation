/**
 * @file ArcReader.hpp
 * @brief Reader for BIOVIA/Accelrys ARC trajectory files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#pragma once

#include "BaseReader.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <functional>
#include <string>
#include <vector>

namespace correlation::readers {

/**
 * @brief Reads an ARC file and returns a vector of correlation::core::Cell
 * objects.
 */
class ArcReader : public BaseReader {
public:
  std::string getName() const override { return "Accelrys ARC"; }
  std::vector<std::string> getExtensions() const override { return {"arc"}; }
  bool isTrajectory() const override { return true; }

  correlation::core::Cell
  readStructure(const std::string &filename,
                std::function<void(float, const std::string &)> progress_callback = nullptr) override;

  correlation::core::Trajectory
  readTrajectory(const std::string &filename,
                 std::function<void(float, const std::string &)> progress_callback = nullptr) override;

  /**
   * @brief Low-level parser for the specific file format.
   * @param file_name Path to the source file.
   * @param progress_callback Optional callback for progress updates.
   * @return A vector of parsed frames.
   */
  static std::vector<correlation::core::Cell>
  read(const std::string &file_name,
       const std::function<void(float, const std::string &)> &progress_callback = nullptr);

private:
  static void updateProgress(std::streampos current_pos, std::streampos file_size, std::streampos &last_progress_pos,
                             size_t update_interval,
                             const std::function<void(float, const std::string &)> &progress_callback);

  static void parseLine(const std::string &line, correlation::core::Cell &tempCell,
                        std::vector<correlation::core::Cell> &frames);
};

} // namespace correlation::readers
