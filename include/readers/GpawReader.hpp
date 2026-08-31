/**
 * @file GpawReader.hpp
 * @brief Reader for GPAW DFT simulation output and structure files (.gpaw, .txt, .log).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseReader.hpp"

namespace correlation::readers {

/**
 * @class GpawReader
 * @brief Handles parsing of GPAW text log and output files containing unit cells and coordinates.
 */
class GpawReader : public BaseReader {
public:
  [[nodiscard]] std::string getName() const override { return "GPAW Reader"; }
  [[nodiscard]] std::vector<std::string> getExtensions() const override { return {".gpaw", ".txt", ".log"}; }
  [[nodiscard]] bool isTrajectory() const override { return true; }

  correlation::core::Cell
  readStructure(const std::string &filename,
                std::function<void(float, const std::string &)> progress_callback = nullptr) override;

  correlation::core::Trajectory
  readTrajectory(const std::string &filename,
                 std::function<void(float, const std::string &)> progress_callback = nullptr) override;
};

} // namespace correlation::readers
