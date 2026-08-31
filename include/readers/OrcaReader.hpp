/**
 * @file OrcaReader.hpp
 * @brief Reader for ORCA quantum chemistry output files (.orca, .out, .log).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseReader.hpp"

namespace correlation::readers {

/**
 * @class OrcaReader
 * @brief Handles parsing of ORCA output and log files for structures and geometry optimization trajectories.
 */
class OrcaReader : public BaseReader {
public:
  [[nodiscard]] std::string getName() const override { return "ORCA Reader"; }
  [[nodiscard]] std::vector<std::string> getExtensions() const override { return {".orca", ".out", ".log"}; }
  [[nodiscard]] bool isTrajectory() const override { return true; }

  correlation::core::Cell
  readStructure(const std::string &filename,
                std::function<void(float, const std::string &)> progress_callback = nullptr) override;

  correlation::core::Trajectory
  readTrajectory(const std::string &filename,
                 std::function<void(float, const std::string &)> progress_callback = nullptr) override;
};

} // namespace correlation::readers
