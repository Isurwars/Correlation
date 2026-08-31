/**
 * @file AbinitReader.hpp
 * @brief Reader for ABINIT input, output, and trajectory files (.abi, .abo, .abinit, .out).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseReader.hpp"

namespace correlation::readers {

/**
 * @class AbinitReader
 * @brief Handles parsing of ABINIT input/output files containing acell, rprim, and atomic positions.
 */
class AbinitReader : public BaseReader {
public:
  [[nodiscard]] std::string getName() const override { return "ABINIT Reader"; }
  [[nodiscard]] std::vector<std::string> getExtensions() const override { return {".abi", ".abo", ".abinit", ".out"}; }
  [[nodiscard]] bool isTrajectory() const override { return true; }

  correlation::core::Cell
  readStructure(const std::string &filename,
                std::function<void(float, const std::string &)> progress_callback = nullptr) override;

  correlation::core::Trajectory
  readTrajectory(const std::string &filename,
                 std::function<void(float, const std::string &)> progress_callback = nullptr) override;
};

} // namespace correlation::readers
