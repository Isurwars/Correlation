/**
 * @file DftbReader.hpp
 * @brief Reader for DFTB+ geometry files (.gen, .dftb, .hsd).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseReader.hpp"

namespace correlation::readers {

/**
 * @class DftbReader
 * @brief Handles parsing of DFTB+ GenFormat (.gen) and HSD input geometry structures.
 */
class DftbReader : public BaseReader {
public:
  [[nodiscard]] std::string getName() const override { return "DFTB+ Reader"; }
  [[nodiscard]] std::vector<std::string> getExtensions() const override { return {".gen", ".dftb", ".hsd"}; }
  [[nodiscard]] bool isTrajectory() const override { return false; }

  correlation::core::Cell
  readStructure(const std::string &filename,
                std::function<void(float, const std::string &)> progress_callback = nullptr) override;

  correlation::core::Trajectory
  readTrajectory(const std::string &filename,
                 std::function<void(float, const std::string &)> progress_callback = nullptr) override;
};

} // namespace correlation::readers
