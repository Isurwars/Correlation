/**
 * @file CP2KReader.hpp
 * @brief Reader for CP2K (.inp, .restart, .out) files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#pragma once

#include "BaseReader.hpp"

namespace correlation::readers {

/**
 * @class CP2KReader
 * @brief Handles parsing of CP2K restart and input files.
 * @note Shares ".out" extension with Quantum ESPRESSO. Ambiguity is resolved
 * using content sniffing.
 */
class CP2KReader : public BaseReader {
public:
  std::string getName() const override { return "CP2K Reader"; }
  std::vector<std::string> getExtensions() const override { return {".inp", ".restart", ".out", ".cp2k"}; }
  bool isTrajectory() const override { return true; }

  correlation::core::Cell
  readStructure(const std::string &filename,
                std::function<void(float, const std::string &)> progress_callback = nullptr) override;

  correlation::core::Trajectory
  readTrajectory(const std::string &filename,
                 std::function<void(float, const std::string &)> progress_callback = nullptr) override;
};

} // namespace correlation::readers
