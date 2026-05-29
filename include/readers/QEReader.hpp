/**
 * @file QEReader.hpp
 * @brief Reader for Quantum ESPRESSO (.pwi, .pwo, .in, .out) files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#pragma once

#include "BaseReader.hpp"

namespace correlation::readers {

/**
 * @class QEReader
 * @brief Handles parsing of Quantum ESPRESSO output and input files.
 * @note Shares ".out" and ".in" extensions with CP2K. Ambiguity is resolved
 * using content sniffing.
 */
class QEReader : public BaseReader {
public:
  std::string getName() const override { return "Quantum ESPRESSO Reader"; }
  std::vector<std::string> getExtensions() const override { return {".pwi", ".pwo", ".in", ".out"}; }
  bool isTrajectory() const override { return true; }

  correlation::core::Cell
  readStructure(const std::string &filename,
                std::function<void(float, const std::string &)> progress_callback = nullptr) override;

  correlation::core::Trajectory
  readTrajectory(const std::string &filename,
                 std::function<void(float, const std::string &)> progress_callback = nullptr) override;
};

} // namespace correlation::readers
