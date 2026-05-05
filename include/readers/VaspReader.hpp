/**
 * @file VaspReader.hpp
 * @brief Reader for VASP POSCAR/CONTCAR structure files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#pragma once

#include "BaseReader.hpp"

#include <string>

namespace correlation::readers {

/**
 * @brief Reads VASP POSCAR and CONTCAR structure files.
 *
 * Supports VASP 5+ format with explicit species names.
 * Handles both Direct (fractional) and Cartesian coordinate types,
 * scaling factors (positive = multiplier, negative = target volume),
 * and optional Selective Dynamics lines.
 */
class VaspReader : public BaseReader {
public:
  std::string getName() const override { return "VASP POSCAR"; }
  std::vector<std::string> getExtensions() const override {
    return {"poscar", "contcar", "vasp"};
  }
  bool isTrajectory() const override { return false; }

  correlation::core::Cell readStructure(
      const std::string &filename,
      std::function<void(float, const std::string &)> progress_callback =
          nullptr) override;

  correlation::core::Trajectory readTrajectory(
      const std::string &filename,
      std::function<void(float, const std::string &)> progress_callback =
          nullptr) override;

  /**
   * @brief Low-level parser for POSCAR/CONTCAR files.
   * @param file_name Path to the source file.
   * @return The parsed simulation cell.
   */
  static correlation::core::Cell read(const std::string &file_name);
};

} // namespace correlation::readers
