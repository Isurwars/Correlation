/**
 * @file PoscarWriter.hpp
 * @brief Writer for VASP POSCAR/CONTCAR structure files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "writers/structure/BaseStructureWriter.hpp"

namespace correlation::writers::structure {

/**
 * @brief Writer for VASP POSCAR and CONTCAR files.
 *
 * Implements VASP 5+ format with explicit species symbols, contiguous
 * species grouping, and support for Direct (fractional) or Cartesian coordinates.
 */
class PoscarWriter : public BaseStructureWriter {
public:
  [[nodiscard]] std::string getName() const override { return "VASP POSCAR"; }
  [[nodiscard]] std::vector<std::string> getExtensions() const override {
    return {"vasp", "poscar", "contcar"};
  }

  void write(const std::string &filepath, const correlation::core::Cell &cell,
             const StructureWriterOptions &options = {}) const override;
};

} // namespace correlation::writers::structure
