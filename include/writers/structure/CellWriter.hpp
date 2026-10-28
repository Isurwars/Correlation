/**
 * @file CellWriter.hpp
 * @brief Writer for CASTEP .cell crystal structure files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "writers/structure/BaseStructureWriter.hpp"

namespace correlation::writers::structure {

/**
 * @brief Writer for CASTEP .cell structure files.
 *
 * Emits `%BLOCK LATTICE_CART` for lattice vectors and either
 * `%BLOCK POSITIONS_FRAC` or `%BLOCK POSITIONS_ABS` for coordinates.
 */
class CellWriter : public BaseStructureWriter {
public:
  [[nodiscard]] std::string getName() const override { return "CASTEP Cell"; }
  [[nodiscard]] std::vector<std::string> getExtensions() const override { return {"cell"}; }

  void write(const std::string &filepath, const correlation::core::Cell &cell,
             const StructureWriterOptions &options = {}) const override;
};

} // namespace correlation::writers::structure
