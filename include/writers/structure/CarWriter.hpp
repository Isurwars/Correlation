/**
 * @file CarWriter.hpp
 * @brief Writer for Accelrys / Materials Studio CAR structure files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "writers/structure/BaseStructureWriter.hpp"

namespace correlation::writers::structure {

/**
 * @brief Writer for Accelrys .car files with periodic boundary conditions.
 *
 * Emits `!BIOSYM archive 3`, `PBC=ON`, the lattice parameter block, and
 * 8-column atom records with Cartesian coordinates.
 */
class CarWriter : public BaseStructureWriter {
public:
  [[nodiscard]] std::string getName() const override { return "Accelrys CAR"; }
  [[nodiscard]] std::vector<std::string> getExtensions() const override { return {"car"}; }

  void write(const std::string &filepath, const correlation::core::Cell &cell,
             const StructureWriterOptions &options = {}) const override;
};

} // namespace correlation::writers::structure
