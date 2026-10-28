/**
 * @file CifWriter.hpp
 * @brief Writer for Crystallographic Information Files (.cif).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "writers/structure/BaseStructureWriter.hpp"

namespace correlation::writers::structure {

/**
 * @brief Writer for Crystallographic Information Files (CIF) in P1 symmetry.
 *
 * Serializes cell lengths, cell angles, cell volume, and fractional atom site
 * coordinates adhering to IUCr core CIF standards.
 */
class CifWriter : public BaseStructureWriter {
public:
  [[nodiscard]] std::string getName() const override { return "CIF"; }
  [[nodiscard]] std::vector<std::string> getExtensions() const override { return {"cif"}; }

  void write(const std::string &filepath, const correlation::core::Cell &cell,
             const StructureWriterOptions &options = {}) const override;
};

} // namespace correlation::writers::structure
