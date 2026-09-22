/**
 * @file BondCutoffMapper.hpp
 * @brief Decoupled mapper for atomic bond cutoffs between UI models and domain matrices.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
#include "core/Atom.hpp"

#include <span>
#include <string>
#include <vector>

namespace correlation::app {

/**
 * @struct CutoffEntry
 * @brief Plain data representation of a bond cutoff pair decoupled from UI frameworks.
 */
struct CutoffEntry {
  std::string element1;     ///< First element symbol.
  std::string element2;     ///< Second element symbol.
  std::string min_distance; ///< Minimum cutoff distance string.
  std::string max_distance; ///< Maximum cutoff distance string.
};

/**
 * @class BondCutoffMapper
 * @brief Stateless transformation service mapping atomic elements to cutoff configurations.
 */
class BondCutoffMapper {
public:
  /**
   * @brief Generates default covalent cutoff entries for a collection of elements.
   * @param[in] elements Collection of unique elements in the system.
   * @return Vector of formatted CutoffEntry objects with estimated covalent bounds.
   */
  [[nodiscard]] static std::vector<CutoffEntry>
  createDefaultCutoffEntries(std::span<const correlation::core::Element> elements);

  /**
   * @brief Parses and validates cutoff entries into a symmetric BondCutoffMatrix.
   * @param[in] entries Collection of cutoff entries to parse.
   * @param[in] elements Collection of elements defining matrix dimension and index order.
   * @return Symmetrized BondCutoffMatrix with squared distance ranges.
   */
  [[nodiscard]] static correlation::analysis::BondCutoffMatrix
  parseCutoffMatrix(std::span<const CutoffEntry> entries,
                    std::span<const correlation::core::Element> elements);

private:
  /**
   * @brief Safely parses a distance string to real_t with a zero fallback.
   * @param[in] str Distance string to parse.
   * @return Parsed distance value or 0.0 on conversion failure.
   */
  [[nodiscard]] static real_t parseDistanceSafe(const std::string &str);

  /**
   * @brief Finds the index of an element symbol within an element span.
   * @param[in] elements Collection of elements.
   * @param[in] symbol Symbol to look up.
   * @return Index if found, or -1 if not found.
   */
  [[nodiscard]] static int findElementIndex(std::span<const correlation::core::Element> elements,
                                           const std::string &symbol);
};

} // namespace correlation::app
