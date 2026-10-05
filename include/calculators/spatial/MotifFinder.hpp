/**
 * @file MotifFinder.hpp
 * @brief Ring and structural motif finder using neighbour graph search.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include <map>
#include <string_view>
#include <vector>

namespace correlation::calculators {

using correlation::analysis::RingProjectionMode;
using correlation::analysis::RingType;

/**
 * @class MotifFinder
 * @brief Utility class for finding structural motifs (rings) in a system.
 *
 * This class provides static methods to search for and count chordless rings
 * within a molecular system, leveraging the neighbor graph structure.
 * Supports both King shortest-distance rings and Franzblau unique-geodesic primitive rings.
 */
class MotifFinder {
public:
  /**
   * @brief Finds and counts all chordless rings up to a maximum size using King's criterion.
   *
   * @param[in] graph The neighbor graph representing atomic bonds.
   * @param[in] max_size The maximum ring size (number of atoms) to search for.
   * @return A map where the key is the ring size and the value is the total count.
   */
  static std::map<int, size_t> findRings(const correlation::core::NeighborGraph &graph,
                                         size_t max_size = 6);

  /**
   * @brief Finds and counts all Franzblau shortest-path primitive rings up to a maximum size.
   *
   * @param[in] graph The neighbor graph representing atomic bonds.
   * @param[in] max_size The maximum ring size (number of atoms) to search for.
   * @return A map where the key is the ring size and the value is the total count.
   */
  static std::map<int, size_t> findFranzblauRings(const correlation::core::NeighborGraph &graph,
                                                  size_t max_size = 6);

  /**
   * @brief Finds and counts all rings according to the specified ring criterion.
   *
   * @param[in] graph The neighbor graph representing atomic bonds.
   * @param[in] max_size The maximum ring size to search for.
   * @param[in] ring_type The ring criterion (King or Franzblau).
   * @return A map where the key is the ring size and the value is the total count.
   */
  static std::map<int, size_t> findRings(const correlation::core::NeighborGraph &graph,
                                         size_t max_size, RingType ring_type);

  /**
   * @brief Extracts all exact cycles of a specific target size using King's criterion.
   *
   * @param[in] graph The neighbor graph to search.
   * @param[in] target_size The exact size of the rings to extract.
   * @return A vector of rings, where each ring is represented as a vector of AtomIDs in order.
   */
  static std::vector<std::vector<correlation::core::AtomID>>
  extractCycles(const correlation::core::NeighborGraph &graph, size_t target_size);

  /**
   * @brief Extracts all exact cycles of a specific target size using Franzblau's criterion.
   *
   * @param[in] graph The neighbor graph to search.
   * @param[in] target_size The exact size of the rings to extract.
   * @return A vector of rings, where each ring is represented as a vector of AtomIDs in order.
   */
  static std::vector<std::vector<correlation::core::AtomID>>
  extractFranzblauCycles(const correlation::core::NeighborGraph &graph, size_t target_size);

  /**
   * @brief Extracts all exact cycles of a specific target size according to the ring criterion.
   *
   * @param[in] graph The neighbor graph to search.
   * @param[in] target_size The exact size of the rings to extract.
   * @param[in] ring_type The ring criterion (King or Franzblau).
   * @return A vector of rings.
   */
  static std::vector<std::vector<correlation::core::AtomID>>
  extractCycles(const correlation::core::NeighborGraph &graph, size_t target_size,
                RingType ring_type);

  /**
   * @brief Extracts all exact cycles up to a maximum size according to the ring criterion.
   *
   * @param[in] graph The neighbor graph to search.
   * @param[in] max_size The maximum ring size to extract.
   * @param[in] ring_type The ring criterion (King or Franzblau).
   * @return A vector of all detected rings.
   */
  static std::vector<std::vector<correlation::core::AtomID>>
  extractAllCycles(const correlation::core::NeighborGraph &graph, size_t max_size,
                   RingType ring_type = RingType::King);

  /**
   * @brief Builds a projected neighbor graph collapsing bridging atoms (e.g. Si-O-Si -> Si-Si).
   *
   * In network glasses such as SiO2, this constructs an adjacency graph purely among network-former
   * atoms that share a common bridging atom.
   *
   * @param[in] graph The underlying atomic neighbor graph.
   * @param[in] cell The atomic cell containing element species and positions.
   * @param[in] former_element Chemical symbol of network-former (e.g. "Si").
   * @param[in] bridging_element Chemical symbol of bridging species (e.g. "O").
   * @return A new NeighborGraph linking former atoms.
   */
  static correlation::core::NeighborGraph
  buildBridgedGraph(const correlation::core::NeighborGraph &graph,
                    const correlation::core::Cell &cell, std::string_view former_element,
                    std::string_view bridging_element);

  /**
   * @brief Filters cycle list to retain only cycles strictly alternating between two elements.
   *
   * For example in SiO2, retains rings whose sequence is Si-O-Si-O...
   *
   * @param[in] cycles List of extracted atomic cycles.
   * @param[in] cell The atomic cell containing element species.
   * @param[in] element_a First element in alternating sequence.
   * @param[in] element_b Second element in alternating sequence.
   * @return Filtered vector of alternating cycles.
   */
  static std::vector<std::vector<correlation::core::AtomID>>
  filterAlternatingCycles(const std::vector<std::vector<correlation::core::AtomID>> &cycles,
                          const correlation::core::Cell &cell, std::string_view element_a,
                          std::string_view element_b);
};

} // namespace correlation::calculators
