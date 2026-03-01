// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "../NeighborGraph.hpp"
#include <map>
#include <vector>

namespace calculators {

class MotifFinder {
public:
  /**
   * @brief Finds and counts all chordless rings up to a maximum size.
   *
   * @param graph The neighbor graph representing atomic bonds.
   * @param max_size The maximum ring size (number of atoms) to search for.
   * @return A map where the key is the ring size and the value is the total
   * count of such rings.
   */
  static std::map<int, size_t> findRings(const NeighborGraph &graph,
                                         size_t max_size = 6);

  /**
   * @brief Extracts the explicit sequences of AtomIDs that form cycles of
   * `target_size`.
   * @param graph The underlying neighbor graph.
   * @param target_size The exact size of cycle to extract.
   * @return A vector of sequences, where each sequence is a list of AtomIDs
   * forming the cycle.
   */
  /**
   * @brief Extracts all exact cycles of a specific target size.
   *
   * @param graph The neighbor graph to search.
   * @param target_size The exact size of the rings to extract.
   * @return A vector of rings, where each ring is represented as a vector of
   * AtomIDs in order.
   */
  static std::vector<std::vector<AtomID>>
  extractCycles(const NeighborGraph &graph, size_t target_size);
};

} // namespace calculators
