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
   * @brief Uses Bounded Depth-First Search (DFS) to find all unique simple
   * cycles up to `max_size`.
   * @param graph The underlying neighbor graph.
   * @param max_size The maximum cycle length to search for (e.g. 6 for
   * hexagons).
   * @return A map of cycle length to the number of unique occurrences.
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
  static std::vector<std::vector<AtomID>>
  extractCycles(const NeighborGraph &graph, size_t target_size);
};

} // namespace calculators
