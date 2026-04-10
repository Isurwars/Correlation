/**
 * @file NeighborGraph.hpp
 * @brief Neighbor graph representation of atomic bonds and adjacency.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "core/Atom.hpp"
#include "math/LinearAlgebra.hpp"

#include <vector>

namespace correlation::core {

/**
 * @brief Represents a neighbor atom in a neighbor list.
 */
struct Neighbor {
  AtomID index;               ///< Index of the neighbor in the atom list
  double distance;            ///< Distance to the neighbor
  math::Vector3<double> r_ij; ///< Vector from central atom to neighbor
};

class NeighborGraph {
public:
  explicit NeighborGraph() = default;
  explicit NeighborGraph(size_t node_count);

  void addDirectedEdge(size_t from, size_t to, double distance,
                       const math::Vector3<double> &r_ij);

  [[nodiscard]] const std::vector<Neighbor> &
  getNeighbors(size_t atom_index) const;

  [[nodiscard]] bool areConnected(size_t i, size_t j) const;

  [[nodiscard]] std::vector<bool> getDenseAdjacencyMatrix() const;

  [[nodiscard]] size_t nodeCount() const;

private:
  std::vector<std::vector<Neighbor>> adj_list_;
};

} // namespace correlation::core
