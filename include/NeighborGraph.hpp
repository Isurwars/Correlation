// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "Atom.hpp"
#include <vector>

/**
 * @brief Represents a neighbor atom in a neighbor list.
 */
struct Neighbor {
  AtomID index;                 ///< Index of the neighbor in the atom list
  double distance;              ///< Distance to the neighbor
  linalg::Vector3<double> r_ij; ///< Vector from central atom to neighbor
};

class NeighborGraph {
public:
  explicit NeighborGraph() = default;
  explicit NeighborGraph(size_t node_count);

  void addDirectedEdge(size_t from, size_t to, double distance,
                       const linalg::Vector3<double> &r_ij);

  [[nodiscard]] const std::vector<Neighbor> &
  getNeighbors(size_t atom_index) const;

  [[nodiscard]] bool areConnected(size_t i, size_t j) const;

  [[nodiscard]] std::vector<bool> getDenseAdjacencyMatrix() const;

  [[nodiscard]] size_t nodeCount() const;

private:
  std::vector<std::vector<Neighbor>> adj_list_;
};
