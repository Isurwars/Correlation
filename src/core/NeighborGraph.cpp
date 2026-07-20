/**
 * @file NeighborGraph.cpp
 * @brief Implementation of the atomic neighbor graph.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/NeighborGraph.hpp"
#include "math/LinearAlgebra.hpp"

#include <algorithm>

namespace correlation::core {

NeighborGraph::NeighborGraph(size_t node_count) : adj_list_(node_count) {}

void NeighborGraph::addDirectedEdge(size_t source, size_t target, real_t distance,
                                    const math::Vector3<real_t> &r_ij) { // NOLINT(bugprone-easily-swappable-parameters)
  if (source >= adj_list_.size()) {
    return;
  }
  adj_list_[source].push_back({.index = static_cast<AtomID>(target), .distance = distance, .r_ij = r_ij});
}

const std::vector<Neighbor> &NeighborGraph::getNeighbors(size_t atom_index) const {
  static const std::vector<Neighbor> empty_list;
  if (atom_index >= adj_list_.size()) {
    return empty_list;
  }
  return adj_list_[atom_index];
}

bool NeighborGraph::areConnected(AtomIndex first_atom, AtomIndex second_atom) const {
  if (first_atom.id >= adj_list_.size()) {
    return false;
  }
  auto target_id = static_cast<AtomID>(second_atom.id);
  return std::ranges::any_of(adj_list_[first_atom.id],
                             [target_id](const auto &neighbor) { return neighbor.index == target_id; });
}

std::vector<bool> NeighborGraph::getDenseAdjacencyMatrix() const {
  size_t const num_nodes = adj_list_.size();
  std::vector<bool> matrix(num_nodes * num_nodes, false);
  for (size_t i = 0; i < num_nodes; ++i) {
    for (const auto &neighbor : adj_list_[i]) {
      if (neighbor.index < num_nodes) {
        matrix[i * num_nodes + neighbor.index] = true;
      }
    }
  }
  return matrix;
}

size_t NeighborGraph::nodeCount() const { return adj_list_.size(); }

} // namespace correlation::core
