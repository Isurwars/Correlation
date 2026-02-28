// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "NeighborGraph.hpp"

NeighborGraph::NeighborGraph(size_t node_count) : adj_list_(node_count) {}

void NeighborGraph::addDirectedEdge(size_t from, size_t to, double distance,
                                    const linalg::Vector3<double> &r_ij) {
  if (from >= adj_list_.size()) {
    return;
  }
  adj_list_[from].push_back({static_cast<AtomID>(to), distance, r_ij});
}

const std::vector<Neighbor> &
NeighborGraph::getNeighbors(size_t atom_index) const {
  static const std::vector<Neighbor> empty_list;
  if (atom_index >= adj_list_.size()) {
    return empty_list;
  }
  return adj_list_[atom_index];
}

bool NeighborGraph::areConnected(size_t i, size_t j) const {
  if (i >= adj_list_.size())
    return false;
  for (const auto &neighbor : adj_list_[i]) {
    if (neighbor.index == j) {
      return true;
    }
  }
  return false;
}

std::vector<bool> NeighborGraph::getDenseAdjacencyMatrix() const {
  size_t n = adj_list_.size();
  std::vector<bool> matrix(n * n, false);
  for (size_t i = 0; i < n; ++i) {
    for (const auto &neighbor : adj_list_[i]) {
      if (neighbor.index < n) {
        matrix[i * n + neighbor.index] = true;
      }
    }
  }
  return matrix;
}

size_t NeighborGraph::nodeCount() const { return adj_list_.size(); }
