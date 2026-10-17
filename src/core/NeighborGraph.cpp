/**
 * @file NeighborGraph.cpp
 * @brief Implementation of the Compressed Sparse Row (CSR) atomic neighbor graph.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/NeighborGraph.hpp"

#include <algorithm>
#include <cassert>

namespace correlation::core {

namespace {

/**
 * @brief Compares two Neighbor instances deterministically by (index, distance, r_ij).
 */
[[nodiscard]] constexpr bool compareNeighbors(const Neighbor &lhs, const Neighbor &rhs) noexcept {
  if (lhs.index != rhs.index) {
    return lhs.index < rhs.index;
  }
  if (lhs.distance != rhs.distance) {
    return lhs.distance < rhs.distance;
  }
  if (lhs.r_ij.x() != rhs.r_ij.x()) {
    return lhs.r_ij.x() < rhs.r_ij.x();
  }
  if (lhs.r_ij.y() != rhs.r_ij.y()) {
    return lhs.r_ij.y() < rhs.r_ij.y();
  }
  return lhs.r_ij.z() < rhs.r_ij.z();
}

} // namespace

NeighborGraph::NeighborGraph(size_t node_count)
    : node_count_(node_count), is_finalized_(true), offsets_(node_count + 1, 0) {}

void NeighborGraph::addDirectedEdge(const Edge &edge) {
  assert(edge.source < node_count_ &&
         "Source node index out of bounds in NeighborGraph::addDirectedEdge");
  if (edge.source >= node_count_) {
    return;
  }

  if (is_finalized_ && !neighbors_.empty()) {
    for (size_t u = 0; u < node_count_; ++u) {
      const size_t start = offsets_[u];
      const size_t end = offsets_[u + 1];
      for (size_t k = start; k < end; ++k) {
        unfinalized_edges_.push_back(Edge{
            .source = static_cast<AtomID>(u),
            .target = neighbors_[k].index,
            .distance = neighbors_[k].distance,
            .r_ij = neighbors_[k].r_ij,
        });
      }
    }
    neighbors_.clear();
    targets_.clear();
  }

  is_finalized_ = false;
  unfinalized_edges_.push_back(edge);
}

void NeighborGraph::finalize() const {
  if (is_finalized_) {
    return;
  }

  std::vector<size_t> row_counts(node_count_, 0);
  for (const auto &edge : unfinalized_edges_) {
    if (edge.source < node_count_) {
      row_counts[edge.source]++;
    }
  }

  offsets_.assign(node_count_ + 1, 0);
  for (size_t i = 0; i < node_count_; ++i) {
    offsets_[i + 1] = offsets_[i] + row_counts[i];
  }

  const size_t total_edges = offsets_[node_count_];
  neighbors_.resize(total_edges);
  targets_.resize(total_edges);

  std::vector<size_t> write_pos = offsets_;
  for (const auto &edge : unfinalized_edges_) {
    if (edge.source < node_count_) {
      const size_t slot = write_pos[edge.source]++;
      neighbors_[slot] = Neighbor{
          .index = edge.target,
          .distance = edge.distance,
          .r_ij = edge.r_ij,
      };
    }
  }

  for (size_t u = 0; u < node_count_; ++u) {
    const size_t start = offsets_[u];
    const size_t end = offsets_[u + 1];
    std::sort(neighbors_.begin() + static_cast<std::ptrdiff_t>(start),
              neighbors_.begin() + static_cast<std::ptrdiff_t>(end), compareNeighbors);
    for (size_t k = start; k < end; ++k) {
      targets_[k] = neighbors_[k].index;
    }
  }

  unfinalized_edges_.clear();
  unfinalized_edges_.shrink_to_fit();
  is_finalized_ = true;
}

void NeighborGraph::ensureFinalized() const {
  if (!is_finalized_) {
    finalize();
  }
}

void NeighborGraph::assignCSR(size_t node_count, std::vector<size_t> offsets,
                              std::vector<Neighbor> neighbors, std::vector<AtomID> targets) {
  node_count_ = node_count;
  offsets_ = std::move(offsets);
  neighbors_ = std::move(neighbors);
  targets_ = std::move(targets);
  unfinalized_edges_.clear();
  unfinalized_edges_.shrink_to_fit();
  is_finalized_ = true;
}

std::span<const Neighbor> NeighborGraph::getNeighbors(size_t atom_index) const {
  ensureFinalized();
  if (atom_index >= node_count_) {
    return {};
  }
  const size_t start = offsets_[atom_index];
  const size_t count = offsets_[atom_index + 1] - start;
  return std::span<const Neighbor>(neighbors_.data() + start, count);
}

std::span<const AtomID> NeighborGraph::getNeighborIndices(size_t atom_index) const {
  ensureFinalized();
  if (atom_index >= node_count_) {
    return {};
  }
  const size_t start = offsets_[atom_index];
  const size_t count = offsets_[atom_index + 1] - start;
  return std::span<const AtomID>(targets_.data() + start, count);
}

bool NeighborGraph::areConnected(AtomIndex first_atom, AtomIndex second_atom) const {
  ensureFinalized();
  if (first_atom.id >= node_count_) {
    return false;
  }
  const auto target_id = static_cast<AtomID>(second_atom.id);
  const size_t start = offsets_[first_atom.id];
  const size_t end = offsets_[first_atom.id + 1];
  return std::binary_search(targets_.begin() + static_cast<std::ptrdiff_t>(start),
                            targets_.begin() + static_cast<std::ptrdiff_t>(end), target_id);
}

std::vector<bool> NeighborGraph::getDenseAdjacencyMatrix() const {
  ensureFinalized();
  const size_t num_nodes = node_count_;
  constexpr size_t K_MAX_DENSE_NODES = 20000;
  if (num_nodes > K_MAX_DENSE_NODES) {
    throw std::runtime_error("NeighborGraph::getDenseAdjacencyMatrix: graph size (" +
                             std::to_string(num_nodes) +
                             ") exceeds maximum safe limit for dense matrix allocation (" +
                             std::to_string(K_MAX_DENSE_NODES) + ").");
  }
  std::vector<bool> matrix(num_nodes * num_nodes, false);
  for (size_t i = 0; i < num_nodes; ++i) {
    const size_t start = offsets_[i];
    const size_t end = offsets_[i + 1];
    for (size_t k = start; k < end; ++k) {
      const AtomID target = targets_[k];
      if (target < num_nodes) {
        matrix[i * num_nodes + target] = true;
      }
    }
  }
  return matrix;
}

size_t NeighborGraph::nodeCount() const noexcept { return node_count_; }

size_t NeighborGraph::edgeCount() const noexcept {
  ensureFinalized();
  return targets_.size();
}

std::span<const size_t> NeighborGraph::offsets() const noexcept {
  ensureFinalized();
  return offsets_;
}

std::span<const AtomID> NeighborGraph::targets() const noexcept {
  ensureFinalized();
  return targets_;
}

std::span<const Neighbor> NeighborGraph::neighbors() const noexcept {
  ensureFinalized();
  return neighbors_;
}

} // namespace correlation::core
