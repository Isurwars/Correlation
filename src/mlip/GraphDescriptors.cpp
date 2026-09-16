/**
 * @file GraphDescriptors.cpp
 * @brief Implementation of topological, structural, and spectral graph descriptors.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "mlip/GraphDescriptors.hpp"
#include "calculators/MotifFinder.hpp"
#include "core/NeighborGraph.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace correlation::mlip {

namespace {

correlation::core::NeighborGraph buildNeighborGraphFromBuffers(const PeriodicGraphData &graph) {
  correlation::core::NeighborGraph neighbor_graph(graph.atom_count);
  const size_t total_edges = graph.edge_count;
  const bool has_distances = graph.edge_distances.size() == total_edges;
  const bool has_vectors = graph.edge_vectors_flat.size() == total_edges * 3;

  for (size_t edge_idx = 0; edge_idx < total_edges; ++edge_idx) {
    const auto src = static_cast<correlation::core::AtomID>(graph.edge_index_flat[edge_idx]);
    const auto dst =
        static_cast<correlation::core::AtomID>(graph.edge_index_flat[total_edges + edge_idx]);
    const real_t dist = has_distances ? graph.edge_distances[edge_idx] : static_cast<real_t>(0.0);

    correlation::math::Vector3<real_t> r_ij{0.0, 0.0, 0.0};
    if (has_vectors) {
      r_ij = correlation::math::Vector3<real_t>{graph.edge_vectors_flat[edge_idx * 3],
                                                graph.edge_vectors_flat[edge_idx * 3 + 1],
                                                graph.edge_vectors_flat[edge_idx * 3 + 2]};
    }
    neighbor_graph.addDirectedEdge(src, dst, dist, r_ij);
  }
  return neighbor_graph;
}

std::vector<std::set<size_t>> buildAdjacencySets(const PeriodicGraphData &graph) {
  std::vector<std::set<size_t>> adj(graph.atom_count);
  const size_t total_edges = graph.edge_count;
  for (size_t edge_idx = 0; edge_idx < total_edges; ++edge_idx) {
    const auto src = static_cast<size_t>(graph.edge_index_flat[edge_idx]);
    const auto dst = static_cast<size_t>(graph.edge_index_flat[total_edges + edge_idx]);
    if (src < graph.atom_count && dst < graph.atom_count && src != dst) {
      adj[src].insert(dst);
      adj[dst].insert(src);
    }
  }
  return adj;
}

struct DfsStackFrame {
  size_t node;
  size_t neighbor_idx;
  size_t max_child_len;
};

// Stack-based iterative DFS to avoid recursion
size_t dfsLongestPathIterative(size_t start_node,
                               const std::map<size_t, std::vector<size_t>> &adj) {
  std::set<size_t> visited;
  visited.insert(start_node);

  std::vector<DfsStackFrame> stack;
  stack.reserve(16);
  stack.push_back({
      .node = start_node,
      .neighbor_idx = 0,
      .max_child_len = 0,
  });

  size_t final_best = 0;

  while (!stack.empty()) {
    const size_t current_node = stack.back().node;
    const size_t neighbor_idx = stack.back().neighbor_idx;
    const auto adj_iter = adj.find(current_node);

    if (adj_iter == adj.end() || neighbor_idx >= adj_iter->second.size()) {
      const size_t max_child_len = stack.back().max_child_len;
      stack.pop_back();
      visited.erase(current_node);

      if (!stack.empty()) {
        stack.back().max_child_len = std::max(stack.back().max_child_len, 1 + max_child_len);
      } else {
        final_best = max_child_len;
      }
    } else {
      const size_t neighbor = adj_iter->second[neighbor_idx];
      stack.back().neighbor_idx++;

      if (!visited.contains(neighbor)) {
        visited.insert(neighbor);
        stack.push_back({
            .node = neighbor,
            .neighbor_idx = 0,
            .max_child_len = 0,
        });
      }
    }
  }

  return final_best;
}

size_t computeLongestChain(const std::vector<size_t> &common,
                           const std::map<size_t, std::vector<size_t>> &adj) {
  size_t best = 0;
  for (const size_t start : common) {
    best = std::max(best, dfsLongestPathIterative(start, adj));
  }
  return best;
}

struct CNAPairSig {
  size_t n_common{0};
  size_t n_bonds{0};
  size_t n_chain{0};

  bool operator==(const CNAPairSig &other) const noexcept {
    return n_common == other.n_common && n_bonds == other.n_bonds && n_chain == other.n_chain;
  }
};

CNAPairSig evaluatePairSignature(size_t atom_i, size_t atom_j,
                                 const std::vector<std::set<size_t>> &adj) {
  std::vector<size_t> common;
  for (const size_t nbr : adj[atom_i]) {
    if (adj[atom_j].contains(nbr)) {
      common.push_back(nbr);
    }
  }

  size_t bonds = 0;
  std::map<size_t, std::vector<size_t>> sub_adj;
  const size_t n_common = common.size();
  for (size_t idx_a = 0; idx_a < n_common; ++idx_a) {
    for (size_t idx_b = idx_a + 1; idx_b < n_common; ++idx_b) {
      if (adj[common[idx_a]].contains(common[idx_b])) {
        ++bonds;
        sub_adj[common[idx_a]].push_back(common[idx_b]);
        sub_adj[common[idx_b]].push_back(common[idx_a]);
      }
    }
  }

  size_t chain = 0;
  if (bonds > 0) {
    chain = computeLongestChain(common, sub_adj);
  }
  return CNAPairSig{.n_common = n_common, .n_bonds = bonds, .n_chain = chain};
}

struct MotifSignatureCounts {
  size_t coordination{0};
  size_t count_421{0};
  size_t count_422{0};
  size_t count_666{0};
  size_t count_444{0};
  size_t count_555{0};
};

CNALabel classifyEnvironment(const MotifSignatureCounts &counts) noexcept {
  if (counts.coordination == 12 && counts.count_421 == 12) {
    return CNALabel::FCC;
  }
  if (counts.coordination == 12 && counts.count_421 == 6 && counts.count_422 == 6) {
    return CNALabel::HCP;
  }
  if ((counts.coordination == 14 && counts.count_666 == 8) ||
      (counts.coordination == 8 && counts.count_666 == 8)) {
    return CNALabel::BCC;
  }
  if (counts.coordination == 12 && counts.count_555 == 12) {
    return CNALabel::ICO;
  }
  return CNALabel::Other;
}

CNALabel determineAtomMotif(size_t atom_idx, const std::vector<std::set<size_t>> &adj) {
  const size_t coord = adj[atom_idx].size();
  if (coord < 8) {
    return CNALabel::Other;
  }

  MotifSignatureCounts counts{.coordination = coord};

  for (const size_t nbr : adj[atom_idx]) {
    const auto sig = evaluatePairSignature(atom_idx, nbr, adj);
    if (sig.n_common == 4 && sig.n_bonds == 2 && sig.n_chain <= 1) {
      ++counts.count_421;
    } else if (sig.n_common == 4 && sig.n_bonds == 2 && sig.n_chain >= 2) {
      ++counts.count_422;
    } else if (sig.n_common == 6 && sig.n_bonds == 6) {
      ++counts.count_666;
    } else if (sig.n_common == 4 && sig.n_bonds == 4) {
      ++counts.count_444;
    } else if (sig.n_common == 5 && sig.n_bonds == 5) {
      ++counts.count_555;
    }
  }

  return classifyEnvironment(counts);
}

// Symmetric matrix diagonalization via Givens / Jacobi rotations
std::vector<real_t> computeSymmetricEigenvalues(std::vector<real_t> mat, size_t matrix_dim) {
  std::vector<real_t> eigenvalues(matrix_dim, 0.0);
  for (size_t idx = 0; idx < matrix_dim; ++idx) {
    eigenvalues[idx] = mat[idx * matrix_dim + idx];
  }

  const size_t max_iterations = 50;
  for (size_t iter = 0; iter < max_iterations; ++iter) {
    real_t max_off_diag = 0.0;
    size_t idx_p = 0;
    size_t idx_q = 1;
    for (size_t i = 0; i < matrix_dim; ++i) {
      for (size_t j = i + 1; j < matrix_dim; ++j) {
        const real_t off = std::abs(mat[i * matrix_dim + j]);
        if (off > max_off_diag) {
          max_off_diag = off;
          idx_p = i;
          idx_q = j;
        }
      }
    }

    if (max_off_diag < static_cast<real_t>(1e-9)) {
      break;
    }

    const real_t app = mat[idx_p * matrix_dim + idx_p];
    const real_t aqq = mat[idx_q * matrix_dim + idx_q];
    const real_t apq = mat[idx_p * matrix_dim + idx_q];
    const real_t theta = static_cast<real_t>(0.5) * (aqq - app) / apq;
    const real_t t_val = (theta >= 0.0)
                             ? static_cast<real_t>(1.0) /
                                   (theta + std::sqrt(theta * theta + static_cast<real_t>(1.0)))
                             : static_cast<real_t>(-1.0) /
                                   (-theta + std::sqrt(theta * theta + static_cast<real_t>(1.0)));
    const real_t cos_val =
        static_cast<real_t>(1.0) / std::sqrt(t_val * t_val + static_cast<real_t>(1.0));
    const real_t sin_val = t_val * cos_val;
    const real_t tau_val = sin_val / (static_cast<real_t>(1.0) + cos_val);

    mat[idx_p * matrix_dim + idx_p] = app - t_val * apq;
    mat[idx_q * matrix_dim + idx_q] = aqq + t_val * apq;
    mat[idx_p * matrix_dim + idx_q] = 0.0;
    mat[idx_q * matrix_dim + idx_p] = 0.0;

    for (size_t idx_r = 0; idx_r < matrix_dim; ++idx_r) {
      if (idx_r != idx_p && idx_r != idx_q) {
        const real_t arp = mat[idx_r * matrix_dim + idx_p];
        const real_t arq = mat[idx_r * matrix_dim + idx_q];
        mat[idx_r * matrix_dim + idx_p] = arp - sin_val * (arq + tau_val * arp);
        mat[idx_p * matrix_dim + idx_r] = mat[idx_r * matrix_dim + idx_p];
        mat[idx_r * matrix_dim + idx_q] = arq + sin_val * (arp - tau_val * arq);
        mat[idx_q * matrix_dim + idx_r] = mat[idx_r * matrix_dim + idx_q];
      }
    }
  }

  for (size_t idx = 0; idx < matrix_dim; ++idx) {
    eigenvalues[idx] = mat[idx * matrix_dim + idx];
  }
  std::ranges::sort(eigenvalues, std::greater<>{});
  return eigenvalues;
}

} // anonymous namespace

std::vector<real_t>
GraphDescriptors::computeRingStatisticsDescriptor(const PeriodicGraphData &graph, size_t max_size) {
  if (graph.atom_count == 0 || max_size < 3) {
    return {};
  }

  std::vector<real_t> ring_desc(graph.atom_count * max_size, static_cast<real_t>(0.0));
  const auto neighbor_graph = buildNeighborGraphFromBuffers(graph);

  for (size_t size = 3; size <= max_size; ++size) {
    const auto cycles = correlation::calculators::MotifFinder::extractCycles(neighbor_graph, size);
    for (const auto &cycle : cycles) {
      for (const auto atom_id : cycle) {
        if (static_cast<size_t>(atom_id) < graph.atom_count) {
          ring_desc[static_cast<size_t>(atom_id) * max_size + (size - 1)] +=
              static_cast<real_t>(1.0);
        }
      }
    }
  }
  return ring_desc;
}

std::vector<int> GraphDescriptors::computeCNADescriptor(const PeriodicGraphData &graph) {
  if (graph.atom_count == 0) {
    return {};
  }

  const auto adj = buildAdjacencySets(graph);
  std::vector<int> labels(graph.atom_count, 0);

  for (size_t i = 0; i < graph.atom_count; ++i) {
    labels[i] = static_cast<int>(determineAtomMotif(i, adj));
  }
  return labels;
}

std::vector<real_t> GraphDescriptors::computeCoordinationEmbedding(const PeriodicGraphData &graph) {
  if (graph.atom_count == 0) {
    return {};
  }

  std::vector<real_t> coord(graph.atom_count, static_cast<real_t>(0.0));
  const size_t total_edges = graph.edge_count;
  for (size_t edge_idx = 0; edge_idx < total_edges; ++edge_idx) {
    const auto src = static_cast<size_t>(graph.edge_index_flat[edge_idx]);
    if (src < graph.atom_count) {
      coord[src] += static_cast<real_t>(1.0);
    }
  }
  return coord;
}

std::vector<real_t> GraphDescriptors::computeGraphSpectrum(const PeriodicGraphData &graph,
                                                           size_t k_eigenvalues) {
  const size_t num_atoms = graph.atom_count;
  if (num_atoms == 0 || k_eigenvalues == 0) {
    return {};
  }

  std::vector<real_t> mat(num_atoms * num_atoms, static_cast<real_t>(0.0));
  const size_t total_edges = graph.edge_count;
  for (size_t edge_idx = 0; edge_idx < total_edges; ++edge_idx) {
    const auto src = static_cast<size_t>(graph.edge_index_flat[edge_idx]);
    const auto dst = static_cast<size_t>(graph.edge_index_flat[total_edges + edge_idx]);
    if (src < num_atoms && dst < num_atoms && src != dst) {
      mat[src * num_atoms + dst] = static_cast<real_t>(1.0);
      mat[dst * num_atoms + src] = static_cast<real_t>(1.0);
    }
  }

  auto eigenvalues = computeSymmetricEigenvalues(std::move(mat), num_atoms);
  const size_t result_k = std::min(k_eigenvalues, eigenvalues.size());
  eigenvalues.resize(result_k);
  return eigenvalues;
}

void GraphDescriptors::populateDescriptors(PeriodicGraphData &graph, size_t max_ring_size) {
  graph.cna_labels = computeCNADescriptor(graph);
  graph.coordination_desc = computeCoordinationEmbedding(graph);
  graph.ring_desc = computeRingStatisticsDescriptor(graph, max_ring_size);
}

} // namespace correlation::mlip
