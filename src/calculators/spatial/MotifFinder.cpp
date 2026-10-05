/**
 * @file MotifFinder.cpp
 * @brief Implementation of ring and structural motif search.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/spatial/MotifFinder.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <map>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <vector>

namespace correlation::calculators {

namespace {

// ---------------------------------------------------------------------------
// Per-thread scratch state.  One instance per TBB worker thread, kept alive
// for the lifetime of getAllShortestRings and reused across root iterations
// assigned to that thread.  Vectors are sized to num_nodes on construction
// and grow no further; they are reset selectively (only visited nodes) after
// each root iteration so reset cost is O(visited) rather than O(num_nodes).
// ---------------------------------------------------------------------------
struct BFSScratch {
  // BFS state
  std::vector<int> dist;
  std::vector<std::vector<correlation::core::AtomID>> parents;
  std::vector<size_t> q;
  std::vector<size_t> visited;
  std::vector<std::pair<size_t, size_t>> cross_edges;

  // King-ring and Franzblau check state (reused inside isRingValid)
  std::vector<int> dist_king;
  std::vector<uint32_t> path_counts;
  std::vector<size_t> q_king;
  std::vector<size_t> visited_king;

  // Output accumulated by this thread
  std::vector<std::vector<correlation::core::AtomID>> local_cycles;

  explicit BFSScratch(size_t n) : dist(n, -1), parents(n), dist_king(n, -1), path_counts(n, 0) {
    q.reserve(n);
    q_king.reserve(n);
    visited.reserve(n);
    visited_king.reserve(n);
    cross_edges.reserve(1024);
  }

  struct PathEndpoints {
    size_t start;
    size_t root;
  };

  struct KingBFSSettings {
    size_t start_node;
    int max_check_dist;
  };

  struct RootSearchSettings {
    size_t root;
    size_t max_size;
    correlation::analysis::RingType ring_type{correlation::analysis::RingType::King};
  };
};

// ---------------------------------------------------------------------------
// Reconstruct paths from a BFS node back to root (unchanged from original).
// ---------------------------------------------------------------------------
void getPaths(BFSScratch::PathEndpoints endpoints,
              const std::vector<std::vector<correlation::core::AtomID>> &parents,
              std::vector<correlation::core::AtomID> &current_path,
              std::vector<std::vector<correlation::core::AtomID>> &all_paths) {
  struct StackFrame {
    size_t node;
    size_t parent_idx;
  };

  std::vector<StackFrame> stack;
  stack.reserve(parents.empty() ? 16 : std::min<size_t>(parents.size(), 32));
  stack.push_back({
      .node = endpoints.start,
      .parent_idx = 0,
  });
  current_path.push_back(static_cast<correlation::core::AtomID>(endpoints.start));

  while (!stack.empty()) {
    auto &frame = stack.back();
    size_t const node = frame.node;

    if (node == endpoints.root) {
      all_paths.push_back(current_path);
      stack.pop_back();
      current_path.pop_back();
      continue;
    }

    const auto &node_parents = parents[node];
    if (frame.parent_idx < node_parents.size()) {
      size_t parent_node = node_parents[frame.parent_idx];
      frame.parent_idx++;

      stack.push_back({
          .node = parent_node,
          .parent_idx = 0,
      });
      current_path.push_back(static_cast<correlation::core::AtomID>(parent_node));
    } else {
      stack.pop_back();
      current_path.pop_back();
    }
  }
}

// ---------------------------------------------------------------------------
// King-ring and Franzblau primitive ring check.
// Uses dist_king / path_counts / q_king / visited_king from the caller-supplied
// scratch so no heap allocations occur here.
// ---------------------------------------------------------------------------
void runRingBFS(const correlation::core::NeighborGraph &graph, BFSScratch::KingBFSSettings settings,
                std::vector<int> &dist_king, std::vector<uint32_t> &path_counts,
                std::vector<size_t> &q_king, std::vector<size_t> &visited_nodes, bool count_paths) {
  q_king.clear();
  visited_nodes.clear();

  dist_king[settings.start_node] = 0;
  if (count_paths) {
    path_counts[settings.start_node] = 1;
  }
  q_king.push_back(settings.start_node);
  visited_nodes.push_back(settings.start_node);

  size_t q_head = 0;
  while (q_head < q_king.size()) {
    size_t const node = q_king[q_head++];
    int const current_dist = dist_king[node];
    uint32_t const current_paths = count_paths ? path_counts[node] : 0;

    if (current_dist >= settings.max_check_dist) {
      continue;
    }

    for (const auto &neighbor : graph.getNeighbors(node)) {
      size_t const neighbor_node = neighbor.index;
      if (dist_king[neighbor_node] == -1) {
        dist_king[neighbor_node] = current_dist + 1;
        if (count_paths) {
          path_counts[neighbor_node] = current_paths;
        }
        q_king.push_back(neighbor_node);
        visited_nodes.push_back(neighbor_node);
      } else if (count_paths && dist_king[neighbor_node] == current_dist + 1) {
        path_counts[neighbor_node] += current_paths;
      }
    }
  }
}

bool checkCyclePairDistances(const std::vector<correlation::core::AtomID> &cycle, size_t source_idx,
                             const std::vector<int> &dist_king,
                             const std::vector<uint32_t> &path_counts, bool count_paths) {
  size_t const size = cycle.size();
  for (size_t target_idx = 0; target_idx < size; ++target_idx) {
    if (source_idx == target_idx) {
      continue;
    }
    size_t const target_node = cycle[target_idx];
    size_t const diff =
        (target_idx > source_idx) ? (target_idx - source_idx) : (source_idx - target_idx);
    size_t const dist_in_cycle = std::min(diff, size - diff);
    int const d_g = dist_king[target_node];

    // 1. King chordless criterion: no shortcut path strictly shorter than cycle perimeter
    if (d_g != -1 && std::cmp_less(d_g, dist_in_cycle)) {
      return false;
    }

    // 2. Franzblau unique-geodesic condition: sub-path must be unique shortest path
    if (count_paths && std::cmp_equal(dist_in_cycle, d_g)) {
      uint32_t const expected_paths = (size % 2 == 0 && dist_in_cycle == size / 2) ? 2U : 1U;
      if (path_counts[target_node] != expected_paths) {
        return false;
      }
    }
  }
  return true;
}

bool isRingValid(const correlation::core::NeighborGraph &graph,
                 const std::vector<correlation::core::AtomID> &cycle, std::vector<int> &dist_king,
                 std::vector<uint32_t> &path_counts, std::vector<size_t> &q_king,
                 std::vector<size_t> &visited_king, correlation::analysis::RingType ring_type) {
  size_t const size = cycle.size();
  if (size < 3) {
    return false;
  }

  bool const count_paths = (ring_type == correlation::analysis::RingType::Franzblau);

  for (size_t i = 0; i < size; ++i) {
    size_t const start_node = cycle[i];

    runRingBFS(graph,
               {
                   .start_node = start_node,
                   .max_check_dist = static_cast<int>(size / 2),
               },
               dist_king, path_counts, q_king, visited_king, count_paths);

    bool const valid = checkCyclePairDistances(cycle, i, dist_king, path_counts, count_paths);

    for (size_t const visited_node : visited_king) {
      dist_king[visited_node] = -1;
      if (count_paths) {
        path_counts[visited_node] = 0;
      }
    }

    if (!valid) {
      return false;
    }
  }

  return true;
}

// ---------------------------------------------------------------------------
// BFS from a single root.  All state lives in `sc`; found rings are appended
// to `sc.local_cycles`.
//
// Safety note: the `v >= root` guard in the edge-exploration loop ensures each
// ring is discoverable ONLY from its minimum-index node.  Different threads
// therefore cannot produce the same ring — no shared deduplication set is
// needed during the parallel section.  A final sort+unique after the
// parallel_for handles any remaining orientation duplicates.
// ---------------------------------------------------------------------------
void processNeighbor(size_t curr_node, size_t neighbor_node,
                     BFSScratch::RootSearchSettings settings, BFSScratch &bsc) {
  if (bsc.dist[neighbor_node] == -1) {
    bsc.dist[neighbor_node] = bsc.dist[curr_node] + 1;
    bsc.parents[neighbor_node].push_back(static_cast<correlation::core::AtomID>(curr_node));
    bsc.visited.push_back(neighbor_node);
    if (2 * bsc.dist[neighbor_node] + 1 <= static_cast<int>(settings.max_size)) {
      bsc.q.push_back(neighbor_node);
    }
  } else if (bsc.dist[neighbor_node] == bsc.dist[curr_node]) {
    if (curr_node < neighbor_node) {
      bsc.cross_edges.emplace_back(curr_node, neighbor_node);
    }
  } else if (bsc.dist[neighbor_node] == bsc.dist[curr_node] + 1) {
    if (std::find(bsc.parents[neighbor_node].begin(), bsc.parents[neighbor_node].end(),
                  static_cast<correlation::core::AtomID>(curr_node)) ==
        bsc.parents[neighbor_node].end()) {
      bsc.parents[neighbor_node].push_back(static_cast<correlation::core::AtomID>(curr_node));
      bsc.cross_edges.emplace_back(curr_node, neighbor_node);
    }
  }
}

void findCrossEdges(const correlation::core::NeighborGraph &graph,
                    BFSScratch::RootSearchSettings settings, BFSScratch &bsc) {
  size_t q_head = 0;

  while (q_head < bsc.q.size()) {
    size_t const curr_node = bsc.q[q_head++];

    for (const auto &neighbor : graph.getNeighbors(curr_node)) {
      size_t const neighbor_node = neighbor.index;

      if (neighbor_node == curr_node) {
        continue;
      }
      if (neighbor_node < settings.root) {
        continue;
      }

      processNeighbor(curr_node, neighbor_node, settings, bsc);
    }
  }
}

bool pathsIntersect(const std::vector<correlation::core::AtomID> &path_u,
                    const std::vector<correlation::core::AtomID> &path_v) {
  for (size_t i = 0; i < path_u.size() - 1; ++i) {
    for (size_t j = 0; j < path_v.size() - 1; ++j) {
      if (path_u[i] == path_v[j]) {
        return true;
      }
    }
  }
  return false;
}

void processCrossEdge(const correlation::core::NeighborGraph &graph,
                      const std::pair<size_t, size_t> &edge,
                      BFSScratch::RootSearchSettings settings, BFSScratch &bsc) {
  size_t const first_node = edge.first;
  size_t const second_node = edge.second;

  if (bsc.dist[first_node] + bsc.dist[second_node] + 1 > static_cast<int>(settings.max_size)) {
    return;
  }

  std::vector<std::vector<correlation::core::AtomID>> paths_u;
  std::vector<std::vector<correlation::core::AtomID>> paths_v;
  std::vector<correlation::core::AtomID> cur_u;
  std::vector<correlation::core::AtomID> cur_v;
  cur_u.reserve(settings.max_size);
  cur_v.reserve(settings.max_size);

  getPaths(
      {
          .start = first_node,
          .root = settings.root,
      },
      bsc.parents, cur_u, paths_u);
  getPaths(
      {
          .start = second_node,
          .root = settings.root,
      },
      bsc.parents, cur_v, paths_v);

  for (const auto &path_u : paths_u) {
    for (const auto &path_v : paths_v) {
      if (pathsIntersect(path_u, path_v)) {
        continue;
      }

      std::vector<correlation::core::AtomID> cycle;
      cycle.reserve(path_u.size() + path_v.size() - 1);
      cycle.push_back(static_cast<correlation::core::AtomID>(settings.root));
      for (int i = static_cast<int>(path_u.size()) - 2; i >= 0; --i) {
        cycle.push_back(path_u[i]);
      }
      for (size_t i = 0; i < path_v.size() - 1; ++i) {
        cycle.push_back(path_v[i]);
      }

      if (cycle.size() < 3 || cycle.size() > settings.max_size) {
        continue;
      }

      // Normalise direction
      if (cycle[1] > cycle.back()) {
        std::reverse(cycle.begin() + 1, cycle.end());
      }

      if (isRingValid(graph, cycle, bsc.dist_king, bsc.path_counts, bsc.q_king, bsc.visited_king,
                      settings.ring_type)) {
        bsc.local_cycles.push_back(std::move(cycle));
      }
    }
  }
}

void processRoot(const correlation::core::NeighborGraph &graph,
                 BFSScratch::RootSearchSettings settings, BFSScratch &bsc) {
  bsc.visited.clear();
  bsc.cross_edges.clear();
  bsc.q.clear();

  bsc.dist[settings.root] = 0;
  bsc.q.push_back(settings.root);
  bsc.visited.push_back(settings.root);

  findCrossEdges(graph, settings, bsc);

  for (const auto &edge : bsc.cross_edges) {
    processCrossEdge(graph, edge, settings, bsc);
  }

  // Selective reset: only touch the nodes we visited
  for (size_t const visited_node : bsc.visited) {
    bsc.dist[visited_node] = -1;
    bsc.parents[visited_node].clear();
  }
}

// ---------------------------------------------------------------------------
// Main ring-finding function — parallel over roots.
// ---------------------------------------------------------------------------
std::vector<std::vector<correlation::core::AtomID>> getAllShortestRings(
    const correlation::core::NeighborGraph &graph, size_t max_size,
    correlation::analysis::RingType ring_type = correlation::analysis::RingType::King) {
  std::vector<std::vector<correlation::core::AtomID>> all_cycles;
  if (max_size < 3) {
    return all_cycles;
  }

  const size_t num_nodes = graph.nodeCount();
  if (num_nodes == 0) {
    return all_cycles;
  }

  // Each TBB thread owns one BFSScratch, sized at construction and reused
  // across all root iterations assigned to that thread.
  tbb::enumerable_thread_specific<BFSScratch> ets([num_nodes] { return BFSScratch(num_nodes); });

  // Grain size 16: balances TBB overhead (~µs per task) against load
  // imbalance (root-0 does far more work than root-N-1).
  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, num_nodes, /*grain=*/16),
      [&](const tbb::blocked_range<size_t> &range) {
        BFSScratch &bsc = ets.local();
        for (size_t root = range.begin(); root != range.end(); ++root) {
          processRoot(graph,
                      {
                          .root = root,
                          .max_size = max_size,
                          .ring_type = ring_type,
                      },
                      bsc);
        }
      },
      tbb::auto_partitioner{});

  // Serial merge of all per-thread cycle lists
  for (auto &bsc : ets) {
    all_cycles.insert(all_cycles.end(), std::make_move_iterator(bsc.local_cycles.begin()),
                      std::make_move_iterator(bsc.local_cycles.end()));
  }

  // Final deduplication
  std::ranges::sort(all_cycles);
  auto const [first, last] = std::ranges::unique(all_cycles);
  all_cycles.erase(first, last);

  return all_cycles;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------
std::map<int, size_t> MotifFinder::findRings(const correlation::core::NeighborGraph &graph,
                                             size_t max_size) {
  return findRings(graph, max_size, RingType::King);
}

std::map<int, size_t> MotifFinder::findFranzblauRings(const correlation::core::NeighborGraph &graph,
                                                      size_t max_size) {
  return findRings(graph, max_size, RingType::Franzblau);
}

std::map<int, size_t> MotifFinder::findRings(const correlation::core::NeighborGraph &graph,
                                             size_t max_size, RingType ring_type) {
  auto all_cycles = getAllShortestRings(graph, max_size, ring_type);
  std::map<int, size_t> ring_counts;
  for (const auto &cycle : all_cycles) {
    ring_counts[static_cast<int>(cycle.size())]++;
  }
  return ring_counts;
}

std::vector<std::vector<correlation::core::AtomID>>
MotifFinder::extractCycles(const correlation::core::NeighborGraph &graph, size_t target_size) {
  return extractCycles(graph, target_size, RingType::King);
}

std::vector<std::vector<correlation::core::AtomID>>
MotifFinder::extractFranzblauCycles(const correlation::core::NeighborGraph &graph,
                                    size_t target_size) {
  return extractCycles(graph, target_size, RingType::Franzblau);
}

std::vector<std::vector<correlation::core::AtomID>>
MotifFinder::extractCycles(const correlation::core::NeighborGraph &graph, size_t target_size,
                           RingType ring_type) {
  auto all_cycles = getAllShortestRings(graph, target_size, ring_type);
  std::vector<std::vector<correlation::core::AtomID>> exact_cycles;
  exact_cycles.reserve(all_cycles.size());
  for (auto &cycle : all_cycles) {
    if (cycle.size() == target_size) {
      exact_cycles.push_back(std::move(cycle));
    }
  }
  return exact_cycles;
}

std::vector<std::vector<correlation::core::AtomID>>
MotifFinder::extractAllCycles(const correlation::core::NeighborGraph &graph, size_t max_size,
                              RingType ring_type) {
  return getAllShortestRings(graph, max_size, ring_type);
}

correlation::core::NeighborGraph
MotifFinder::buildBridgedGraph(const correlation::core::NeighborGraph &graph,
                               const correlation::core::Cell &cell, std::string_view former_element,
                               std::string_view bridging_element) {
  const auto &atoms = cell.atoms();
  size_t const num_nodes = graph.nodeCount();
  correlation::core::NeighborGraph bridged_graph(num_nodes);

  if (former_element.empty() || bridging_element.empty() || num_nodes != atoms.size()) {
    return bridged_graph;
  }

  std::vector<bool> connected(num_nodes, false);

  for (size_t i = 0; i < num_nodes; ++i) {
    if (atoms[i].element().symbol != former_element) {
      continue;
    }

    connected.assign(num_nodes, false);
    connected[i] = true;

    // 1. Second-hop neighbors through bridging element
    for (const auto &bridge_nbr : graph.getNeighbors(i)) {
      size_t const bridge_idx = bridge_nbr.index;
      if (bridge_idx >= num_nodes || atoms[bridge_idx].element().symbol != bridging_element) {
        continue;
      }

      for (const auto &second_nbr : graph.getNeighbors(bridge_idx)) {
        size_t const second_idx = second_nbr.index;
        if (second_idx < num_nodes && atoms[second_idx].element().symbol == former_element &&
            !connected[second_idx]) {
          connected[second_idx] = true;
          auto const disp = cell.minimumImage(atoms[second_idx].position() - atoms[i].position());
          real_t const dist =
              std::sqrt(disp.x() * disp.x() + disp.y() * disp.y() + disp.z() * disp.z());
          bridged_graph.addDirectedEdge(i, second_idx, dist, disp);
        }
      }
    }

    // 2. Direct homopolar former-former bonds
    for (const auto &nbr : graph.getNeighbors(i)) {
      size_t const nbr_idx = nbr.index;
      if (nbr_idx < num_nodes && atoms[nbr_idx].element().symbol == former_element &&
          !connected[nbr_idx]) {
        connected[nbr_idx] = true;
        bridged_graph.addDirectedEdge(i, nbr_idx, nbr.distance, nbr.r_ij);
      }
    }
  }

  return bridged_graph;
}

std::vector<std::vector<correlation::core::AtomID>> MotifFinder::filterAlternatingCycles(
    const std::vector<std::vector<correlation::core::AtomID>> &cycles,
    const correlation::core::Cell &cell, std::string_view element_a, std::string_view element_b) {
  const auto &atoms = cell.atoms();
  std::vector<std::vector<correlation::core::AtomID>> filtered;

  for (const auto &cycle : cycles) {
    size_t const len = cycle.size();
    if (len < 4 || len % 2 != 0) {
      continue;
    }

    std::string_view const first_sym = atoms[cycle[0]].element().symbol;
    std::string_view expected_even;
    std::string_view expected_odd;

    if (first_sym == element_a) {
      expected_even = element_a;
      expected_odd = element_b;
    } else if (first_sym == element_b) {
      expected_even = element_b;
      expected_odd = element_a;
    } else {
      continue;
    }

    bool is_alternating = true;
    for (size_t k = 0; k < len; ++k) {
      std::string_view const sym = atoms[cycle[k]].element().symbol;
      std::string_view const expected = (k % 2 == 0) ? expected_even : expected_odd;
      if (sym != expected) {
        is_alternating = false;
        break;
      }
    }

    if (is_alternating) {
      filtered.push_back(cycle);
    }
  }

  return filtered;
}

} // namespace correlation::calculators
