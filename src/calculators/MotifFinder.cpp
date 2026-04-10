/**
 * @file MotifFinder.cpp
 * @brief Implementation of ring and structural motif search.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "calculators/MotifFinder.hpp"

#include <algorithm>
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

  // King-ring check state (reused inside isKingRing)
  std::vector<int> dist_king;
  std::vector<size_t> q_king;

  // Output accumulated by this thread
  std::vector<std::vector<correlation::core::AtomID>> local_cycles;

  explicit BFSScratch(size_t n) : dist(n, -1), parents(n), dist_king(n, -1) {
    q.reserve(n);
    q_king.reserve(n);
    visited.reserve(n);
    cross_edges.reserve(1024);
  }
};

// ---------------------------------------------------------------------------
// Reconstruct paths from a BFS node back to root (unchanged from original).
// ---------------------------------------------------------------------------
void get_paths(
    size_t node, size_t root,
    const std::vector<std::vector<correlation::core::AtomID>> &parents,
    std::vector<correlation::core::AtomID> &current_path,
    std::vector<std::vector<correlation::core::AtomID>> &all_paths) {
  if (node == root) {
    current_path.push_back(static_cast<correlation::core::AtomID>(root));
    all_paths.push_back(current_path);
    current_path.pop_back();
    return;
  }
  current_path.push_back(static_cast<correlation::core::AtomID>(node));
  for (correlation::core::AtomID p : parents[node]) {
    get_paths(p, root, parents, current_path, all_paths);
  }
  current_path.pop_back();
}

// ---------------------------------------------------------------------------
// King-ring check.  Uses dist_king / q_king from the caller-supplied scratch
// so no heap allocations occur here (unchanged algorithm, scalar state only).
// ---------------------------------------------------------------------------
bool isKingRing(const correlation::core::NeighborGraph &graph,
                const std::vector<correlation::core::AtomID> &cycle,
                std::vector<int> &dist_king, std::vector<size_t> &q_king) {
  int n = cycle.size();
  if (n < 3)
    return false;

  std::vector<size_t> visited_nodes;
  visited_nodes.reserve(n * 10);

  for (int i = 0; i < n; ++i) {
    size_t start_node = cycle[i];

    q_king.clear();
    visited_nodes.clear();

    dist_king[start_node] = 0;
    q_king.push_back(start_node);
    visited_nodes.push_back(start_node);

    size_t q_head = 0;
    int max_check_dist = n / 2;

    while (q_head < q_king.size()) {
      size_t u = q_king[q_head++];
      int current_dist = dist_king[u];

      if (current_dist >= max_check_dist)
        continue;

      for (const auto &neighbor : graph.getNeighbors(u)) {
        size_t v = neighbor.index;
        if (dist_king[v] == -1) {
          dist_king[v] = current_dist + 1;
          q_king.push_back(v);
          visited_nodes.push_back(v);
        }
      }
    }

    bool is_king = true;
    for (int j = 0; j < n; ++j) {
      if (i == j)
        continue;
      size_t target_node = cycle[j];
      int dist_in_cycle = std::min(std::abs(j - i), n - std::abs(j - i));
      if (dist_king[target_node] != -1 &&
          dist_king[target_node] < dist_in_cycle) {
        is_king = false;
        break;
      }
    }

    for (size_t v : visited_nodes)
      dist_king[v] = -1;

    if (!is_king)
      return false;
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
void process_root(const correlation::core::NeighborGraph &graph, size_t root,
                  size_t max_size, BFSScratch &sc) {
  sc.visited.clear();
  sc.cross_edges.clear();
  sc.q.clear();

  sc.dist[root] = 0;
  sc.q.push_back(root);
  sc.visited.push_back(root);

  size_t q_head = 0;

  while (q_head < sc.q.size()) {
    size_t u = sc.q[q_head++];

    for (const auto &neighbor : graph.getNeighbors(u)) {
      size_t v = neighbor.index;

      if (v == u)
        continue;
      if (v < root)
        continue;

      if (sc.dist[v] == -1) {
        sc.dist[v] = sc.dist[u] + 1;
        sc.parents[v].push_back(static_cast<correlation::core::AtomID>(u));
        sc.visited.push_back(v);
        if (2 * sc.dist[v] + 1 <= static_cast<int>(max_size))
          sc.q.push_back(v);
      } else if (sc.dist[v] == sc.dist[u]) {
        if (u < v)
          sc.cross_edges.push_back({u, v});
      } else if (sc.dist[v] == sc.dist[u] + 1) {
        if (std::find(sc.parents[v].begin(), sc.parents[v].end(),
                      static_cast<correlation::core::AtomID>(u)) ==
            sc.parents[v].end()) {
          sc.parents[v].push_back(static_cast<correlation::core::AtomID>(u));
          sc.cross_edges.push_back({u, v});
        }
      }
    }
  }

  for (const auto &edge : sc.cross_edges) {
    size_t u = edge.first;
    size_t v = edge.second;

    if (sc.dist[u] + sc.dist[v] + 1 > static_cast<int>(max_size))
      continue;

    std::vector<std::vector<correlation::core::AtomID>> paths_u, paths_v;
    std::vector<correlation::core::AtomID> cur_u, cur_v;
    cur_u.reserve(max_size);
    cur_v.reserve(max_size);

    get_paths(u, root, sc.parents, cur_u, paths_u);
    get_paths(v, root, sc.parents, cur_v, paths_v);

    for (const auto &pu : paths_u) {
      for (const auto &pv : paths_v) {
        bool intersect = false;
        for (size_t i = 0; i < pu.size() - 1 && !intersect; ++i)
          for (size_t j = 0; j < pv.size() - 1; ++j)
            if (pu[i] == pv[j]) {
              intersect = true;
              break;
            }

        if (intersect)
          continue;

        std::vector<correlation::core::AtomID> cycle;
        cycle.reserve(pu.size() + pv.size() - 1);
        cycle.push_back(static_cast<correlation::core::AtomID>(root));
        for (int i = static_cast<int>(pu.size()) - 2; i >= 0; --i)
          cycle.push_back(pu[i]);
        for (size_t i = 0; i < pv.size() - 1; ++i)
          cycle.push_back(pv[i]);

        if (cycle.size() < 3 || cycle.size() > max_size)
          continue;

        // Normalise direction
        if (cycle[1] > cycle.back())
          std::reverse(cycle.begin() + 1, cycle.end());

        if (isKingRing(graph, cycle, sc.dist_king, sc.q_king))
          sc.local_cycles.push_back(std::move(cycle));
      }
    }
  }

  // Selective reset: only touch the nodes we visited
  for (size_t v : sc.visited) {
    sc.dist[v] = -1;
    sc.parents[v].clear();
  }
}

// ---------------------------------------------------------------------------
// Main ring-finding function — now parallel over roots.
// ---------------------------------------------------------------------------
std::vector<std::vector<correlation::core::AtomID>>
getAllShortestRings(const correlation::core::NeighborGraph &graph,
                    size_t max_size) {
  std::vector<std::vector<correlation::core::AtomID>> all_cycles;
  if (max_size < 3)
    return all_cycles;

  const size_t num_nodes = graph.nodeCount();
  if (num_nodes == 0)
    return all_cycles;

  // Each TBB thread owns one BFSScratch, sized at construction and reused
  // across all root iterations assigned to that thread.
  tbb::enumerable_thread_specific<BFSScratch> ets(
      [num_nodes] { return BFSScratch(num_nodes); });

  // Grain size 16: balances TBB overhead (~µs per task) against load
  // imbalance (root-0 does far more work than root-N-1).
  // auto_partitioner further subdivides at runtime if needed.
  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, num_nodes, /*grain=*/16),
      [&](const tbb::blocked_range<size_t> &range) {
        BFSScratch &sc = ets.local();
        for (size_t root = range.begin(); root != range.end(); ++root)
          process_root(graph, root, max_size, sc);
      },
      tbb::auto_partitioner{});

  // Serial merge of all per-thread cycle lists
  for (auto &sc : ets) {
    all_cycles.insert(all_cycles.end(),
                      std::make_move_iterator(sc.local_cycles.begin()),
                      std::make_move_iterator(sc.local_cycles.end()));
  }

  // Final deduplication (oriented cycles may still have orientation variants
  // produced by different paths within the same root's BFS)
  std::sort(all_cycles.begin(), all_cycles.end());
  all_cycles.erase(std::unique(all_cycles.begin(), all_cycles.end()),
                   all_cycles.end());

  return all_cycles;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Public API — unchanged
// ---------------------------------------------------------------------------
std::map<int, size_t>
MotifFinder::findRings(const correlation::core::NeighborGraph &graph,
                       size_t max_size) {
  auto all_cycles = getAllShortestRings(graph, max_size);
  std::map<int, size_t> ring_counts;
  for (const auto &cycle : all_cycles)
    ring_counts[static_cast<int>(cycle.size())]++;
  return ring_counts;
}

std::vector<std::vector<correlation::core::AtomID>>
MotifFinder::extractCycles(const correlation::core::NeighborGraph &graph,
                           size_t target_size) {
  auto all_cycles = getAllShortestRings(graph, target_size);
  std::vector<std::vector<correlation::core::AtomID>> exact_cycles;
  exact_cycles.reserve(all_cycles.size());
  for (auto &cycle : all_cycles)
    if (cycle.size() == target_size)
      exact_cycles.push_back(std::move(cycle));
  return exact_cycles;
}

} // namespace correlation::calculators
