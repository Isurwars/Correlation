/**
 * @file MotifFinder.cpp
 * @brief Implementation of ring and structural motif search.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/MotifFinder.hpp"

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
  };
};

// ---------------------------------------------------------------------------
// Reconstruct paths from a BFS node back to root (unchanged from original).
// ---------------------------------------------------------------------------
void get_paths(BFSScratch::PathEndpoints endpoints, const std::vector<std::vector<correlation::core::AtomID>> &parents,
               std::vector<correlation::core::AtomID> &current_path,
               std::vector<std::vector<correlation::core::AtomID>> &all_paths) {
  struct StackFrame {
    size_t node;
    size_t parent_idx;
  };

  std::vector<StackFrame> stack;
  stack.push_back({.node = endpoints.start, .parent_idx = 0});
  current_path.push_back(static_cast<correlation::core::AtomID>(endpoints.start));

  while (!stack.empty()) {
    auto &frame = stack.back();
    size_t node = frame.node;

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

      stack.push_back({.node = parent_node, .parent_idx = 0});
      current_path.push_back(static_cast<correlation::core::AtomID>(parent_node));
    } else {
      stack.pop_back();
      current_path.pop_back();
    }
  }
}

// ---------------------------------------------------------------------------
// King-ring check.  Uses dist_king / q_king from the caller-supplied scratch
// so no heap allocations occur here (unchanged algorithm, scalar state only).
// ---------------------------------------------------------------------------
void runKingBFS(const correlation::core::NeighborGraph &graph, BFSScratch::KingBFSSettings settings,
                std::vector<int> &dist_king, std::vector<size_t> &q_king, std::vector<size_t> &visited_nodes) {
  q_king.clear();
  visited_nodes.clear();

  dist_king[settings.start_node] = 0;
  q_king.push_back(settings.start_node);
  visited_nodes.push_back(settings.start_node);

  size_t q_head = 0;
  while (q_head < q_king.size()) {
    size_t const node = q_king[q_head++];
    int const current_dist = dist_king[node];

    if (current_dist >= settings.max_check_dist) {
      continue;
    }

    for (const auto &neighbor : graph.getNeighbors(node)) {
      size_t const neighbor_node = neighbor.index;
      if (dist_king[neighbor_node] == -1) {
        dist_king[neighbor_node] = current_dist + 1;
        q_king.push_back(neighbor_node);
        visited_nodes.push_back(neighbor_node);
      }
    }
  }
}

bool isKingRing(const correlation::core::NeighborGraph &graph, const std::vector<correlation::core::AtomID> &cycle,
                std::vector<int> &dist_king, std::vector<size_t> &q_king) {
  size_t const size = cycle.size();
  if (size < 3) {
    return false;
  }

  std::vector<size_t> visited_nodes;
  visited_nodes.reserve(size * 10);

  for (size_t i = 0; i < size; ++i) {
    size_t const start_node = cycle[i];

    runKingBFS(graph, {.start_node = start_node, .max_check_dist = static_cast<int>(size / 2)}, dist_king, q_king,
               visited_nodes);

    bool is_king = true;
    for (size_t j = 0; j < size; ++j) {
      if (i == j) {
        continue;
      }
      size_t const target_node = cycle[j];
      size_t const diff = (j > i) ? (j - i) : (i - j);
      size_t const dist_in_cycle = std::min(diff, size - diff);
      if (dist_king[target_node] != -1 && static_cast<size_t>(dist_king[target_node]) < dist_in_cycle) {
        is_king = false;
        break;
      }
    }

    for (size_t const visited_node : visited_nodes) {
      dist_king[visited_node] = -1;
    }

    if (!is_king) {
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
void processNeighbor(size_t curr_node, size_t neighbor_node, BFSScratch::RootSearchSettings settings, BFSScratch &bsc) {
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
                  static_cast<correlation::core::AtomID>(curr_node)) == bsc.parents[neighbor_node].end()) {
      bsc.parents[neighbor_node].push_back(static_cast<correlation::core::AtomID>(curr_node));
      bsc.cross_edges.emplace_back(curr_node, neighbor_node);
    }
  }
}

void findCrossEdges(const correlation::core::NeighborGraph &graph, BFSScratch::RootSearchSettings settings,
                    BFSScratch &bsc) {
  size_t q_head = 0;

  while (q_head < bsc.q.size()) {
    size_t curr_node = bsc.q[q_head++];

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

void processCrossEdge(const correlation::core::NeighborGraph &graph, const std::pair<size_t, size_t> &edge,
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

  get_paths({.start = first_node, .root = settings.root}, bsc.parents, cur_u, paths_u);
  get_paths({.start = second_node, .root = settings.root}, bsc.parents, cur_v, paths_v);

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

      if (isKingRing(graph, cycle, bsc.dist_king, bsc.q_king)) {
        bsc.local_cycles.push_back(std::move(cycle));
      }
    }
  }
}

void process_root(const correlation::core::NeighborGraph &graph, BFSScratch::RootSearchSettings settings,
                  BFSScratch &bsc) {
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
// Main ring-finding function — now parallel over roots.
// ---------------------------------------------------------------------------
std::vector<std::vector<correlation::core::AtomID>> getAllShortestRings(const correlation::core::NeighborGraph &graph,
                                                                        size_t max_size) {
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
  // auto_partitioner further subdivides at runtime if needed.
  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, num_nodes, /*grain=*/16),
      [&](const tbb::blocked_range<size_t> &range) {
        BFSScratch &bsc = ets.local();
        for (size_t root = range.begin(); root != range.end(); ++root) {
          process_root(graph, {.root = root, .max_size = max_size}, bsc);
        }
      },
      tbb::auto_partitioner{});

  // Serial merge of all per-thread cycle lists
  for (auto &bsc : ets) {
    all_cycles.insert(all_cycles.end(), std::make_move_iterator(bsc.local_cycles.begin()),
                      std::make_move_iterator(bsc.local_cycles.end()));
  }

  // Final deduplication (oriented cycles may still have orientation variants
  // produced by different paths within the same root's BFS)
  std::sort(all_cycles.begin(), all_cycles.end());
  auto last = std::unique(all_cycles.begin(), all_cycles.end());
  all_cycles.erase(last, all_cycles.end());

  return all_cycles;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Public API — unchanged
// ---------------------------------------------------------------------------
std::map<int, size_t> MotifFinder::findRings(const correlation::core::NeighborGraph &graph, size_t max_size) {
  auto all_cycles = getAllShortestRings(graph, max_size);
  std::map<int, size_t> ring_counts;
  for (const auto &cycle : all_cycles) {
    ring_counts[static_cast<int>(cycle.size())]++;
  }
  return ring_counts;
}

std::vector<std::vector<correlation::core::AtomID>>
MotifFinder::extractCycles(const correlation::core::NeighborGraph &graph, size_t target_size) {
  auto all_cycles = getAllShortestRings(graph, target_size);
  std::vector<std::vector<correlation::core::AtomID>> exact_cycles;
  exact_cycles.reserve(all_cycles.size());
  for (auto &cycle : all_cycles) {
    if (cycle.size() == target_size) {
      exact_cycles.push_back(std::move(cycle));
    }
  }
  return exact_cycles;
}

} // namespace correlation::calculators
