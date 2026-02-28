// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/MotifFinder.hpp"
#include <algorithm>
#include <set>

namespace calculators {

namespace {

// Helper to normalize a cycle representation to easily avoid duplicates
// E.g. cycle 3-1-2 is equivalent to 1-2-3 (forward) and 3-2-1 is equivalent to
// 1-3-2 A unique ring representation starts with the smallest node ID, and the
// second node is the smaller of its two neighbors.
std::vector<AtomID> normalizeCycle(const std::vector<AtomID> &cycle) {
  if (cycle.empty())
    return cycle;

  auto min_it = std::min_element(cycle.begin(), cycle.end());
  size_t min_idx = std::distance(cycle.begin(), min_it);

  std::vector<AtomID> forward(cycle.size());
  std::vector<AtomID> backward(cycle.size());

  for (size_t i = 0; i < cycle.size(); ++i) {
    forward[i] = cycle[(min_idx + i) % cycle.size()];
    backward[i] = cycle[(min_idx + cycle.size() - i) % cycle.size()];
  }

  if (forward[1] < backward[1]) {
    return forward;
  }
  return backward;
}

void dfs(size_t current_node, size_t start_node, size_t depth, size_t max_depth,
         const NeighborGraph &graph, std::vector<AtomID> &current_path,
         std::vector<bool> &visited,
         std::set<std::vector<AtomID>> &found_cycles) {

  if (depth > max_depth)
    return;

  visited[current_node] = true;
  current_path.push_back(static_cast<AtomID>(current_node));

  for (const auto &neighbor : graph.getNeighbors(current_node)) {
    size_t next_node = neighbor.index;

    // Don't go right back to the node we just came from
    if (current_path.size() > 1 &&
        next_node == current_path[current_path.size() - 2]) {
      continue;
    }

    if (next_node == start_node) {
      // Found a cycle! (only if depth >= 3, which is implicitly true because we
      // don't go backwards)
      if (depth >= 3) {
        found_cycles.insert(normalizeCycle(current_path));
      }
    } else if (!visited[next_node] && depth < max_depth) {
      dfs(next_node, start_node, depth + 1, max_depth, graph, current_path,
          visited, found_cycles);
    }
  }

  // Backtrack
  current_path.pop_back();
  visited[current_node] = false;
}

} // anonymous namespace

std::map<int, size_t> MotifFinder::findRings(const NeighborGraph &graph,
                                             size_t max_size) {
  std::map<int, size_t> ring_counts;

  // We cannot reliably search below size 3 since size 1 is a self loop and 2 is
  // just an edge.
  if (max_size < 3)
    return ring_counts;

  std::set<std::vector<AtomID>> all_unique_cycles;

  std::vector<bool> visited(graph.nodeCount(), false);
  std::vector<AtomID> current_path;

  for (size_t i = 0; i < graph.nodeCount(); ++i) {
    // Start DFS from each node
    dfs(i, i, 1, max_size, graph, current_path, visited, all_unique_cycles);
  }

  for (const auto &cycle : all_unique_cycles) {
    ring_counts[static_cast<int>(cycle.size())]++;
  }

  return ring_counts;
}

std::vector<std::vector<AtomID>>
MotifFinder::extractCycles(const NeighborGraph &graph, size_t target_size) {
  std::vector<std::vector<AtomID>> exact_cycles;

  if (target_size < 3)
    return exact_cycles;

  std::set<std::vector<AtomID>> all_unique_cycles;
  std::vector<bool> visited(graph.nodeCount(), false);
  std::vector<AtomID> current_path;

  for (size_t i = 0; i < graph.nodeCount(); ++i) {
    dfs(i, i, 1, target_size, graph, current_path, visited, all_unique_cycles);
  }

  for (const auto &cycle : all_unique_cycles) {
    if (cycle.size() == target_size) {
      exact_cycles.push_back(cycle);
    }
  }

  return exact_cycles;
}

} // namespace calculators
