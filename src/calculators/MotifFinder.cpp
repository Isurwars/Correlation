// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/MotifFinder.hpp"
#include <algorithm>
#include <set>
#include <vector>

namespace calculators {

namespace {

/**
 * @brief A helper to recursively trace paths back from a target node to the
 * root.
 *
 * Used during the Breadth-First Search (BFS) after a cyclic cross-edge is found
 * to reconstruct the two halves of the ring.
 */
void get_paths(size_t node, size_t root,
               const std::vector<std::vector<AtomID>> &parents,
               std::vector<AtomID> &current_path,
               std::vector<std::vector<AtomID>> &all_paths) {
  if (node == root) {
    current_path.push_back(static_cast<AtomID>(root));
    all_paths.push_back(current_path);
    current_path.pop_back();
    return;
  }
  current_path.push_back(static_cast<AtomID>(node));
  for (AtomID p : parents[node]) {
    get_paths(p, root, parents, current_path, all_paths);
  }
  current_path.pop_back();
}

/**
 * @brief Checks if a given cycle is a King's ring (shortest-path ring) in the
 * graph.
 *
 * A cycle is a King's ring if the shortest path distance in the entire graph
 * between any two nodes in the cycle is equal to their shortest distance along
 * the cycle. This prevents identifying large macro-cycles that wrap around
 * smaller internal structures.
 */
bool isKingRing(const NeighborGraph &graph, const std::vector<AtomID> &cycle,
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

      if (current_dist >= max_check_dist) {
        continue;
      }

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

    for (size_t v : visited_nodes) {
      dist_king[v] = -1;
    }

    if (!is_king)
      return false;
  }

  return true;
}

/**
 * @brief Breadth-First Search algorithm to find all shortest chordless cycles
 * (rings).
 *
 * The algorithm iterates through each node as a root, performing a level-order
 * traversal. When the BFS encounters an edge connecting two nodes at the same
 * depth (or depth + 1) without an immediate parent relationship, it identifies
 * a "cross-edge" completing a ring.
 *
 * @param graph The atomic neighbor graph.
 * @param max_size The maximum chordless ring size.
 */
std::vector<std::vector<AtomID>> getAllShortestRings(const NeighborGraph &graph,
                                                     size_t max_size) {
  std::vector<std::vector<AtomID>> all_cycles;
  if (max_size < 3)
    return all_cycles;

  std::set<std::vector<AtomID>> unique_ring_nodes;

  size_t num_nodes = graph.nodeCount();
  // We hoist memory allocations out of the search loops to drastically improve
  // performance.
  std::vector<int> dist(num_nodes, -1);
  std::vector<std::vector<AtomID>> parents(num_nodes);
  std::vector<size_t> visited;
  visited.reserve(num_nodes); // Prevent reallocation

  std::vector<int> dist_king(num_nodes, -1);
  std::vector<size_t> q_king;
  q_king.reserve(num_nodes);

  std::vector<std::pair<size_t, size_t>> cross_edges;
  cross_edges.reserve(1024);

  // Use a vector as a queue to avoid reallocations
  std::vector<size_t> q;
  q.reserve(num_nodes);

  for (size_t root = 0; root < num_nodes; ++root) {
    visited.clear();
    cross_edges.clear();
    q.clear();

    dist[root] = 0;
    q.push_back(root);
    visited.push_back(root);

    size_t q_head = 0;

    while (q_head < q.size()) {
      size_t u = q[q_head++];

      for (const auto &neighbor : graph.getNeighbors(u)) {
        size_t v = neighbor.index;

        if (v == u)
          continue; // Skip self-loops

        if (v < root)
          continue;

        if (dist[v] == -1) {
          dist[v] = dist[u] + 1;
          parents[v].push_back(static_cast<AtomID>(u));
          visited.push_back(v);
          if (2 * dist[v] + 1 <= static_cast<int>(max_size)) {
            q.push_back(v);
          }
        } else if (dist[v] == dist[u]) {
          if (u < v) {
            cross_edges.push_back({u, v});
          }
        } else if (dist[v] == dist[u] + 1) {
          // Prevent duplicate traversal of cross-edges by checking if already a
          // parent
          if (std::find(parents[v].begin(), parents[v].end(),
                        static_cast<AtomID>(u)) == parents[v].end()) {
            parents[v].push_back(static_cast<AtomID>(u));
            cross_edges.push_back({u, v});
          }
        }
      }
    }

    for (const auto &edge : cross_edges) {
      size_t u = edge.first;
      size_t v = edge.second;

      if (dist[u] + dist[v] + 1 > static_cast<int>(max_size))
        continue;

      std::vector<std::vector<AtomID>> paths_u;
      std::vector<AtomID> current_u;
      current_u.reserve(max_size);
      get_paths(u, root, parents, current_u, paths_u);

      std::vector<std::vector<AtomID>> paths_v;
      std::vector<AtomID> current_v;
      current_v.reserve(max_size);
      get_paths(v, root, parents, current_v, paths_v);

      for (const auto &pu : paths_u) {
        for (const auto &pv : paths_v) {
          bool intersect = false;
          for (size_t i = 0; i < pu.size() - 1; ++i) {
            for (size_t j = 0; j < pv.size() - 1; ++j) {
              if (pu[i] == pv[j]) {
                intersect = true;
                break;
              }
            }
            if (intersect)
              break;
          }
          if (!intersect) {
            std::vector<AtomID> cycle;
            cycle.reserve(pu.size() + pv.size() - 1);
            cycle.push_back(static_cast<AtomID>(root));

            for (int i = static_cast<int>(pu.size()) - 2; i >= 0; --i) {
              cycle.push_back(pu[i]);
            }
            for (size_t i = 0; i < pv.size() - 1; ++i) {
              cycle.push_back(pv[i]);
            }

            if (cycle.size() >= 3 && cycle.size() <= max_size) {
              if (cycle[1] > cycle.back()) {
                std::reverse(cycle.begin() + 1, cycle.end());
              }

              if (isKingRing(graph, cycle, dist_king, q_king)) {
                std::vector<AtomID> sorted_nodes = cycle;
                std::sort(sorted_nodes.begin(), sorted_nodes.end());

                if (unique_ring_nodes.insert(sorted_nodes).second) {
                  all_cycles.push_back(cycle);
                }
              }
            }
          }
        }
      }
    }

    // Reset only the visited nodes for the next root iteration
    for (size_t v : visited) {
      dist[v] = -1;
      parents[v].clear();
    }
  }

  std::sort(all_cycles.begin(), all_cycles.end());
  all_cycles.erase(std::unique(all_cycles.begin(), all_cycles.end()),
                   all_cycles.end());

  return all_cycles;
}

} // anonymous namespace

std::map<int, size_t> MotifFinder::findRings(const NeighborGraph &graph,
                                             size_t max_size) {
  auto all_cycles = getAllShortestRings(graph, max_size);
  std::map<int, size_t> ring_counts;

  for (const auto &cycle : all_cycles) {
    ring_counts[static_cast<int>(cycle.size())]++;
  }

  return ring_counts;
}

std::vector<std::vector<AtomID>>
MotifFinder::extractCycles(const NeighborGraph &graph, size_t target_size) {
  auto all_cycles = getAllShortestRings(graph, target_size);
  std::vector<std::vector<AtomID>> exact_cycles;
  exact_cycles.reserve(all_cycles.size());

  for (auto &cycle : all_cycles) {
    if (cycle.size() == target_size) {
      exact_cycles.push_back(std::move(cycle));
    }
  }

  return exact_cycles;
}

} // namespace calculators
