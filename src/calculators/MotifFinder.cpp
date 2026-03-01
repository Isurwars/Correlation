// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/MotifFinder.hpp"
#include <algorithm>
#include <queue>

namespace calculators {

namespace {

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

std::vector<std::vector<AtomID>> getAllShortestRings(const NeighborGraph &graph,
                                                     size_t max_size) {
  std::vector<std::vector<AtomID>> all_cycles;
  if (max_size < 3)
    return all_cycles;

  for (size_t root = 0; root < graph.nodeCount(); ++root) {
    std::vector<int> dist(graph.nodeCount(), -1);
    std::vector<std::vector<AtomID>> parents(graph.nodeCount());
    std::queue<size_t> q;

    dist[root] = 0;
    q.push(root);

    std::vector<std::pair<size_t, size_t>> cross_edges;

    while (!q.empty()) {
      size_t u = q.front();
      q.pop();

      for (const auto &neighbor : graph.getNeighbors(u)) {
        size_t v = neighbor.index;
        if (v < root)
          continue;

        if (dist[v] == -1) {
          dist[v] = dist[u] + 1;
          parents[v].push_back(static_cast<AtomID>(u));
          if (2 * dist[v] + 1 <= static_cast<int>(max_size)) {
            q.push(v);
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
      get_paths(u, root, parents, current_u, paths_u);

      std::vector<std::vector<AtomID>> paths_v;
      std::vector<AtomID> current_v;
      get_paths(v, root, parents, current_v, paths_v);

      for (const auto &pu : paths_u) {
        for (const auto &pv : paths_v) {
          bool intersect = false;
          for (size_t i = 0; i < pu.size() - 1; ++i) {
            if (std::find(pv.begin(), pv.end() - 1, pu[i]) != pv.end() - 1) {
              intersect = true;
              break;
            }
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
              all_cycles.push_back(cycle);
            }
          }
        }
      }
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
