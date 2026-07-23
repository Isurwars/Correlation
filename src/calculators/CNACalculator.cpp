/**
 * @file CNACalculator.cpp
 * @brief Implementation of the Common Neighbor Analysis (CNA) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/CNACalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include <algorithm>
#include <map>
#include <set>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
const bool registered = CalculatorFactory::registerTypeSafe<CNACalculator>("CNACalculator");

/**
 * @brief DFS helper to find the longest path in the common neighbor subgraph.
 *
 * Since common neighbor subgraphs are extremely small (< 12 atoms),
 * an exhaustive DFS is both correct and performant.
 *
 * @param start_node  Current node in the DFS traversal.
 * @param adj         Adjacency list for the common neighbor subgraph.
 * @return The length of the longest path (in edges) starting from @p start_node.
 */
size_t dfsLongestPath(size_t start_node, const std::map<size_t, std::vector<size_t>> &adj) {
  struct StackFrame {
    size_t node;
    size_t neighbor_idx;
    size_t max_child_len;
  };

  std::set<size_t> visited;
  visited.insert(start_node);

  std::vector<StackFrame> stack;
  stack.reserve(12);
  stack.push_back({.node = start_node, .neighbor_idx = 0, .max_child_len = 0});

  size_t final_best = 0;

  while (!stack.empty()) {
    size_t const current_node = stack.back().node;
    size_t const neighbor_idx = stack.back().neighbor_idx;
    auto const adj_iterator = adj.find(current_node);

    if (adj_iterator == adj.end() || neighbor_idx >= adj_iterator->second.size()) {
      size_t const max_child_len = stack.back().max_child_len;
      stack.pop_back();
      visited.erase(current_node);

      if (!stack.empty()) {
        stack.back().max_child_len = std::max(stack.back().max_child_len, 1 + max_child_len);
      } else {
        final_best = max_child_len;
      }
    } else {
      size_t const neighbor = adj_iterator->second[neighbor_idx];
      stack.back().neighbor_idx++;

      if (!visited.contains(neighbor)) {
        visited.insert(neighbor);
        stack.push_back({.node = neighbor, .neighbor_idx = 0, .max_child_len = 0});
      }
    }
  }

  return final_best;
}

/**
 * @brief Finds the longest continuous chain of bonds in the common neighbor
 *        subgraph by trying every node as a starting point.
 *
 * @param common_neighbors  Indices of the common neighbors.
 * @param adj               Adjacency list among common neighbors.
 * @return The longest path length (in edges).
 */
size_t findLongestChain(const std::vector<size_t> &common_neighbors, const std::map<size_t, std::vector<size_t>> &adj) {
  size_t longest = 0;
  for (size_t const start : common_neighbors) {
    size_t const len = dfsLongestPath(start, adj);
    longest = std::max(len, longest);
  }
  return longest;
}

/**
 * @brief Helper to find the common neighbors between atom_i and atom_j.
 */
std::vector<size_t> findCommonNeighbors(const std::set<size_t> &neighbors_set_i,
                                        const correlation::core::NeighborGraph &neighbor_graph,
                                        size_t atom_j) {
  const auto &neighbors_j = neighbor_graph.getNeighbors(atom_j);
  std::vector<size_t> common_neighbors;
  for (const auto &neighbor : neighbors_j) {
    if (neighbors_set_i.contains(neighbor.index)) {
      common_neighbors.push_back(neighbor.index);
    }
  }
  return common_neighbors;
}

struct CommonNeighborAdjacency {
  size_t bond_count;
  std::map<size_t, std::vector<size_t>> adjacency_list;
};

/**
 * @brief Helper to build the adjacency list and count the bonds between common neighbors.
 */
CommonNeighborAdjacency buildCommonNeighborAdjacency(
    const std::vector<size_t> &common_neighbors,
    const correlation::core::NeighborGraph &neighbor_graph) {
  size_t bond_count = 0;
  std::map<size_t, std::vector<size_t>> adjacency_list;
  size_t const n_common = common_neighbors.size();

  for (size_t idx_a_pos = 0; idx_a_pos < n_common; ++idx_a_pos) {
    size_t const idx_a = common_neighbors[idx_a_pos];
    const auto &neighbors_a = neighbor_graph.getNeighbors(idx_a);
    for (size_t idx_b_pos = idx_a_pos + 1; idx_b_pos < n_common; ++idx_b_pos) {
      size_t const idx_b = common_neighbors[idx_b_pos];
      for (const auto &neighbor : neighbors_a) {
        if (neighbor.index == idx_b) {
          bond_count++;
          adjacency_list[idx_a].push_back(idx_b);
          adjacency_list[idx_b].push_back(idx_a);
          break;
        }
      }
    }
  }
  return {.bond_count = bond_count, .adjacency_list = std::move(adjacency_list)};
}

/**
 * @brief Helper to build the final CNA histogram from counted configurations.
 */
correlation::analysis::Histogram buildCNAHistogram(
    const std::map<std::string, real_t> &cna_counts,
    real_t total_pairs) {
  correlation::analysis::Histogram hist;
  hist.x_label = "CNA Index";
  hist.title = "Common Neighbor Analysis";
  hist.y_label = "Frequency";
  hist.x_unit = "index";
  hist.y_unit = "fraction";
  hist.description = "Local structure classification via CNA.";
  hist.file_suffix = "_CNA";

  if (total_pairs > 0) {
    std::vector<std::string> keys;
    keys.reserve(cna_counts.size());
    for (const auto &[key, val] : cna_counts) {
      keys.push_back(key);
    }

    const size_t n_bins = keys.size();

    hist.bins.resize(n_bins);
    for (size_t idx = 0; idx < n_bins; ++idx) {
      hist.bins[idx] = static_cast<real_t>(idx);
    }

    for (size_t idx = 0; idx < n_bins; ++idx) {
      std::vector<real_t> values(n_bins, 0.0);
      values[idx] = cna_counts.at(keys[idx]) / total_pairs;
      hist.partials[keys[idx]] = std::move(values);
    }

    std::vector<real_t> total(n_bins, 0.0);
    for (size_t idx = 0; idx < n_bins; ++idx) {
      total[idx] = cna_counts.at(keys[idx]) / total_pairs;
    }
    hist.partials["Total"] = std::move(total);
  }

  return hist;
}
} // anonymous namespace

void CNACalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                   const correlation::analysis::AnalysisSettings & /*settings*/) const {
  dists.addHistogram("CNA", calculate(dists.cell(), dists.neighbors()));
}

correlation::analysis::Histogram CNACalculator::calculate(const correlation::core::Cell &cell,
                                                           const correlation::analysis::StructureAnalyzer *neighbors) {
  if (neighbors == nullptr) {
    return {};
  }

  const auto &neighbor_graph = neighbors->neighborGraph();
  const size_t num_atoms = cell.atomCount();

  struct ThreadLocalCNA {
    std::map<std::string, real_t> counts;
    real_t total_pairs = 0.0;
  };

  tbb::enumerable_thread_specific<ThreadLocalCNA> ets;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_atoms), [&](const tbb::blocked_range<size_t> &range) {
    auto &local = ets.local();

    for (size_t atom_i = range.begin(); atom_i != range.end(); ++atom_i) {
      const auto &neighbors_i = neighbor_graph.getNeighbors(atom_i);
      std::set<size_t> neighbors_set_i;
      for (const auto &neighbor : neighbors_i) {
        neighbors_set_i.insert(neighbor.index);
      }

      for (const auto &neighbor_j : neighbors_i) {
        size_t const atom_j = neighbor_j.index;
        if (atom_i >= atom_j) {
          continue;
        }

        std::vector<size_t> const common_neighbors = findCommonNeighbors(neighbors_set_i, neighbor_graph, atom_j);
        auto const [n_bonds, common_adj] = buildCommonNeighborAdjacency(common_neighbors, neighbor_graph);

        size_t n_longest = 0;
        if (n_bonds > 0) {
          n_longest = findLongestChain(common_neighbors, common_adj);
        }

        std::string const index = "1" + std::to_string(common_neighbors.size()) + std::to_string(n_bonds) + std::to_string(n_longest);
        local.counts[index]++;
        local.total_pairs++;
      }
    }
  });

  // Reduce thread-local results
  std::map<std::string, real_t> cna_counts;
  real_t total_pairs = 0.0;
  for (const auto &local : ets) {
    for (const auto &[key, val] : local.counts) {
      cna_counts[key] += val;
    }
    total_pairs += local.total_pairs;
  }

  return buildCNAHistogram(cna_counts, total_pairs);
}

} // namespace correlation::calculators
