/**
 * @file CNACalculator.cpp
 * @brief Implementation of the Common Neighbor Analysis (CNA) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/spatial/CNACalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <map>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
const bool REGISTERED = CalculatorFactory::registerTypeSafe<CNACalculator>("CNACalculator");

constexpr size_t K_MAX_DP_NODES = 16;
constexpr size_t K_MAX_BITMASK_NODES = 32;

/**
 * @brief Finds the longest simple path (in edges) in a small subgraph via bitmask DP O(k^2 * 2^k).
 */
[[nodiscard]] size_t longestPathBitmaskDP(size_t n_nodes,
                                          const std::array<uint32_t, K_MAX_BITMASK_NODES> &adj_mask,
                                          std::vector<uint32_t> &dp_reach) {
  const size_t num_masks = size_t{1} << n_nodes;
  dp_reach.assign(num_masks, 0U);
  for (size_t i = 0; i < n_nodes; ++i) {
    dp_reach[size_t{1} << i] = uint32_t{1} << i;
  }

  size_t best_edges = 0;
  for (size_t mask = 1; mask < num_masks; ++mask) {
    const uint32_t ends = dp_reach[mask];
    if (ends == 0U) {
      continue;
    }
    const auto edges = static_cast<size_t>(std::popcount(mask) - 1);
    best_edges = std::max(best_edges, edges);

    uint32_t rem = ends;
    while (rem != 0U) {
      const auto u_idx = static_cast<size_t>(std::countr_zero(rem));
      rem &= (rem - 1U);
      uint32_t next_candidates = adj_mask.at(u_idx) & ~static_cast<uint32_t>(mask);
      while (next_candidates != 0U) {
        const auto v_idx = static_cast<size_t>(std::countr_zero(next_candidates));
        next_candidates &= (next_candidates - 1U);
        dp_reach[mask | (size_t{1} << v_idx)] |= (uint32_t{1} << v_idx);
      }
    }
  }
  return best_edges;
}

/**
 * @brief Bounded iterative bitmask DFS fallback when a dense cutoff yields > 16 common neighbors.
 */
[[nodiscard]] size_t
longestPathBitmaskDFS(size_t n_nodes, const std::array<uint32_t, K_MAX_BITMASK_NODES> &adj_mask) {
  const size_t capped = std::min(n_nodes, K_MAX_BITMASK_NODES);
  size_t best = 0;
  size_t step_budget = 65536;

  struct StackFrame {
    uint32_t visited;
    uint32_t remaining;
    size_t depth;
  };
  std::array<StackFrame, K_MAX_BITMASK_NODES> stack{};

  for (size_t start = 0; start < capped && step_budget > 0; ++start) {
    const uint32_t start_bit = uint32_t{1} << start;
    size_t top = 0;
    stack.at(top++) = StackFrame{
        .visited = start_bit,
        .remaining = adj_mask.at(start) & ~start_bit,
        .depth = 0,
    };

    while (top > 0 && step_budget > 0) {
      auto &frame = stack.at(top - 1);
      best = std::max(best, frame.depth);
      if (frame.remaining == 0U || frame.depth + 1 >= capped) {
        --top;
        continue;
      }
      const auto nxt = static_cast<size_t>(std::countr_zero(frame.remaining));
      frame.remaining &= (frame.remaining - 1U);
      --step_budget;

      const uint32_t next_visited = frame.visited | (uint32_t{1} << nxt);
      stack.at(top++) = StackFrame{
          .visited = next_visited,
          .remaining = adj_mask.at(nxt) & ~next_visited,
          .depth = frame.depth + 1,
      };
    }
  }
  return best;
}

/**
 * @brief Builds compact bitmask adjacency among common neighbors and counts internal bonds.
 */
[[nodiscard]] size_t
buildCommonNeighborMasks(const std::vector<size_t> &common_neighbors,
                         const correlation::core::NeighborGraph &neighbor_graph,
                         std::array<uint32_t, K_MAX_BITMASK_NODES> &adj_mask) {
  const size_t n_common = common_neighbors.size();
  const size_t masked = std::min(n_common, K_MAX_BITMASK_NODES);
  adj_mask.fill(0U);
  size_t bond_count = 0;

  for (size_t idx_a_pos = 0; idx_a_pos < n_common; ++idx_a_pos) {
    const auto &neighbors_a = neighbor_graph.getNeighbors(common_neighbors[idx_a_pos]);
    for (size_t idx_b_pos = idx_a_pos + 1; idx_b_pos < n_common; ++idx_b_pos) {
      const size_t idx_b = common_neighbors[idx_b_pos];
      for (const auto &neighbor : neighbors_a) {
        if (neighbor.index == idx_b) {
          ++bond_count;
          if (idx_a_pos < masked && idx_b_pos < masked) {
            adj_mask.at(idx_a_pos) |= (uint32_t{1} << idx_b_pos);
            adj_mask.at(idx_b_pos) |= (uint32_t{1} << idx_a_pos);
          }
          break;
        }
      }
    }
  }
  return bond_count;
}

/**
 * @brief Helper to build the final CNA histogram from counted configurations.
 */
correlation::analysis::Histogram buildCNAHistogram(const std::map<std::string, real_t> &cna_counts,
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
struct ThreadLocalCNA {
  std::map<std::string, real_t> counts;
  std::vector<uint32_t> mark;
  std::vector<size_t> common;
  std::vector<uint32_t> dp_reach;
  uint32_t epoch{0};
  real_t total_pairs{0.0};
};

void processAtomCNAPairs(size_t atom_i, size_t num_atoms,
                         const correlation::core::NeighborGraph &neighbor_graph,
                         ThreadLocalCNA &local,
                         std::array<uint32_t, K_MAX_BITMASK_NODES> &adj_mask) {
  const auto &neighbors_i = neighbor_graph.getNeighbors(atom_i);
  if (++local.epoch == 0U) {
    std::ranges::fill(local.mark, 0U);
    local.epoch = 1U;
  }
  for (const auto &neighbor : neighbors_i) {
    if (neighbor.index < num_atoms) {
      local.mark[neighbor.index] = local.epoch;
    }
  }

  for (const auto &neighbor_j : neighbors_i) {
    const size_t atom_j = neighbor_j.index;
    if (atom_i >= atom_j) {
      continue;
    }

    local.common.clear();
    for (const auto &nb_j : neighbor_graph.getNeighbors(atom_j)) {
      if (nb_j.index < num_atoms && local.mark[nb_j.index] == local.epoch) {
        local.common.push_back(nb_j.index);
      }
    }

    const size_t n_bonds = buildCommonNeighborMasks(local.common, neighbor_graph, adj_mask);

    size_t n_longest = 0;
    if (n_bonds > 0) {
      n_longest = (local.common.size() <= K_MAX_DP_NODES)
                      ? longestPathBitmaskDP(local.common.size(), adj_mask, local.dp_reach)
                      : longestPathBitmaskDFS(local.common.size(), adj_mask);
    }

    std::string const index = "1" + std::to_string(local.common.size()) +
                              std::to_string(n_bonds) + std::to_string(n_longest);
    local.counts[index]++;
    local.total_pairs++;
  }
}
} // anonymous namespace

void CNACalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &dists,
    const correlation::analysis::AnalysisSettings & /*settings*/) const {
  dists.addHistogram("CNA", calculate(dists.cell(), dists.neighbors()));
}

correlation::analysis::Histogram
CNACalculator::calculate(const correlation::core::Cell &cell,
                         const correlation::analysis::StructureAnalyzer *neighbors) {
  if (neighbors == nullptr) {
    return {};
  }

  const auto &neighbor_graph = neighbors->neighborGraph();
  const size_t num_atoms = cell.atomCount();

  tbb::enumerable_thread_specific<ThreadLocalCNA> ets([num_atoms] {
    ThreadLocalCNA local;
    local.mark.assign(num_atoms, 0U);
    local.common.reserve(32);
    return local;
  });

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, num_atoms), [&](const tbb::blocked_range<size_t> &range) {
        auto &local = ets.local();
        std::array<uint32_t, K_MAX_BITMASK_NODES> adj_mask{};
        for (size_t atom_i = range.begin(); atom_i != range.end(); ++atom_i) {
          processAtomCNAPairs(atom_i, num_atoms, neighbor_graph, local, adj_mask);
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
