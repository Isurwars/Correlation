/**
 * @file CNACalculator.cpp
 * @brief Implementation of the Common Neighbor Analysis (CNA) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "calculators/CNACalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include <map>
#include <set>
#include <vector>

namespace correlation::calculators {

namespace {
bool registered = CalculatorFactory::instance().registerCalculator(
    std::make_unique<CNACalculator>());

/**
 * @brief DFS helper to find the longest path in the common neighbor subgraph.
 *
 * Since common neighbor subgraphs are extremely small (< 12 atoms),
 * an exhaustive DFS is both correct and performant.
 *
 * @param node     Current node in the DFS traversal.
 * @param visited  Set of already-visited nodes on the current path.
 * @param adj      Adjacency list for the common neighbor subgraph.
 * @return The length of the longest path (in edges) starting from @p node.
 */
size_t dfsLongestPath(size_t node,
                      std::set<size_t> &visited,
                      const std::map<size_t, std::vector<size_t>> &adj) {
    size_t best = 0;
    auto it = adj.find(node);
    if (it == adj.end()) return 0;

    for (size_t neighbor : it->second) {
        if (visited.count(neighbor)) continue;
        visited.insert(neighbor);
        size_t len = 1 + dfsLongestPath(neighbor, visited, adj);
        if (len > best) best = len;
        visited.erase(neighbor);
    }
    return best;
}

/**
 * @brief Finds the longest continuous chain of bonds in the common neighbor
 *        subgraph by trying every node as a starting point.
 *
 * @param common_neighbors  Indices of the common neighbors.
 * @param adj               Adjacency list among common neighbors.
 * @return The longest path length (in edges).
 */
size_t findLongestChain(
    const std::vector<size_t> &common_neighbors,
    const std::map<size_t, std::vector<size_t>> &adj) {
    size_t longest = 0;
    for (size_t start : common_neighbors) {
        std::set<size_t> visited;
        visited.insert(start);
        size_t len = dfsLongestPath(start, visited, adj);
        if (len > longest) longest = len;
    }
    return longest;
}
} // anonymous namespace

void CNACalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &df,
    const correlation::analysis::AnalysisSettings &settings) const {
    df.addHistogram("CNA", calculate(df.cell(), df.neighbors()));
}

correlation::analysis::Histogram CNACalculator::calculate(
    const correlation::core::Cell& cell,
    const correlation::analysis::StructureAnalyzer* neighbors) {
    
    if (!neighbors) return {};

    const auto& neighbor_graph = neighbors->neighborGraph();
    const size_t num_atoms = cell.atomCount();

    // Map to store CNA index counts (e.g., "1421" -> count)
    std::map<std::string, double> cna_counts;
    double total_pairs = 0;

    for (size_t i = 0; i < num_atoms; ++i) {
        const auto& neighbors_i = neighbor_graph.getNeighbors(i);
        std::set<size_t> set_i;
        for (const auto& n : neighbors_i) set_i.insert(n.index);

        for (const auto& neighbor_j : neighbors_i) {
            size_t j = neighbor_j.index;
            if (i >= j) continue; // Only count each pair once

            const auto& neighbors_j = neighbor_graph.getNeighbors(j);
            std::vector<size_t> common_neighbors;
            for (const auto& n : neighbors_j) {
                if (set_i.count(n.index)) {
                    common_neighbors.push_back(n.index);
                }
            }

            size_t n_common = common_neighbors.size();
            if (n_common == 0) continue;

            // Count bonds between common neighbors and build adjacency list
            size_t n_bonds = 0;
            std::map<size_t, std::vector<size_t>> common_adj;
            for (size_t a = 0; a < n_common; ++a) {
                size_t idx_a = common_neighbors[a];
                const auto& n_a = neighbor_graph.getNeighbors(idx_a);
                for (size_t b = a + 1; b < n_common; ++b) {
                    size_t idx_b = common_neighbors[b];
                    for (const auto& nb : n_a) {
                        if (nb.index == idx_b) {
                            n_bonds++;
                            common_adj[idx_a].push_back(idx_b);
                            common_adj[idx_b].push_back(idx_a);
                            break;
                        }
                    }
                }
            }

            // Find longest continuous chain via DFS on the common neighbor subgraph
            size_t n_longest = 0;
            if (n_bonds > 0) {
                n_longest = findLongestChain(common_neighbors, common_adj);
            }

            // Standard CNA index: (1, n_common, n_bonds, n_longest)
            std::string index = "1" + std::to_string(n_common)
                                    + std::to_string(n_bonds)
                                    + std::to_string(n_longest);
            cna_counts[index]++;
            total_pairs++;
        }
    }

    // --- Build histogram with consistent bin/partial dimensions ---
    // Each unique CNA index becomes one categorical bin position.
    // All partial vectors have exactly the same length as bins.
    correlation::analysis::Histogram hist;
    hist.x_label = "CNA Index";
    hist.title = "Common Neighbor Analysis";
    hist.y_label = "Frequency";
    hist.x_unit = "index";
    hist.y_unit = "fraction";
    hist.description = "Local structure classification via CNA.";
    hist.file_suffix = "_CNA";

    if (total_pairs > 0) {
        // Collect sorted CNA index keys for deterministic ordering
        std::vector<std::string> keys;
        keys.reserve(cna_counts.size());
        for (const auto& [k, _] : cna_counts) {
            keys.push_back(k);
        }

        const size_t n_bins = keys.size();

        // Assign each CNA index a sequential categorical bin position
        hist.bins.resize(n_bins);
        for (size_t i = 0; i < n_bins; ++i) {
            hist.bins[i] = static_cast<double>(i);
        }

        // Each CNA index partial has exactly n_bins entries (zero everywhere
        // except at its own position) so that bins.size() == partial.size()
        for (size_t idx = 0; idx < n_bins; ++idx) {
            std::vector<double> values(n_bins, 0.0);
            values[idx] = cna_counts[keys[idx]] / total_pairs;
            hist.partials[keys[idx]] = std::move(values);
        }

        // "Total" partial: fraction at each bin position (sums to 1.0)
        std::vector<double> total(n_bins, 0.0);
        for (size_t idx = 0; idx < n_bins; ++idx) {
            total[idx] = cna_counts[keys[idx]] / total_pairs;
        }
        hist.partials["Total"] = std::move(total);
    }

    return hist;
}

} // namespace correlation::calculators
