/**
 * @file ClusterCalculator.cpp
 * @brief Implementation of the cluster analysis calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/ClusterCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <algorithm>
#include <numeric>
#include <vector>

namespace correlation::calculators {

namespace {

// Static registration of the calculator in the factory
const bool registered = CalculatorFactory::registerTypeSafe<ClusterCalculator>("ClusterCalculator");

// Disjoint-Set (Union-Find) data structure for identifying connected components
class UnionFind {
public:
  explicit UnionFind(size_t num_nodes) : parent(num_nodes), sz(num_nodes, 1) { std::ranges::iota(parent, 0); }

  size_t find(size_t node) {
    size_t root = node;
    while (root != parent[root]) {
      root = parent[root];
    }
    // Path compression
    size_t curr = node;
    while (curr != root) {
      size_t const nxt = parent[curr];
      parent[curr] = root;
      curr = nxt;
    }
    return root;
  }

  void unite(size_t node_i, size_t node_j) {
    size_t root_i = find(node_i);
    size_t root_j = find(node_j);
    if (root_i != root_j) {
      // Union by size
      if (sz[root_i] < sz[root_j]) {
        std::swap(root_i, root_j);
      }
      parent[root_j] = root_i;
      sz[root_i] += sz[root_j];
    }
  }

  size_t getSize(size_t node) { return sz[find(node)]; }

private:
  std::vector<size_t> parent;
  std::vector<size_t> sz;
};

} // namespace

void ClusterCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                       const correlation::analysis::AnalysisSettings & /*settings*/) const {
  const auto *analyzer = dists.neighbors();
  if (analyzer == nullptr) {
    return;
  }

  const auto &graph = analyzer->neighborGraph();
  size_t const n_atoms = graph.nodeCount();
  if (n_atoms == 0) {
    return;
  }

  // Initialize Union-Find
  UnionFind union_find(n_atoms);

  // Traverse the graph and union connected atoms
  for (size_t i = 0; i < n_atoms; ++i) {
    for (const auto &neighbor : graph.getNeighbors(i)) {
      union_find.unite(i, neighbor.index);
    }
  }

  // Identify root sizes
  std::vector<size_t> root_sizes(n_atoms, 0);
  for (size_t i = 0; i < n_atoms; ++i) {
    // We only count sizes for actual roots to avoid real_t counting
    if (union_find.find(i) == i) {
      root_sizes[i] = union_find.getSize(i);
    }
  }

  // Determine the maximum cluster size for histogram bins
  size_t max_size = 0;
  for (size_t const size : root_sizes) {
    max_size = std::max(size, max_size);
  }

  if (max_size == 0) {
    return;
  }

  // Create the histogram
  correlation::analysis::Histogram hist;
  hist.bins.resize(max_size);
  for (size_t i = 0; i < max_size; ++i) {
    hist.bins[i] = static_cast<real_t>(i + 1);
  }
  hist.title = "Cluster Size Distribution";
  hist.x_label = "Cluster Size";
  hist.y_label = "Count";
  hist.x_unit = "atoms";
  hist.y_unit = "clusters";
  hist.description = "Distribution of connected cluster sizes.";
  hist.file_suffix = "_clusters";

  auto &partial = hist.partials["Total"];
  partial.assign(max_size, 0.0);

  for (size_t const size : root_sizes) {
    if (size > 0) {
      // Cluster size 's' corresponds to bin index 's - 1'
      partial[size - 1] += 1.0;
    }
  }

  // Add the result safely to the distribution functions manager
  dists.addHistogram("Cluster Size", std::move(hist));
}

} // namespace correlation::calculators
