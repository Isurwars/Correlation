/**
 * @file ClusterCalculator.cpp
 * @brief Implementation of the cluster analysis calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "calculators/ClusterCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <numeric>
#include <vector>

namespace correlation::calculators {

namespace {

// Static registration of the calculator in the factory
bool registered = CalculatorFactory::instance().registerCalculator(std::make_unique<ClusterCalculator>());

// Disjoint-Set (Union-Find) data structure for identifying connected components
class UnionFind {
public:
  explicit UnionFind(size_t n) : parent(n), sz(n, 1) { std::iota(parent.begin(), parent.end(), 0); }

  size_t find(size_t i) {
    size_t root = i;
    while (root != parent[root]) {
      root = parent[root];
    }
    // Path compression
    size_t curr = i;
    while (curr != root) {
      size_t nxt = parent[curr];
      parent[curr] = root;
      curr = nxt;
    }
    return root;
  }

  void unite(size_t i, size_t j) {
    size_t root_i = find(i);
    size_t root_j = find(j);
    if (root_i != root_j) {
      // Union by size
      if (sz[root_i] < sz[root_j]) {
        std::swap(root_i, root_j);
      }
      parent[root_j] = root_i;
      sz[root_i] += sz[root_j];
    }
  }

  size_t getSize(size_t i) { return sz[find(i)]; }

private:
  std::vector<size_t> parent;
  std::vector<size_t> sz;
};

} // namespace

void ClusterCalculator::calculateFrame(correlation::analysis::DistributionFunctions &df,
                                       const correlation::analysis::AnalysisSettings &settings) const {
  const auto *analyzer = df.neighbors();
  if (!analyzer)
    return;

  const auto &graph = analyzer->neighborGraph();
  size_t n_atoms = graph.nodeCount();
  if (n_atoms == 0)
    return;

  // Initialize Union-Find
  UnionFind uf(n_atoms);

  // Traverse the graph and union connected atoms
  for (size_t i = 0; i < n_atoms; ++i) {
    for (const auto &neighbor : graph.getNeighbors(i)) {
      uf.unite(i, neighbor.index);
    }
  }

  // Identify root sizes
  std::vector<size_t> root_sizes(n_atoms, 0);
  for (size_t i = 0; i < n_atoms; ++i) {
    // We only count sizes for actual roots to avoid double counting
    if (uf.find(i) == i) {
      root_sizes[i] = uf.getSize(i);
    }
  }

  // Determine the maximum cluster size for histogram bins
  size_t max_size = 0;
  for (size_t s : root_sizes) {
    if (s > max_size) {
      max_size = s;
    }
  }

  if (max_size == 0)
    return;

  // Create the histogram
  correlation::analysis::Histogram hist;
  hist.bins.resize(max_size);
  for (size_t i = 0; i < max_size; ++i) {
    hist.bins[i] = static_cast<double>(i + 1);
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

  for (size_t s : root_sizes) {
    if (s > 0) {
      // Cluster size 's' corresponds to bin index 's - 1'
      partial[s - 1] += 1.0;
    }
  }

  // Add the result safely to the distribution functions manager
  df.addHistogram("Cluster Size", std::move(hist));
}

} // namespace correlation::calculators
