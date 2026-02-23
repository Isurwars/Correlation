// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/CNCalculator.hpp"
#include <numeric>
#include <stdexcept>

Histogram CNCalculator::calculate(const Cell &cell,
                                  const StructureAnalyzer *neighbors) {
  if (!neighbors) {
    throw std::logic_error(
        "Cannot calculate Coordination Number. Neighbor list has not been "
        "computed.");
  }

  const auto &atoms = cell.atoms();
  const auto &all_neighbors = neighbors->neighbors();

  std::map<std::string, std::vector<int>> partial_dists;
  size_t max_cn = 0;

  for (size_t i = 0; i < atoms.size(); ++i) {
    const auto &central_atom = atoms[i];
    const std::string &central_symbol = central_atom.element().symbol;

    std::map<std::string, int> neighbor_counts_for_this_atom;
    for (const auto &neighbor : all_neighbors[i]) {
      const auto &neighbor_atom = atoms[neighbor.index];
      neighbor_counts_for_this_atom[neighbor_atom.element().symbol]++;
    }

    for (const auto &[neighbor_symbol, count] : neighbor_counts_for_this_atom) {
      std::string key = central_symbol + "-" + neighbor_symbol;
      if (static_cast<size_t>(count) >= partial_dists[key].size()) {
        partial_dists[key].resize(count + 1, 0);
      }
      partial_dists[key][count]++;
      if (static_cast<size_t>(count) > max_cn) {
        max_cn = count;
      }
    }
  }

  Histogram cn_histogram;
  const size_t num_bins = max_cn + 3;
  cn_histogram.bin_label = "# neighbors";

  cn_histogram.bins.resize(num_bins);
  std::iota(cn_histogram.bins.begin(), cn_histogram.bins.end(), 0.0);

  for (auto &[key, dist_vector] : partial_dists) {
    dist_vector.resize(num_bins, 0);
    cn_histogram.partials[key].assign(dist_vector.begin(), dist_vector.end());
  }

  const auto &elements = cell.elements();
  std::vector<double> any_any_dist(num_bins, 0.0);

  for (const auto &elem : elements) {
    std::string element_any_key = elem.symbol + "-Any";
    std::vector<double> element_any_dist(num_bins, 0.0);
    bool found_any = false;

    for (const auto &[key, dist_vector] : cn_histogram.partials) {
      if (key.rfind(elem.symbol + "-", 0) == 0) {
        found_any = true;
        for (size_t b = 0; b < num_bins; ++b) {
          element_any_dist[b] += dist_vector[b];
        }
      }
    }

    if (found_any) {
      cn_histogram.partials[element_any_key] = element_any_dist;
      for (size_t b = 0; b < num_bins; ++b) {
        any_any_dist[b] += element_any_dist[b];
      }
    }
  }

  cn_histogram.partials["Any-Any"] = any_any_dist;

  return cn_histogram;
}
