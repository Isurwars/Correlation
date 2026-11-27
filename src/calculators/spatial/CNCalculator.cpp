/**
 * @file CNCalculator.cpp
 * @brief Implementation of the coordination number calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/spatial/CNCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <stdexcept>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
const bool REGISTERED = CalculatorFactory::registerTypeSafe<CNCalculator>("CNCalculator");

/**
 * @brief Helper to increment a coordination-number histogram bin, growing the vector as needed.
 */
void incrementCoordinationBin(std::map<std::string, std::vector<int>> &dists,
                              const std::string &key, size_t count) {
  auto &vec = dists[key];
  if (count >= vec.size()) {
    vec.resize(count + 1, 0);
  }
  vec[count]++;
}
} // namespace

void CNCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &dists,
    const correlation::analysis::AnalysisSettings & /*settings*/) const {
  dists.addHistogram("CN", calculate(dists.cell(), dists.neighbors()));
}

correlation::analysis::Histogram
CNCalculator::calculate(const correlation::core::Cell &cell,
                        const correlation::analysis::StructureAnalyzer *neighbors) {
  if (neighbors == nullptr) {
    throw std::logic_error("Cannot calculate Coordination Number. Neighbor list has not been "
                           "computed.");
  }

  const auto &atoms = cell.atoms();
  const auto &elements = cell.elements();
  const auto &neighbor_graph = neighbors->neighborGraph();

  std::map<std::string, std::vector<int>> partial_dists;
  size_t max_cn = 0;

  for (size_t i = 0; i < atoms.size(); ++i) {
    const std::string &central_symbol = atoms[i].element().symbol;
    const auto neighbor_indices = neighbor_graph.getNeighborIndices(i);

    std::map<std::string, size_t> neighbor_counts_for_this_atom;
    for (const auto neighbor_idx : neighbor_indices) {
      neighbor_counts_for_this_atom[atoms[neighbor_idx].element().symbol]++;
    }

    for (const auto &elem : elements) {
      const auto iter = neighbor_counts_for_this_atom.find(elem.symbol);
      const size_t count = (iter != neighbor_counts_for_this_atom.end()) ? iter->second : 0;
      std::string key = central_symbol;
      key += '-';
      key += elem.symbol;
      incrementCoordinationBin(partial_dists, key, count);
    }

    const size_t total_cn = neighbor_indices.size();
    incrementCoordinationBin(partial_dists, central_symbol + "-Any", total_cn);
    incrementCoordinationBin(partial_dists, "Any-Any", total_cn);
    max_cn = std::max(max_cn, total_cn);
  }

  correlation::analysis::Histogram cn_histogram;
  const size_t num_bins = max_cn + 3;
  cn_histogram.x_label = "CN";
  cn_histogram.title = "Coordination Number Distribution";
  cn_histogram.y_label = "Count";
  cn_histogram.x_unit = "neighbors";
  cn_histogram.y_unit = "Count";
  cn_histogram.description = "Coordination Number";
  cn_histogram.file_suffix = "_CN";

  cn_histogram.bins.resize(num_bins);
  for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
    cn_histogram.bins[bin_idx] = static_cast<real_t>(bin_idx);
  }

  if (!partial_dists.contains("Any-Any")) {
    partial_dists["Any-Any"].assign(num_bins, 0);
  }

  for (auto &[key, dist_vector] : partial_dists) {
    dist_vector.resize(num_bins, 0);
    auto &target_vec = cn_histogram.partials[key];
    target_vec.clear();
    target_vec.reserve(num_bins);
    for (int val : dist_vector) {
      target_vec.push_back(static_cast<real_t>(val));
    }
  }

  cn_histogram.partials["Total"] = cn_histogram.partials["Any-Any"];

  return cn_histogram;
}

} // namespace correlation::calculators
