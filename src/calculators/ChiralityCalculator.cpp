/**
 * @file ChiralityCalculator.cpp
 * @brief Implementation of the Chiral Order Parameter (COP) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/ChiralityCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
const bool registered = CalculatorFactory::registerTypeSafe<ChiralityCalculator>("ChiralityCalculator");

struct BinningConfig {
  real_t min_val;
  real_t max_val;
  real_t d_val;
};

void initHistogramMap(std::map<std::string, std::vector<real_t>> &partials,
                      const std::vector<std::string> &element_symbols, size_t bins) {
  for (const auto &sym : element_symbols) {
    partials[sym].assign(bins, 0.0);
  }
  partials["Total"].assign(bins, 0.0);
}

void addValueToHistogram(std::map<std::string, std::vector<real_t>> &partials, const std::string &symbol, real_t val,
                         BinningConfig config, size_t bins) {
  if (val >= config.min_val && val <= config.max_val) {
    auto bin_idx = static_cast<size_t>((val - config.min_val) / config.d_val);
    if (bin_idx >= bins) {
      bin_idx = bins - 1;
    }
    partials[symbol][bin_idx] += 1.0;
    partials["Total"][bin_idx] += 1.0;
  }
}

void accumulateHistogramMap(std::map<std::string, std::vector<real_t>> &dest,
                            const std::map<std::string, std::vector<real_t>> &src) {
  for (const auto &[key, vec] : src) {
    auto &dest_vec = dest[key];
    for (size_t bin_idx = 0; bin_idx < vec.size(); ++bin_idx) {
      dest_vec[bin_idx] += vec[bin_idx];
    }
  }
}

void normalizeHistogramMap(std::map<std::string, std::vector<real_t>> &partials, real_t factor) {
  if (factor <= 0.0) {
    return;
  }
  for (auto &[key, vec] : partials) {
    for (auto &val : vec) {
      val /= factor;
    }
  }
}

void copyPartialsToHistogram(correlation::analysis::Histogram &hist,
                             const std::map<std::string, std::vector<real_t>> &partials) {
  for (const auto &[key, vec] : partials) {
    hist.partials[key] = vec;
  }
}
} // namespace

void ChiralityCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                         const correlation::analysis::AnalysisSettings & /*settings*/) const {
  dists.addHistogram("COP", calculate(dists.cell(), dists.neighbors()));
}

real_t ChiralityCalculator::computeSingleAtomChirality(size_t atom_idx, const correlation::core::Cell & /*cell*/,
                                                       const correlation::analysis::StructureAnalyzer *neighbors) {
  if (neighbors == nullptr) {
    return 0.0;
  }

  const auto &neighbor_graph = neighbors->neighborGraph();
  const auto &atom_neighbors = neighbor_graph.getNeighbors(atom_idx);

  // We need at least 3 valid neighbors (distance > 0.0)
  std::vector<correlation::core::Neighbor> valid_neighbors;
  valid_neighbors.reserve(atom_neighbors.size());
  for (const auto &neighbor : atom_neighbors) {
    if (neighbor.distance > 0.0) {
      valid_neighbors.push_back(neighbor);
    }
  }

  if (valid_neighbors.size() < 3) {
    return 0.0;
  }

  // Sort neighbors by distance in ascending order to find the 3 nearest neighbors
  std::ranges::sort(valid_neighbors,
                    [](const correlation::core::Neighbor &neighbor_a, const correlation::core::Neighbor &neighbor_b) {
                      return neighbor_a.distance < neighbor_b.distance;
                    });

  const auto &r_1 = valid_neighbors[0].r_ij;
  const auto &r_2 = valid_neighbors[1].r_ij;
  const auto &r_3 = valid_neighbors[2].r_ij;

  real_t const d_1 = valid_neighbors[0].distance;
  real_t const d_2 = valid_neighbors[1].distance;
  real_t const d_3 = valid_neighbors[2].distance;

  if (d_1 <= 0.0 || d_2 <= 0.0 || d_3 <= 0.0) {
    return 0.0;
  }

  // Calculate scalar triple product: (r1 x r2) . r3
  real_t const triple_product = correlation::math::dot(correlation::math::cross(r_1, r_2), r_3);

  // Normalize by lengths of the vectors
  return triple_product / (d_1 * d_2 * d_3);
}

correlation::analysis::Histogram
ChiralityCalculator::calculate(const correlation::core::Cell &cell,
                               const correlation::analysis::StructureAnalyzer *neighbors) {
  if (neighbors == nullptr) {
    throw std::logic_error("Cannot calculate Chiral Order Parameter. Neighbor list has not been computed.");
  }

  const auto &atoms = cell.atoms();
  size_t const num_atoms = atoms.size();

  // Compute local chirality for all atoms
  std::vector<real_t> chiralities(num_atoms, 0.0);
  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_atoms), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i != range.end(); ++i) {
      chiralities[i] = computeSingleAtomChirality(i, cell, neighbors);
    }
  });

  // Setup histogram configuration
  size_t const bins = 100;
  real_t const min_val = -1.0;
  real_t const max_val = 1.0;
  real_t const d_val = (max_val - min_val) / static_cast<real_t>(bins);

  correlation::analysis::Histogram hist;
  hist.title = "Chiral Order Parameter Distribution";
  hist.x_label = "Chiral Order Parameter (\\chi)";
  hist.y_label = "Probability";
  hist.x_unit = "";
  hist.y_unit = "probability";
  hist.description = "Chiral Order Parameter (COP) distribution";
  hist.file_suffix = "_cop";
  hist.bins.resize(bins);
  for (size_t i = 0; i < bins; ++i) {
    hist.bins[i] = min_val + (static_cast<real_t>(i) + static_cast<real_t>(0.5)) * d_val;
  }

  // Pre-size thread-local maps for element symbols
  std::vector<std::string> element_symbols;
  for (const auto &elem : cell.elements()) {
    element_symbols.push_back(elem.symbol);
  }

  struct ThreadLocalHist {
    std::map<std::string, std::vector<real_t>> partials;
    real_t num_atoms_f = 0.0;
  };

  tbb::enumerable_thread_specific<ThreadLocalHist> ets([&]() {
    ThreadLocalHist local;
    initHistogramMap(local.partials, element_symbols, bins);
    return local;
  });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_atoms), [&](const tbb::blocked_range<size_t> &range) {
    auto &local = ets.local();
    for (size_t i = range.begin(); i != range.end(); ++i) {
      local.num_atoms_f += 1.0;
      const std::string &symbol = atoms[i].element().symbol;
      addValueToHistogram(local.partials, symbol, chiralities[i],
                          {.min_val = min_val, .max_val = max_val, .d_val = d_val}, bins);
    }
  });

  // Reduce thread-local histograms
  std::map<std::string, std::vector<real_t>> partials;
  initHistogramMap(partials, element_symbols, bins);
  real_t total_atoms_f = 0.0;

  for (const auto &local_hist : ets) {
    accumulateHistogramMap(partials, local_hist.partials);
    total_atoms_f += local_hist.num_atoms_f;
  }

  // Normalize histograms to get probability distribution
  normalizeHistogramMap(partials, total_atoms_f);
  copyPartialsToHistogram(hist, partials);

  return hist;
}

} // namespace correlation::calculators
