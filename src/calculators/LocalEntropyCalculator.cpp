/**
 * @file LocalEntropyCalculator.cpp
 * @brief Implementation of the Local Entropic Fingerprint order parameter.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/LocalEntropyCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include "math/Precision.hpp"

#include <cmath>
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
// NOLINTNEXTLINE(cert-err58-cpp)
const bool registered = CalculatorFactory::registerTypeSafe<LocalEntropyCalculator>("LocalEntropyCalculator");

struct LocalEntropyParams {
  real_t cutoff;
  real_t sigma;
};

real_t computeSingleAtomEntropy(size_t atom_idx, const correlation::core::Cell &cell,
                                const correlation::analysis::StructureAnalyzer *neighbors,
                                const LocalEntropyParams &params) {
  real_t const cutoff = params.cutoff;
  real_t const sigma = params.sigma;
  const auto &neighbor_graph = neighbors->neighborGraph();
  const auto &atom_neighbors = neighbor_graph.getNeighbors(atom_idx);

  real_t const vol = cell.volume();
  if (vol <= 0.0) {
    return 0.0;
  }
  real_t const density = static_cast<real_t>(cell.atomCount()) / vol;

  // We integrate from r = dr to Rc with step size dr
  real_t const dr_val = 0.02;
  auto const num_steps = static_cast<size_t>(std::ceil(cutoff / dr_val));

  real_t integral = 0.0;

  // Gaussian prefactor: 1 / sqrt(2 * pi * sigma^2)
  real_t const gaussian_prefactor =
      static_cast<real_t>(1.0) / (sigma * std::sqrt(static_cast<real_t>(2.0) * correlation::math::pi));

  for (size_t step = 1; step <= num_steps; ++step) {
    real_t const r_val = static_cast<real_t>(step) * dr_val;
    if (r_val > cutoff) {
      break;
    }

    real_t g_val = 0.0;
    for (const auto &neighbor : atom_neighbors) {
      real_t const r_ij = neighbor.distance;
      if (r_ij == 0.0) {
        continue;
      }
      // Gaussian kernel with periodic boundary corrections at r=0
      real_t const diff = r_val - r_ij;
      real_t const sum_r = r_val + r_ij;
      real_t const term1 = static_cast<real_t>(std::exp(-(diff * diff) / (2.0 * sigma * sigma)));
      real_t const term2 = static_cast<real_t>(std::exp(-(sum_r * sum_r) / (2.0 * sigma * sigma)));
      g_val += term1 + term2;
    }

    g_val *= gaussian_prefactor;

    // Normalization: 4 * pi * density * r^2
    real_t const norm = static_cast<real_t>(4.0) * correlation::math::pi * density * r_val * r_val;
    if (norm > 0.0) {
      g_val /= norm;
    } else {
      g_val = 0.0;
    }

    real_t integrand = 0.0;
    if (g_val > 1e-10) {
      integrand = (g_val * std::log(g_val) - g_val + static_cast<real_t>(1.0)) * r_val * r_val;
    } else {
      integrand = r_val * r_val;
    }

    // Trapezoidal rule integration
    real_t weight = dr_val;
    if (step == num_steps) {
      weight = static_cast<real_t>(0.5) * dr_val;
    }
    integral += integrand * weight;
  }

  return static_cast<real_t>(-2.0) * correlation::math::pi * density * integral;
}

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
                         BinningConfig config) {
  if (val >= config.min_val && val < config.max_val) {
    auto const bin_idx = static_cast<size_t>((val - config.min_val) / config.d_val);
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

void LocalEntropyCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                            const correlation::analysis::AnalysisSettings &settings) const {
  auto hist = calculate(dists.cell(), dists.neighbors(), settings.lef_cutoff, settings.lef_sigma);
  dists.addHistogram("LEF", std::move(hist));
}

correlation::analysis::Histogram
LocalEntropyCalculator::calculate(const correlation::core::Cell &cell,
                                  const correlation::analysis::StructureAnalyzer *neighbors, real_t cutoff,
                                  real_t sigma) {
  if (neighbors == nullptr) {
    throw std::logic_error("Cannot calculate Local Entropy. Neighbor list has not been computed.");
  }

  const auto &atoms = cell.atoms();
  size_t const num_atoms = atoms.size();

  // Compute local entropy for all atoms
  std::vector<real_t> entropies(num_atoms, 0.0);
  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_atoms), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i != range.end(); ++i) {
      entropies[i] = computeSingleAtomEntropy(i, cell, neighbors, {.cutoff = cutoff, .sigma = sigma});
    }
  });

  // Setup histogram configuration
  size_t const bins = 150;
  real_t const min_val = -15.0;
  real_t const max_val = 0.0;
  real_t const d_val = (max_val - min_val) / static_cast<real_t>(bins);

  correlation::analysis::Histogram hist;
  hist.title = "Local Entropic Fingerprint";
  hist.x_label = "Local Entropy (s_i / k_B)";
  hist.y_label = "Probability";
  hist.x_unit = "k_B";
  hist.y_unit = "probability";
  hist.description = "Local Entropic Fingerprint order parameter distribution";
  hist.file_suffix = "_lef";
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
      addValueToHistogram(local.partials, symbol, entropies[i],
                          {.min_val = min_val, .max_val = max_val, .d_val = d_val});
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
