/**
 * @file TenWoldeFrenkelCalculator.cpp
 * @brief Implementation of Ten Wolde-Frenkel bond coherence and crystalline nucleus identification.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/order/TenWoldeFrenkelCalculator.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace correlation::calculators {

namespace {

const bool REGISTERED =
    CalculatorFactory::registerTypeSafe<TenWoldeFrenkelCalculator>("TenWoldeFrenkelCalculator");

[[nodiscard]] inline real_t computePairCoherence(const std::complex<real_t> *q6m_i,
                                                 const std::complex<real_t> *q6m_j) noexcept {
  real_t chi = static_cast<real_t>(0.0);
  for (size_t m = 0; m < 13; ++m) {
    chi += q6m_i[m].real() * q6m_j[m].real() + q6m_i[m].imag() * q6m_j[m].imag();
  }
  return chi;
}

} // namespace

void TenWoldeFrenkelCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &dists,
    const correlation::analysis::AnalysisSettings &settings) const {
  if (dists.neighbors() == nullptr) {
    return;
  }
  auto hist = calculate(dists.cell(), dists.neighbors(),
                        {
                            .chi_cut = settings.twf_chi_cut,
                            .nconn_min = settings.twf_nconn_min,
                        });
  dists.addHistogram("TWF_Nconn", std::move(hist));
}

std::vector<real_t> TenWoldeFrenkelCalculator::calculatePerAtom(
    const correlation::core::Cell &cell, const correlation::analysis::StructureAnalyzer *neighbors,
    TenWoldeFrenkelParams params, const SteinhardtCalculator::SteinhardtParams *steinhardt_params) {
  if (neighbors == nullptr) {
    throw std::logic_error("Cannot calculate TWF parameters. Neighbor list has not been computed.");
  }

  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();
  if (num_atoms == 0) {
    return {};
  }

  const auto &neighbor_graph = neighbors->neighborGraph();

  std::optional<SteinhardtCalculator::SteinhardtParams> local_steinhardt;
  const std::vector<std::complex<real_t>> *q6m_norm_ptr = nullptr;

  if (steinhardt_params != nullptr && steinhardt_params->q6m_normalized.has_value() &&
      steinhardt_params->q6m_normalized->size() == num_atoms * 13) {
    q6m_norm_ptr = &(*steinhardt_params->q6m_normalized);
  } else {
    local_steinhardt = SteinhardtCalculator::calculatePerAtom(cell, neighbors, {.keep_q6m = true});
    if (local_steinhardt->q6m_normalized.has_value()) {
      q6m_norm_ptr = &(*local_steinhardt->q6m_normalized);
    }
  }

  if (q6m_norm_ptr == nullptr) {
    return std::vector<real_t>(num_atoms, static_cast<real_t>(0.0));
  }

  const auto &q6m_norm = *q6m_norm_ptr;
  std::vector<real_t> n_conn(num_atoms, static_cast<real_t>(0.0));

  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_atoms),
                    [&](const tbb::blocked_range<size_t> &range) {
                      for (size_t i = range.begin(); i != range.end(); ++i) {
                        const auto &atom_neighbors = neighbor_graph.getNeighbors(i);
                        if (atom_neighbors.empty()) {
                          continue;
                        }
                        const auto *q6m_i = &q6m_norm[i * 13];
                        size_t coherent_count = 0;
                        for (const auto &neighbor : atom_neighbors) {
                          const size_t j = neighbor.index;
                          if (j >= num_atoms) {
                            continue;
                          }
                          const auto *q6m_j = &q6m_norm[j * 13];
                          const real_t chi_ij = computePairCoherence(q6m_i, q6m_j);
                          if (chi_ij >= params.chi_cut) {
                            ++coherent_count;
                          }
                        }
                        n_conn[i] = static_cast<real_t>(coherent_count);
                      }
                    });

  return n_conn;
}

correlation::analysis::Histogram TenWoldeFrenkelCalculator::calculate(
    const correlation::core::Cell &cell, const correlation::analysis::StructureAnalyzer *neighbors,
    TenWoldeFrenkelParams params, const SteinhardtCalculator::SteinhardtParams *steinhardt_params) {
  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();
  if (num_atoms == 0) {
    return {};
  }

  const std::vector<real_t> n_conn = calculatePerAtom(cell, neighbors, params, steinhardt_params);
  if (n_conn.empty()) {
    return {};
  }

  constexpr size_t K_BINS = 21; // 0 through 20
  constexpr real_t K_MIN_VAL = static_cast<real_t>(0.0);
  constexpr real_t K_MAX_VAL = static_cast<real_t>(21.0);

  size_t solid_like_count = 0;
  for (const auto val : n_conn) {
    if (val >= static_cast<real_t>(params.nconn_min)) {
      ++solid_like_count;
    }
  }
  const real_t solid_like_fraction =
      static_cast<real_t>(solid_like_count) / static_cast<real_t>(num_atoms);

  std::ostringstream desc_stream;
  desc_stream << "Ten Wolde-Frenkel bond coherence (solid-like fraction: " << std::fixed
              << std::setprecision(4) << solid_like_fraction << ")";

  correlation::analysis::Histogram hist;
  hist.title = "Ten Wolde-Frenkel Bond Coherence";
  hist.x_label = "Coherent Neighbors N_conn";
  hist.y_label = "Probability";
  hist.x_unit = "count";
  hist.y_unit = "probability";
  hist.description = desc_stream.str();
  hist.file_suffix = "_twf_nconn";
  hist.bins.resize(K_BINS);
  for (size_t i = 0; i < K_BINS; ++i) {
    hist.bins[i] = static_cast<real_t>(i);
  }

  std::vector<std::string> element_symbols;
  for (const auto &elem : cell.elements()) {
    element_symbols.push_back(elem.symbol);
  }

  std::map<std::string, std::vector<real_t>> partials;
  for (const auto &sym : element_symbols) {
    partials[sym].assign(K_BINS, static_cast<real_t>(0.0));
  }
  partials["Total"].assign(K_BINS, static_cast<real_t>(0.0));

  for (size_t i = 0; i < num_atoms; ++i) {
    const auto bin_idx = static_cast<size_t>(
        std::clamp(std::floor(n_conn[i]), K_MIN_VAL, K_MAX_VAL - static_cast<real_t>(1.0)));
    const std::string &sym = atoms[i].element().symbol;
    partials[sym][bin_idx] += static_cast<real_t>(1.0);
    partials["Total"][bin_idx] += static_cast<real_t>(1.0);
  }

  const real_t inv_total = static_cast<real_t>(1.0) / static_cast<real_t>(num_atoms);
  for (auto &[sym, counts] : partials) {
    for (auto &count : counts) {
      count *= inv_total;
    }
    hist.partials[sym] = counts;
  }

  return hist;
}

} // namespace correlation::calculators
