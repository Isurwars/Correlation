/**
 * @file SteinhardtCalculator.cpp
 * @brief Implementation of Steinhardt bond-order parameters.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/SteinhardtCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"
#include "math/SpecialFunctions.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
// NOLINTNEXTLINE(cert-err58-cpp)
const bool registered = CalculatorFactory::registerTypeSafe<SteinhardtCalculator>("SteinhardtCalculator");

struct SteinhardtParams {
  std::vector<real_t> Q4;
  std::vector<real_t> Q6;
  std::vector<real_t> W6_hat;
};

struct SingleAtomSteinhardt {
  real_t Q4;
  real_t Q6;
  real_t W6_hat;
};

struct GlobalSteinhardtFactors {
  real_t global_Q4_factor;
  real_t global_Q6_factor;
};

struct HistogramConfigs {
  size_t bins_Q;
  real_t Q_max;
  real_t dQ;
  size_t bins_W;
  real_t W_min;
  real_t W_max;
  real_t dW;
};

struct Wigner6Table {
  std::array<std::array<real_t, 13>, 13> table{};
  Wigner6Table() {
    for (int m_one = -6; m_one <= 6; ++m_one) {
      for (int m_two = -6; m_two <= 6; ++m_two) {
        int const m_three = -(m_one + m_two);
        if (m_three >= -6 && m_three <= 6) {
          table.at(m_one + 6).at(m_two + 6) = SteinhardtCalculator::wigner3j(6, 6, 6, m_one, m_two, m_three);
        } else {
          table.at(m_one + 6).at(m_two + 6) = static_cast<real_t>(0.0);
        }
      }
    }
  }
};

real_t computeW6(const std::vector<std::complex<real_t>> &q6m) {
  static const Wigner6Table wigner6;
  real_t w6_val = static_cast<real_t>(0.0);
  for (int m_one = -6; m_one <= 6; ++m_one) {
    for (int m_two = -6; m_two <= 6; ++m_two) {
      int const m_three = -(m_one + m_two);
      if (m_three < -6 || m_three > 6) {
        continue;
      }
      real_t const w3j = wigner6.table.at(m_one + 6).at(m_two + 6);
      if (w3j == static_cast<real_t>(0.0)) {
        continue;
      }

      std::complex<real_t> const prod = q6m[m_one + 6] * q6m[m_two + 6] * q6m[m_three + 6];
      w6_val += w3j * prod.real();
    }
  }
  return w6_val;
}

SingleAtomSteinhardt computeSingleAtomSteinhardt(size_t atom_idx,
                                                 const correlation::core::NeighborGraph &neighbor_graph,
                                                 GlobalSteinhardtFactors factors) {
  const auto &atom_neighbors = neighbor_graph.getNeighbors(atom_idx);
  size_t const num_neighbors = atom_neighbors.size();
  if (num_neighbors < 2) {
    return {.Q4 = static_cast<real_t>(0.0), .Q6 = static_cast<real_t>(0.0), .W6_hat = static_cast<real_t>(0.0)};
  }

  std::vector<std::complex<real_t>> q4m(9, static_cast<real_t>(0.0));
  std::vector<std::complex<real_t>> q6m(13, static_cast<real_t>(0.0));

  for (const auto &neighbor : atom_neighbors) {
    correlation::math::Vector3<real_t> r_ij = neighbor.r_ij;
    real_t const distance = neighbor.distance;
    if (distance == static_cast<real_t>(0.0)) {
      continue;
    }

    real_t const theta = std::acos(std::clamp(r_ij.z() / distance, static_cast<real_t>(-1.0), static_cast<real_t>(1.0)));
    real_t const phi = std::atan2(r_ij.y(), r_ij.x());

    for (int m_val = -4; m_val <= 4; ++m_val) {
      q4m[m_val + 4] += SteinhardtCalculator::sphericalHarmonic(4, m_val, {.theta = theta, .phi = phi});
    }
    for (int m_val = -6; m_val <= 6; ++m_val) {
      q6m[m_val + 6] += SteinhardtCalculator::sphericalHarmonic(6, m_val, {.theta = theta, .phi = phi});
    }
  }

  for (int m_val = -4; m_val <= 4; ++m_val) {
    q4m[m_val + 4] /= static_cast<real_t>(num_neighbors);
  }
  for (int m_val = -6; m_val <= 6; ++m_val) {
    q6m[m_val + 6] /= static_cast<real_t>(num_neighbors);
  }

  real_t sum_sq_4 = static_cast<real_t>(0.0);
  for (int m_val = -4; m_val <= 4; ++m_val) {
    sum_sq_4 += std::norm(q4m[m_val + 4]);
  }
  real_t const Q4_val = factors.global_Q4_factor * std::sqrt(sum_sq_4);

  real_t sum_sq_6 = static_cast<real_t>(0.0);
  for (int m_val = -6; m_val <= 6; ++m_val) {
    sum_sq_6 += std::norm(q6m[m_val + 6]);
  }
  real_t const Q6_val = factors.global_Q6_factor * std::sqrt(sum_sq_6);

  real_t const w6_val = computeW6(q6m);

  real_t W6_hat_val = static_cast<real_t>(0.0);
  if (sum_sq_6 > static_cast<real_t>(1e-12)) {
    W6_hat_val = w6_val / std::pow(sum_sq_6, static_cast<real_t>(1.5));
  }

  return {.Q4 = Q4_val, .Q6 = Q6_val, .W6_hat = W6_hat_val};
}

void normalizeHistogramMap(std::map<std::string, std::vector<real_t>> &partials, real_t factor) {
  if (factor <= static_cast<real_t>(0.0)) {
    return;
  }
  for (auto &[key, vec] : partials) {
    for (auto &val : vec) {
      val /= factor;
    }
  }
}

struct BinningConfig {
  real_t min_val;
  real_t max_val;
  real_t d_val;
};

void addValueToHistogram(std::map<std::string, std::vector<real_t>> &partials, const std::string &symbol, real_t val,
                         BinningConfig config) {
  if (val >= config.min_val && val < config.max_val) {
    auto const bin_idx = static_cast<size_t>((val - config.min_val) / config.d_val);
    auto &symbol_vec = partials[symbol];
    if (bin_idx < symbol_vec.size()) {
      symbol_vec[bin_idx] += static_cast<real_t>(1.0);
      partials["Total"][bin_idx] += static_cast<real_t>(1.0);
    }
  }
}

void initHistogramMap(std::map<std::string, std::vector<real_t>> &partials,
                      const std::vector<std::string> &element_symbols, size_t bins) {
  for (const auto &sym : element_symbols) {
    partials[sym].assign(bins, static_cast<real_t>(0.0));
  }
  partials["Total"].assign(bins, static_cast<real_t>(0.0));
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

void copyPartialsToHistogram(correlation::analysis::Histogram &hist,
                             const std::map<std::string, std::vector<real_t>> &partials) {
  for (const auto &[key, vec] : partials) {
    hist.partials[key] = vec;
  }
}

void populateHistograms(const correlation::core::Cell &cell, const correlation::analysis::StructureAnalyzer *neighbors,
                        const SteinhardtParams &params, HistogramConfigs configs,
                        correlation::analysis::Histogram &hist_Q4, correlation::analysis::Histogram &hist_Q6,
                        correlation::analysis::Histogram &hist_W6) {
  const auto &atoms = cell.atoms();
  const auto &neighbor_graph = neighbors->neighborGraph();
  size_t const num_atoms = atoms.size();

  // Collect element symbols for pre-sizing thread-local maps
  std::vector<std::string> element_symbols;
  for (const auto &elem : cell.elements()) {
    element_symbols.push_back(elem.symbol);
  }

  struct ThreadLocalHist {
    std::map<std::string, std::vector<real_t>> partials_Q4;
    std::map<std::string, std::vector<real_t>> partials_Q6;
    std::map<std::string, std::vector<real_t>> partials_W6;
    real_t num_atoms_f = static_cast<real_t>(0.0);
  };

  tbb::enumerable_thread_specific<ThreadLocalHist> ets([&]() {
    ThreadLocalHist local;
    initHistogramMap(local.partials_Q4, element_symbols, configs.bins_Q);
    initHistogramMap(local.partials_Q6, element_symbols, configs.bins_Q);
    initHistogramMap(local.partials_W6, element_symbols, configs.bins_W);
    return local;
  });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_atoms), [&](const tbb::blocked_range<size_t> &range) {
    auto &local = ets.local();

    for (size_t i = range.begin(); i != range.end(); ++i) {
      if (neighbor_graph.getNeighbors(i).size() < 2) {
        continue;
      }

      local.num_atoms_f += static_cast<real_t>(1.0);
      const std::string &symbol = atoms[i].element().symbol;

      addValueToHistogram(local.partials_Q4, symbol, params.Q4[i],
                          {.min_val = static_cast<real_t>(0.0), .max_val = configs.Q_max, .d_val = configs.dQ});
      addValueToHistogram(local.partials_Q6, symbol, params.Q6[i],
                          {.min_val = static_cast<real_t>(0.0), .max_val = configs.Q_max, .d_val = configs.dQ});
      addValueToHistogram(local.partials_W6, symbol, params.W6_hat[i],
                          {.min_val = configs.W_min, .max_val = configs.W_max, .d_val = configs.dW});
    }
  });

  // Reduce thread-local histograms
  std::map<std::string, std::vector<real_t>> partials_Q4;
  std::map<std::string, std::vector<real_t>> partials_Q6;
  std::map<std::string, std::vector<real_t>> partials_W6;

  initHistogramMap(partials_Q4, element_symbols, configs.bins_Q);
  initHistogramMap(partials_Q6, element_symbols, configs.bins_Q);
  initHistogramMap(partials_W6, element_symbols, configs.bins_W);

  real_t num_atoms_f = static_cast<real_t>(0.0);
  for (const auto &local : ets) {
    num_atoms_f += local.num_atoms_f;
    accumulateHistogramMap(partials_Q4, local.partials_Q4);
    accumulateHistogramMap(partials_Q6, local.partials_Q6);
    accumulateHistogramMap(partials_W6, local.partials_W6);
  }

  if (num_atoms_f > static_cast<real_t>(0.0)) {
    normalizeHistogramMap(partials_Q4, num_atoms_f * configs.dQ);
    normalizeHistogramMap(partials_Q6, num_atoms_f * configs.dQ);
    normalizeHistogramMap(partials_W6, num_atoms_f * configs.dW);
  }

  copyPartialsToHistogram(hist_Q4, partials_Q4);
  copyPartialsToHistogram(hist_Q6, partials_Q6);
  copyPartialsToHistogram(hist_W6, partials_W6);
}
} // namespace

std::complex<real_t> SteinhardtCalculator::sphericalHarmonic(int degree, int order, SphericalAngles angles) {
  if (order >= 0) {
    real_t const P_lm = correlation::math::sph_legendre({.degree = degree, .order = order}, angles.theta);
    return P_lm * std::polar(static_cast<real_t>(1.0), static_cast<real_t>(order) * angles.phi);
  } // For negative m: Y_l^{-m} = (-1)^m (Y_l^m)*
  int const abs_m = -order;
  real_t const P_lm = correlation::math::sph_legendre({.degree = degree, .order = abs_m}, angles.theta);
  std::complex<real_t> const Y_l_m =
      P_lm * std::polar(static_cast<real_t>(1.0), static_cast<real_t>(abs_m) * angles.phi);
  std::complex<real_t> Y_l_minus_m = std::conj(Y_l_m);
  if (abs_m % 2 != 0) {
    Y_l_minus_m = -Y_l_minus_m;
  }
  return Y_l_minus_m;
}

real_t SteinhardtCalculator::wigner3j(int j_one, int j_two, int j_three, int m_one, int m_two, int m_three) {
  if (m_one + m_two + m_three != 0) {
    return 0.0;
  }
  if (j_three < std::abs(j_one - j_two) || j_three > j_one + j_two) {
    return 0.0;
  }
  if (std::abs(m_one) > j_one || std::abs(m_two) > j_two || std::abs(m_three) > j_three) {
    return 0.0;
  }

  real_t delta = static_cast<real_t>(correlation::math::factorial(j_one + j_two - j_three) *
                                     correlation::math::factorial(j_one - j_two + j_three) *
                                     correlation::math::factorial(-j_one + j_two + j_three) /
                                     correlation::math::factorial(j_one + j_two + j_three + 1));
  delta = std::sqrt(delta);

  real_t comp = static_cast<real_t>(
      correlation::math::factorial(j_one - m_one) * correlation::math::factorial(j_one + m_one) *
      correlation::math::factorial(j_two - m_two) * correlation::math::factorial(j_two + m_two) *
      correlation::math::factorial(j_three - m_three) * correlation::math::factorial(j_three + m_three));
  comp = std::sqrt(comp);

  real_t const phase1 = ((j_one - j_two - m_three) % 2 != 0) ? -1.0 : 1.0;

  int const k_min = std::max({0, j_two - j_three - m_one, j_one - j_three + m_two});
  int const k_max = std::min({j_one + j_two - j_three, j_one - m_one, j_two + m_two});

  real_t sum = 0.0;
  for (int k = k_min; k <= k_max; ++k) {
    real_t const k_phase = (k % 2 != 0) ? -1.0 : 1.0;
    real_t const denom = static_cast<real_t>(
        correlation::math::factorial(k) * correlation::math::factorial(j_one + j_two - j_three - k) *
        correlation::math::factorial(j_one - m_one - k) * correlation::math::factorial(j_two + m_two - k) *
        correlation::math::factorial(j_three - j_two + m_one + k) *
        correlation::math::factorial(j_three - j_one - m_two + k));
    sum += k_phase / denom;
  }

  return phase1 * delta * comp * sum;
}

void SteinhardtCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                          const correlation::analysis::AnalysisSettings & /*settings*/) const {
  auto histograms = calculate(dists.cell(), dists.neighbors());
  for (auto &[name, hist] : histograms) {
    dists.addHistogram(name, std::move(hist));
  }
}

std::map<std::string, correlation::analysis::Histogram>
SteinhardtCalculator::calculate(const correlation::core::Cell &cell,
                                const correlation::analysis::StructureAnalyzer *neighbors) {
  if (neighbors == nullptr) {
    throw std::logic_error("Cannot calculate Steinhardt Parameters. Neighbor list has not been "
                           "computed.");
  }

  const auto &atoms = cell.atoms();
  const auto &neighbor_graph = neighbors->neighborGraph();
  size_t const num_atoms = atoms.size();

  // 1. Compute Q4, Q6, W6_hat for all atoms
  SteinhardtParams params;
  params.Q4.resize(num_atoms, 0.0);
  params.Q6.resize(num_atoms, 0.0);
  params.W6_hat.resize(num_atoms, 0.0);

  real_t const global_Q4_factor = static_cast<real_t>(std::sqrt(correlation::math::four_pi / 9.0));
  real_t const global_Q6_factor = static_cast<real_t>(std::sqrt(correlation::math::four_pi / 13.0));

  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_atoms), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i != range.end(); ++i) {
      auto const res = computeSingleAtomSteinhardt(
          i, neighbor_graph, {.global_Q4_factor = global_Q4_factor, .global_Q6_factor = global_Q6_factor});
      params.Q4[i] = res.Q4;
      params.Q6[i] = res.Q6;
      params.W6_hat[i] = res.W6_hat;
    }
  });

  // 2. Initialize Histograms
  size_t const bins_Q = 100;
  real_t const Q_max = 1.0;
  real_t const d_q = Q_max / bins_Q;

  size_t const bins_W = 100;
  real_t const W_min = -0.2;
  real_t const W_max = 0.2;
  real_t const d_w = (W_max - W_min) / bins_W;

  correlation::analysis::Histogram hist_Q4;
  hist_Q4.x_label = "Q4";
  hist_Q4.title = "Steinhardt Q4 Interface Parameter";
  hist_Q4.y_label = "Probability";
  hist_Q4.x_unit = "arbitrary units";
  hist_Q4.y_unit = "counts";
  hist_Q4.description = "Steinhardt Q4 Bond Orientational Order Parameter";
  hist_Q4.file_suffix = "_Q4";
  hist_Q4.bins.resize(bins_Q);
  for (size_t bin_idx = 0; bin_idx < bins_Q; ++bin_idx) {
    hist_Q4.bins[bin_idx] = static_cast<real_t>(static_cast<real_t>(bin_idx) + 0.5 * d_q);
  }

  correlation::analysis::Histogram hist_Q6;
  hist_Q6.x_label = "Q6";
  hist_Q6.title = "Steinhardt Q6 Interface Parameter";
  hist_Q6.y_label = "Probability";
  hist_Q6.x_unit = "arbitrary units";
  hist_Q6.y_unit = "counts";
  hist_Q6.description = "Steinhardt Q6 Bond Orientational Order Parameter";
  hist_Q6.file_suffix = "_Q6";
  hist_Q6.bins.resize(bins_Q);
  for (size_t bin_idx = 0; bin_idx < bins_Q; ++bin_idx) {
    hist_Q6.bins[bin_idx] = static_cast<real_t>(static_cast<real_t>(bin_idx) + 0.5 * d_q);
  }

  correlation::analysis::Histogram hist_W6;
  hist_W6.x_label = "W6_hat";
  hist_W6.title = "Steinhardt W6_hat Parameter";
  hist_W6.y_label = "Probability";
  hist_W6.x_unit = "arbitrary units";
  hist_W6.y_unit = "counts";
  hist_W6.description = "Steinhardt Normalized W6 Bond Orientational Order Parameter";
  hist_W6.file_suffix = "_W6_hat";
  hist_W6.bins.resize(bins_W);
  for (size_t bin_idx = 0; bin_idx < bins_W; ++bin_idx) {
    hist_W6.bins[bin_idx] = W_min + static_cast<real_t>(static_cast<real_t>(bin_idx) + 0.5 * d_w);
  }

  // 3. Populate Histograms
  populateHistograms(
      cell, neighbors, params,
      {.bins_Q = bins_Q, .Q_max = Q_max, .dQ = d_q, .bins_W = bins_W, .W_min = W_min, .W_max = W_max, .dW = d_w},
      hist_Q4, hist_Q6, hist_W6);

  std::map<std::string, correlation::analysis::Histogram> hists;
  hists["Q4"] = std::move(hist_Q4);
  hists["Q6"] = std::move(hist_Q6);
  hists["W6_hat"] = std::move(hist_W6);

  return hists;
}

} // namespace correlation::calculators
