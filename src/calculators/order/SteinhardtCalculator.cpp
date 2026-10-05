/**
 * @file SteinhardtCalculator.cpp
 * @brief Implementation of Steinhardt bond-order parameters.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/order/SteinhardtCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"
#include "math/SpecialFunctions.hpp"

#if defined(CORRELATION_USE_CUDA) || defined(CORRELATION_USE_HIP)
#include "calculators/gpu/GPUSteinhardtCalculator.hpp"
#endif
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
const bool REGISTERED =
    CalculatorFactory::registerTypeSafe<SteinhardtCalculator>("SteinhardtCalculator");

using SteinhardtParams = SteinhardtCalculator::SteinhardtParams;
using SingleAtomSteinhardt = SteinhardtCalculator::SingleAtomSteinhardt;

struct GlobalSteinhardtFactors {
  real_t global_q4_factor;
  real_t global_q6_factor;
};

struct HistogramConfigs {
  size_t bins_q{100};
  real_t q_max{1.0};
  real_t step_q{0.01};
  size_t bins_w4{100};
  real_t w4_min{-0.5};
  real_t w4_max{0.5};
  real_t step_w4{0.01};
  size_t bins_w6{100};
  real_t w6_min{-0.2};
  real_t w6_max{0.2};
  real_t step_w6{0.004};
};

struct Wigner4Table {
  std::array<std::array<real_t, 9>, 9> table{};
  Wigner4Table() {
    for (int m_one = -4; m_one <= 4; ++m_one) {
      for (int m_two = -4; m_two <= 4; ++m_two) {
        int const m_three = -(m_one + m_two);
        if (m_three >= -4 && m_three <= 4) {
          table.at(m_one + 4).at(m_two + 4) = SteinhardtCalculator::wigner3j({
              .j_one = 4,
              .j_two = 4,
              .j_three = 4,
              .m_one = m_one,
              .m_two = m_two,
              .m_three = m_three,
          });
        } else {
          table.at(m_one + 4).at(m_two + 4) = static_cast<real_t>(0.0);
        }
      }
    }
  }
};

struct Wigner6Table {
  std::array<std::array<real_t, 13>, 13> table{};
  Wigner6Table() {
    for (int m_one = -6; m_one <= 6; ++m_one) {
      for (int m_two = -6; m_two <= 6; ++m_two) {
        int const m_three = -(m_one + m_two);
        if (m_three >= -6 && m_three <= 6) {
          table.at(m_one + 6).at(m_two + 6) = SteinhardtCalculator::wigner3j({
              .j_one = 6,
              .j_two = 6,
              .j_three = 6,
              .m_one = m_one,
              .m_two = m_two,
              .m_three = m_three,
          });
        } else {
          table.at(m_one + 6).at(m_two + 6) = static_cast<real_t>(0.0);
        }
      }
    }
  }
};

real_t computeW4(const std::array<std::complex<real_t>, 9> &q4m) {
  static const Wigner4Table WIGNER4;
  auto w4_val = static_cast<real_t>(0.0);
  for (size_t idx_one = 0; idx_one < 9; ++idx_one) {
    for (size_t idx_two = 0; idx_two < 9; ++idx_two) {
      if (idx_one + idx_two < 4 || idx_one + idx_two > 12) {
        continue;
      }
      size_t const idx_three = 12 - (idx_one + idx_two);
      real_t const w3j = WIGNER4.table.at(idx_one).at(idx_two);
      if (w3j == static_cast<real_t>(0.0)) {
        continue;
      }

      std::complex<real_t> const prod = q4m.at(idx_one) * q4m.at(idx_two) * q4m.at(idx_three);
      w4_val += w3j * prod.real();
    }
  }
  return w4_val;
}

real_t computeW6(const std::array<std::complex<real_t>, 13> &q6m) {
  static const Wigner6Table WIGNER6;
  auto w6_val = static_cast<real_t>(0.0);
  for (size_t idx_one = 0; idx_one < 13; ++idx_one) {
    for (size_t idx_two = 0; idx_two < 13; ++idx_two) {
      if (idx_one + idx_two < 6 || idx_one + idx_two > 18) {
        continue;
      }
      size_t const idx_three = 18 - (idx_one + idx_two);
      real_t const w3j = WIGNER6.table.at(idx_one).at(idx_two);
      if (w3j == static_cast<real_t>(0.0)) {
        continue;
      }

      std::complex<real_t> const prod = q6m.at(idx_one) * q6m.at(idx_two) * q6m.at(idx_three);
      w6_val += w3j * prod.real();
    }
  }
  return w6_val;
}

struct AtomHarmonicVectors {
  std::array<std::complex<real_t>, 9> q4m{};
  std::array<std::complex<real_t>, 13> q6m{};
  bool has_neighbors{false};
};

AtomHarmonicVectors computeHarmonicVectors(size_t atom_idx,
                                           const correlation::core::NeighborGraph &neighbor_graph) {
  const auto &atom_neighbors = neighbor_graph.getNeighbors(atom_idx);
  size_t const num_neighbors = atom_neighbors.size();
  if (num_neighbors < 2) {
    return {};
  }

  AtomHarmonicVectors vecs;
  vecs.has_neighbors = true;

  for (const auto &neighbor : atom_neighbors) {
    correlation::math::Vector3<real_t> const r_ij = neighbor.r_ij;
    real_t const distance = neighbor.distance;
    if (distance == static_cast<real_t>(0.0)) {
      continue;
    }

    real_t const theta = std::acos(
        std::clamp(r_ij.z() / distance, static_cast<real_t>(-1.0), static_cast<real_t>(1.0)));
    real_t const phi = std::atan2(r_ij.y(), r_ij.x());

    for (size_t m_idx = 0; m_idx < 9; ++m_idx) {
      int const m_val = static_cast<int>(m_idx) - 4;
      vecs.q4m.at(m_idx) += SteinhardtCalculator::sphericalHarmonic(4, m_val,
                                                                    {
                                                                        .theta = theta,
                                                                        .phi = phi,
                                                                    });
    }
    for (size_t m_idx = 0; m_idx < 13; ++m_idx) {
      int const m_val = static_cast<int>(m_idx) - 6;
      vecs.q6m.at(m_idx) += SteinhardtCalculator::sphericalHarmonic(6, m_val,
                                                                    {
                                                                        .theta = theta,
                                                                        .phi = phi,
                                                                    });
    }
  }

  for (auto &val : vecs.q4m) {
    val /= static_cast<real_t>(num_neighbors);
  }
  for (auto &val : vecs.q6m) {
    val /= static_cast<real_t>(num_neighbors);
  }
  return vecs;
}

SingleAtomSteinhardt computeSingleAtomInvariants(const AtomHarmonicVectors &vecs,
                                                 GlobalSteinhardtFactors factors) {
  if (!vecs.has_neighbors) {
    return {};
  }

  auto sum_sq_4 = static_cast<real_t>(0.0);
  for (const auto &val : vecs.q4m) {
    sum_sq_4 += std::norm(val);
  }
  real_t const q4_val = factors.global_q4_factor * std::sqrt(sum_sq_4);

  auto sum_sq_6 = static_cast<real_t>(0.0);
  for (const auto &val : vecs.q6m) {
    sum_sq_6 += std::norm(val);
  }
  real_t const q6_val = factors.global_q6_factor * std::sqrt(sum_sq_6);

  real_t const w4_val = computeW4(vecs.q4m);
  auto w4_hat_val = static_cast<real_t>(0.0);
  if (sum_sq_4 > static_cast<real_t>(1e-12)) {
    w4_hat_val = w4_val / std::pow(sum_sq_4, static_cast<real_t>(1.5));
  }

  real_t const w6_val = computeW6(vecs.q6m);
  auto w6_hat_val = static_cast<real_t>(0.0);
  if (sum_sq_6 > static_cast<real_t>(1e-12)) {
    w6_hat_val = w6_val / std::pow(sum_sq_6, static_cast<real_t>(1.5));
  }

  return {
      .Q4 = q4_val,
      .Q6 = q6_val,
      .W4_hat = w4_hat_val,
      .W6_hat = w6_hat_val,
      .Q4_bar = static_cast<real_t>(0.0),
      .Q6_bar = static_cast<real_t>(0.0),
  };
}

std::pair<real_t, real_t> computeLechnerDellagoAveraged(
    size_t atom_idx, const correlation::core::NeighborGraph &neighbor_graph,
    const std::vector<AtomHarmonicVectors> &all_vecs, GlobalSteinhardtFactors factors) {
  const auto &atom_neighbors = neighbor_graph.getNeighbors(atom_idx);
  if (atom_neighbors.size() < 2 || !all_vecs[atom_idx].has_neighbors) {
    return {static_cast<real_t>(0.0), static_cast<real_t>(0.0)};
  }

  std::array<std::complex<real_t>, 9> q4m_bar{};
  std::array<std::complex<real_t>, 13> q6m_bar{};

  for (size_t m_idx = 0; m_idx < 9; ++m_idx) {
    q4m_bar.at(m_idx) = all_vecs[atom_idx].q4m.at(m_idx);
  }
  for (size_t m_idx = 0; m_idx < 13; ++m_idx) {
    q6m_bar.at(m_idx) = all_vecs[atom_idx].q6m.at(m_idx);
  }

  for (const auto &neighbor : atom_neighbors) {
    size_t const neighbor_idx = neighbor.index;
    for (size_t m_idx = 0; m_idx < 9; ++m_idx) {
      q4m_bar.at(m_idx) += all_vecs[neighbor_idx].q4m.at(m_idx);
    }
    for (size_t m_idx = 0; m_idx < 13; ++m_idx) {
      q6m_bar.at(m_idx) += all_vecs[neighbor_idx].q6m.at(m_idx);
    }
  }

  real_t const inv_count =
      static_cast<real_t>(1.0) / static_cast<real_t>(atom_neighbors.size() + 1);

  auto sum_sq_4_bar = static_cast<real_t>(0.0);
  for (size_t m_idx = 0; m_idx < 9; ++m_idx) {
    q4m_bar.at(m_idx) *= inv_count;
    sum_sq_4_bar += std::norm(q4m_bar.at(m_idx));
  }
  real_t const q4_bar = factors.global_q4_factor * std::sqrt(sum_sq_4_bar);

  auto sum_sq_6_bar = static_cast<real_t>(0.0);
  for (size_t m_idx = 0; m_idx < 13; ++m_idx) {
    q6m_bar.at(m_idx) *= inv_count;
    sum_sq_6_bar += std::norm(q6m_bar.at(m_idx));
  }
  real_t const q6_bar = factors.global_q6_factor * std::sqrt(sum_sq_6_bar);

  return {q4_bar, q6_bar};
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

void addValueToHistogram(std::map<std::string, std::vector<real_t>> &partials,
                         const std::string &symbol, real_t val, BinningConfig config) {
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

struct SteinhardtHistograms {
  correlation::analysis::Histogram *hist_q4{nullptr};
  correlation::analysis::Histogram *hist_q6{nullptr};
  correlation::analysis::Histogram *hist_w4{nullptr};
  correlation::analysis::Histogram *hist_w6{nullptr};
  correlation::analysis::Histogram *hist_q4_bar{nullptr};
  correlation::analysis::Histogram *hist_q6_bar{nullptr};
};

struct ThreadLocalHist {
  std::map<std::string, std::vector<real_t>> partials_q4;
  std::map<std::string, std::vector<real_t>> partials_q6;
  std::map<std::string, std::vector<real_t>> partials_w4;
  std::map<std::string, std::vector<real_t>> partials_w6;
  std::map<std::string, std::vector<real_t>> partials_q4_bar;
  std::map<std::string, std::vector<real_t>> partials_q6_bar;
  real_t num_atoms_f{0.0};
};

void initAllPartials(ThreadLocalHist &local, const std::vector<std::string> &symbols,
                     const HistogramConfigs &configs) {
  initHistogramMap(local.partials_q4, symbols, configs.bins_q);
  initHistogramMap(local.partials_q6, symbols, configs.bins_q);
  initHistogramMap(local.partials_w4, symbols, configs.bins_w4);
  initHistogramMap(local.partials_w6, symbols, configs.bins_w6);
  initHistogramMap(local.partials_q4_bar, symbols, configs.bins_q);
  initHistogramMap(local.partials_q6_bar, symbols, configs.bins_q);
}

void accumulateAllPartials(ThreadLocalHist &dest, const ThreadLocalHist &src) {
  accumulateHistogramMap(dest.partials_q4, src.partials_q4);
  accumulateHistogramMap(dest.partials_q6, src.partials_q6);
  accumulateHistogramMap(dest.partials_w4, src.partials_w4);
  accumulateHistogramMap(dest.partials_w6, src.partials_w6);
  accumulateHistogramMap(dest.partials_q4_bar, src.partials_q4_bar);
  accumulateHistogramMap(dest.partials_q6_bar, src.partials_q6_bar);
}

void copyAllPartials(const SteinhardtHistograms &hists, const ThreadLocalHist &accum) {
  if (hists.hist_q4 != nullptr) {
    copyPartialsToHistogram(*hists.hist_q4, accum.partials_q4);
  }
  if (hists.hist_q6 != nullptr) {
    copyPartialsToHistogram(*hists.hist_q6, accum.partials_q6);
  }
  if (hists.hist_w4 != nullptr) {
    copyPartialsToHistogram(*hists.hist_w4, accum.partials_w4);
  }
  if (hists.hist_w6 != nullptr) {
    copyPartialsToHistogram(*hists.hist_w6, accum.partials_w6);
  }
  if (hists.hist_q4_bar != nullptr) {
    copyPartialsToHistogram(*hists.hist_q4_bar, accum.partials_q4_bar);
  }
  if (hists.hist_q6_bar != nullptr) {
    copyPartialsToHistogram(*hists.hist_q6_bar, accum.partials_q6_bar);
  }
}

void populateHistograms(const correlation::core::Cell &cell,
                        const correlation::analysis::StructureAnalyzer *neighbors,
                        const SteinhardtParams &params, HistogramConfigs configs,
                        SteinhardtHistograms &hists) {
  const auto &atoms = cell.atoms();
  const auto &neighbor_graph = neighbors->neighborGraph();
  size_t const num_atoms = atoms.size();

  // Collect element symbols for pre-sizing thread-local maps
  std::vector<std::string> element_symbols;
  for (const auto &elem : cell.elements()) {
    element_symbols.push_back(elem.symbol);
  }

  BinningConfig const config_q{
      .min_val = static_cast<real_t>(0.0),
      .max_val = configs.q_max,
      .d_val = configs.step_q,
  };
  BinningConfig const config_w4{
      .min_val = configs.w4_min,
      .max_val = configs.w4_max,
      .d_val = configs.step_w4,
  };
  BinningConfig const config_w6{
      .min_val = configs.w6_min,
      .max_val = configs.w6_max,
      .d_val = configs.step_w6,
  };

  tbb::enumerable_thread_specific<ThreadLocalHist> ets([&]() {
    ThreadLocalHist local;
    initAllPartials(local, element_symbols, configs);
    return local;
  });

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, num_atoms), [&](const tbb::blocked_range<size_t> &range) {
        auto &local = ets.local();

        for (size_t i = range.begin(); i != range.end(); ++i) {
          if (neighbor_graph.getNeighbors(i).size() < 2) {
            continue;
          }

          local.num_atoms_f += static_cast<real_t>(1.0);
          const std::string &symbol = atoms[i].element().symbol;

          addValueToHistogram(local.partials_q4, symbol, params.Q4[i], config_q);
          addValueToHistogram(local.partials_q6, symbol, params.Q6[i], config_q);
          addValueToHistogram(local.partials_w4, symbol, params.W4_hat[i], config_w4);
          addValueToHistogram(local.partials_w6, symbol, params.W6_hat[i], config_w6);
          addValueToHistogram(local.partials_q4_bar, symbol, params.Q4_bar[i], config_q);
          addValueToHistogram(local.partials_q6_bar, symbol, params.Q6_bar[i], config_q);
        }
      });

  // Reduce thread-local histograms
  ThreadLocalHist reduced;
  initAllPartials(reduced, element_symbols, configs);

  auto num_atoms_f = static_cast<real_t>(0.0);
  for (const auto &local : ets) {
    num_atoms_f += local.num_atoms_f;
    accumulateAllPartials(reduced, local);
  }

  if (num_atoms_f > static_cast<real_t>(0.0)) {
    normalizeHistogramMap(reduced.partials_q4, num_atoms_f * configs.step_q);
    normalizeHistogramMap(reduced.partials_q6, num_atoms_f * configs.step_q);
    normalizeHistogramMap(reduced.partials_w4, num_atoms_f * configs.step_w4);
    normalizeHistogramMap(reduced.partials_w6, num_atoms_f * configs.step_w6);
    normalizeHistogramMap(reduced.partials_q4_bar, num_atoms_f * configs.step_q);
    normalizeHistogramMap(reduced.partials_q6_bar, num_atoms_f * configs.step_q);
  }

  copyAllPartials(hists, reduced);
}
} // namespace

std::complex<real_t> SteinhardtCalculator::sphericalHarmonic(int degree, int order,
                                                             SphericalAngles angles) {
  if (order >= 0) {
    real_t const p_lm = correlation::math::sphLegendre(
        {
            .degree = degree,
            .order = order,
        },
        angles.theta);
    return p_lm * std::polar(static_cast<real_t>(1.0), static_cast<real_t>(order) * angles.phi);
  } // For negative m: Y_l^{-m} = (-1)^m (Y_l^m)*
  int const abs_m = -order;
  real_t const p_lm = correlation::math::sphLegendre(
      {
          .degree = degree,
          .order = abs_m,
      },
      angles.theta);
  std::complex<real_t> const y_l_m =
      p_lm * std::polar(static_cast<real_t>(1.0), static_cast<real_t>(abs_m) * angles.phi);
  std::complex<real_t> y_l_minus_m = std::conj(y_l_m);
  if (abs_m % 2 != 0) {
    y_l_minus_m = -y_l_minus_m;
  }
  return y_l_minus_m;
}

real_t SteinhardtCalculator::wigner3j(Wigner3jParams params) {
  int const j_one = params.j_one;
  int const j_two = params.j_two;
  int const j_three = params.j_three;
  int const m_one = params.m_one;
  int const m_two = params.m_two;
  int const m_three = params.m_three;

  if (m_one + m_two + m_three != 0) {
    return 0.0;
  }
  if (j_three < std::abs(j_one - j_two) || j_three > j_one + j_two) {
    return 0.0;
  }
  if (std::abs(m_one) > j_one || std::abs(m_two) > j_two || std::abs(m_three) > j_three) {
    return 0.0;
  }

  real_t delta = (correlation::math::factorial(j_one + j_two - j_three) *
                  correlation::math::factorial(j_one - j_two + j_three) *
                  correlation::math::factorial(-j_one + j_two + j_three) /
                  correlation::math::factorial(j_one + j_two + j_three + 1));
  delta = std::sqrt(delta);

  real_t comp =
      (correlation::math::factorial(j_one - m_one) * correlation::math::factorial(j_one + m_one) *
       correlation::math::factorial(j_two - m_two) * correlation::math::factorial(j_two + m_two) *
       correlation::math::factorial(j_three - m_three) *
       correlation::math::factorial(j_three + m_three));
  comp = std::sqrt(comp);

  real_t const phase1 = ((j_one - j_two - m_three) % 2 != 0) ? -1.0 : 1.0;

  int const k_min = std::max({0, j_two - j_three - m_one, j_one - j_three + m_two});
  int const k_max = std::min({j_one + j_two - j_three, j_one - m_one, j_two + m_two});

  real_t sum = 0.0;
  for (int k = k_min; k <= k_max; ++k) {
    real_t const k_phase = (k % 2 != 0) ? -1.0 : 1.0;
    real_t const denom = (correlation::math::factorial(k) *
                          correlation::math::factorial(j_one + j_two - j_three - k) *
                          correlation::math::factorial(j_one - m_one - k) *
                          correlation::math::factorial(j_two + m_two - k) *
                          correlation::math::factorial(j_three - j_two + m_one + k) *
                          correlation::math::factorial(j_three - j_one - m_two + k));
    sum += k_phase / denom;
  }

  return phase1 * delta * comp * sum;
}

void SteinhardtCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &dists,
    const correlation::analysis::AnalysisSettings &settings) const {
#if defined(CORRELATION_USE_CUDA) || defined(CORRELATION_USE_HIP)
  static const GPUSteinhardtCalculator GPU_CALC;
  if (GPU_CALC.hasGPU() && dists.neighbors() != nullptr) {
    GPU_CALC.calculateFrame(dists, settings);
    return;
  }
#else
  (void)settings;
#endif
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

  // 1. Compute Q4, Q6, W4_hat, W6_hat for all atoms (Pass 1)
  SteinhardtParams params;
  params.Q4.resize(num_atoms, static_cast<real_t>(0.0));
  params.Q6.resize(num_atoms, static_cast<real_t>(0.0));
  params.W4_hat.resize(num_atoms, static_cast<real_t>(0.0));
  params.W6_hat.resize(num_atoms, static_cast<real_t>(0.0));
  params.Q4_bar.resize(num_atoms, static_cast<real_t>(0.0));
  params.Q6_bar.resize(num_atoms, static_cast<real_t>(0.0));

  std::vector<AtomHarmonicVectors> all_vecs(num_atoms);

  const auto global_q4_factor = static_cast<real_t>(std::sqrt(correlation::math::four_pi / 9.0));
  const auto global_q6_factor = static_cast<real_t>(std::sqrt(correlation::math::four_pi / 13.0));
  GlobalSteinhardtFactors const factors{
      .global_q4_factor = global_q4_factor,
      .global_q6_factor = global_q6_factor,
  };

  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_atoms),
                    [&](const tbb::blocked_range<size_t> &range) {
                      for (size_t i = range.begin(); i != range.end(); ++i) {
                        all_vecs[i] = computeHarmonicVectors(i, neighbor_graph);
                        auto const res = computeSingleAtomInvariants(all_vecs[i], factors);
                        params.Q4[i] = res.Q4;
                        params.Q6[i] = res.Q6;
                        params.W4_hat[i] = res.W4_hat;
                        params.W6_hat[i] = res.W6_hat;
                      }
                    });

  // 2. Lechner-Dellago local averaging (Pass 2)
  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_atoms),
                    [&](const tbb::blocked_range<size_t> &range) {
                      for (size_t i = range.begin(); i != range.end(); ++i) {
                        auto const [q4_bar, q6_bar] =
                            computeLechnerDellagoAveraged(i, neighbor_graph, all_vecs, factors);
                        params.Q4_bar[i] = q4_bar;
                        params.Q6_bar[i] = q6_bar;
                      }
                    });

  // 3. Initialize Histograms
  size_t const bins_q = 100;
  real_t const q_max = 1.0;
  real_t const d_q = q_max / bins_q;

  size_t const bins_w4 = 100;
  real_t const w4_min = -0.5;
  real_t const w4_max = 0.5;
  real_t const d_w4 = (w4_max - w4_min) / bins_w4;

  size_t const bins_w6 = 100;
  real_t const w6_min = -0.2;
  real_t const w6_max = 0.2;
  real_t const d_w6 = (w6_max - w6_min) / bins_w6;

  correlation::analysis::Histogram hist_q4;
  hist_q4.x_label = "Q4";
  hist_q4.title = "Steinhardt Q4 Interface Parameter";
  hist_q4.y_label = "Probability";
  hist_q4.x_unit = "arbitrary units";
  hist_q4.y_unit = "counts";
  hist_q4.description = "Steinhardt Q4 Bond Orientational Order Parameter";
  hist_q4.file_suffix = "_Q4";
  hist_q4.bins.resize(bins_q);
  for (size_t bin_idx = 0; bin_idx < bins_q; ++bin_idx) {
    hist_q4.bins[bin_idx] = static_cast<real_t>(static_cast<real_t>(bin_idx) + 0.5 * d_q);
  }

  correlation::analysis::Histogram hist_q6;
  hist_q6.x_label = "Q6";
  hist_q6.title = "Steinhardt Q6 Interface Parameter";
  hist_q6.y_label = "Probability";
  hist_q6.x_unit = "arbitrary units";
  hist_q6.y_unit = "counts";
  hist_q6.description = "Steinhardt Q6 Bond Orientational Order Parameter";
  hist_q6.file_suffix = "_Q6";
  hist_q6.bins.resize(bins_q);
  for (size_t bin_idx = 0; bin_idx < bins_q; ++bin_idx) {
    hist_q6.bins[bin_idx] = static_cast<real_t>(static_cast<real_t>(bin_idx) + 0.5 * d_q);
  }

  correlation::analysis::Histogram hist_w4;
  hist_w4.x_label = "W4_hat";
  hist_w4.title = "Steinhardt W4_hat Parameter";
  hist_w4.y_label = "Probability";
  hist_w4.x_unit = "arbitrary units";
  hist_w4.y_unit = "counts";
  hist_w4.description = "Steinhardt Normalized W4 Bond Orientational Order Parameter";
  hist_w4.file_suffix = "_W4_hat";
  hist_w4.bins.resize(bins_w4);
  for (size_t bin_idx = 0; bin_idx < bins_w4; ++bin_idx) {
    hist_w4.bins[bin_idx] = w4_min + static_cast<real_t>(static_cast<real_t>(bin_idx) + 0.5 * d_w4);
  }

  correlation::analysis::Histogram hist_w6;
  hist_w6.x_label = "W6_hat";
  hist_w6.title = "Steinhardt W6_hat Parameter";
  hist_w6.y_label = "Probability";
  hist_w6.x_unit = "arbitrary units";
  hist_w6.y_unit = "counts";
  hist_w6.description = "Steinhardt Normalized W6 Bond Orientational Order Parameter";
  hist_w6.file_suffix = "_W6_hat";
  hist_w6.bins.resize(bins_w6);
  for (size_t bin_idx = 0; bin_idx < bins_w6; ++bin_idx) {
    hist_w6.bins[bin_idx] = w6_min + static_cast<real_t>(static_cast<real_t>(bin_idx) + 0.5 * d_w6);
  }

  correlation::analysis::Histogram hist_q4_bar;
  hist_q4_bar.x_label = "Q4_bar";
  hist_q4_bar.title = "Lechner-Dellago Averaged Q4_bar Parameter";
  hist_q4_bar.y_label = "Probability";
  hist_q4_bar.x_unit = "arbitrary units";
  hist_q4_bar.y_unit = "counts";
  hist_q4_bar.description = "Lechner-Dellago Locally Averaged Q4 Order Parameter";
  hist_q4_bar.file_suffix = "_Q4_bar";
  hist_q4_bar.bins.resize(bins_q);
  for (size_t bin_idx = 0; bin_idx < bins_q; ++bin_idx) {
    hist_q4_bar.bins[bin_idx] = static_cast<real_t>(static_cast<real_t>(bin_idx) + 0.5 * d_q);
  }

  correlation::analysis::Histogram hist_q6_bar;
  hist_q6_bar.x_label = "Q6_bar";
  hist_q6_bar.title = "Lechner-Dellago Averaged Q6_bar Parameter";
  hist_q6_bar.y_label = "Probability";
  hist_q6_bar.x_unit = "arbitrary units";
  hist_q6_bar.y_unit = "counts";
  hist_q6_bar.description = "Lechner-Dellago Locally Averaged Q6 Order Parameter";
  hist_q6_bar.file_suffix = "_Q6_bar";
  hist_q6_bar.bins.resize(bins_q);
  for (size_t bin_idx = 0; bin_idx < bins_q; ++bin_idx) {
    hist_q6_bar.bins[bin_idx] = static_cast<real_t>(static_cast<real_t>(bin_idx) + 0.5 * d_q);
  }

  // 4. Populate Histograms
  HistogramConfigs const configs{
      .bins_q = bins_q,
      .q_max = q_max,
      .step_q = d_q,
      .bins_w4 = bins_w4,
      .w4_min = w4_min,
      .w4_max = w4_max,
      .step_w4 = d_w4,
      .bins_w6 = bins_w6,
      .w6_min = w6_min,
      .w6_max = w6_max,
      .step_w6 = d_w6,
  };
  SteinhardtHistograms hists_bundle{
      .hist_q4 = &hist_q4,
      .hist_q6 = &hist_q6,
      .hist_w4 = &hist_w4,
      .hist_w6 = &hist_w6,
      .hist_q4_bar = &hist_q4_bar,
      .hist_q6_bar = &hist_q6_bar,
  };
  populateHistograms(cell, neighbors, params, configs, hists_bundle);

  std::map<std::string, correlation::analysis::Histogram> hists;
  hists["Q4"] = std::move(hist_q4);
  hists["Q6"] = std::move(hist_q6);
  hists["W4_hat"] = std::move(hist_w4);
  hists["W6_hat"] = std::move(hist_w6);
  hists["Q4_bar"] = std::move(hist_q4_bar);
  hists["Q6_bar"] = std::move(hist_q6_bar);

  return hists;
}

} // namespace correlation::calculators
