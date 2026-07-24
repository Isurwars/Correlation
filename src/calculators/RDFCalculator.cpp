/**
 * @file RDFCalculator.cpp
 * @brief Implementation of the radial distribution function calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/RDFCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include "math/Precision.hpp"
#include "math/SIMDUtils.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace correlation::calculators {

namespace {
std::string getPartialKey(const correlation::core::Cell &cell, size_t type1, size_t type2) {
  const auto &elements = cell.elements();
  if (type1 > type2) {
    std::swap(type1, type2);
  }
  return elements[type1].symbol + "-" + elements[type2].symbol;
}

std::string getInversePartialKey(const correlation::core::Cell &cell, size_t type1, size_t type2) {
  const auto &elements = cell.elements();
  if (type1 < type2) {
    std::swap(type1, type2);
  }
  return elements[type1].symbol + "-" + elements[type2].symbol;
}

// Static registration of the calculator in the factory
const bool registered = CalculatorFactory::registerTypeSafe<RDFCalculator>("RDFCalculator");

struct RDFSettings {
  real_t r_max;
  real_t r_bin_width;
  size_t num_bins;
};

void accumulateRawCounts(const correlation::core::Cell &cell, const correlation::analysis::StructureAnalyzer *neighbors,
                         RDFSettings settings, correlation::analysis::Histogram &H_r) {
  const auto &elements = cell.elements();
  const size_t num_elements = elements.size();

  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      std::string const key = getPartialKey(cell, i, j);
      auto &partial_hist = H_r.partials[key];
      partial_hist.assign(settings.num_bins, 0.0);

      // First Pass: Accumulate the raw distance counts H(r) into histogram
      // bins. Distances are pre-computed in
      // `StructureAnalyzer::distance_tensor_`.
      for (const auto &dist : neighbors->distances()[i][j]) {
        if (dist < settings.r_max) {
          auto const bin = static_cast<size_t>(dist / settings.r_bin_width);
          if (bin < settings.num_bins) {
            partial_hist[bin] += 1.0;
          }
        }
      }

      // For self-pairs (A-A), each pair is counted once in the upper triangular
      // DistanceTensor. We multiply by 2 to account for both A_1 -> A_2 and A_2
      // -> A_1 interactions.
      if (i == j) {
        correlation::math::scale_bins(partial_hist.data(), static_cast<real_t>(2.0), settings.num_bins);
      }
    }
  }
}

struct RDFNormalizationSettings {
  real_t volume;
  real_t bin_width;
  size_t num_bins;
};

void normalizeDistributions(const correlation::core::Cell &cell, const std::map<std::string, real_t> &element_counts,
                            RDFNormalizationSettings settings, const correlation::analysis::Histogram &H_r,
                            correlation::analysis::Histogram &g_r, correlation::analysis::Histogram &G_r,
                            correlation::analysis::Histogram &J_r) {
  const auto &elements = cell.elements();
  const size_t num_elements = elements.size();

  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      std::string const key = getPartialKey(cell, i, j);
      std::string const inversekey = getInversePartialKey(cell, i, j);

      const std::string &sym_i = elements[i].symbol;
      const std::string &sym_j = elements[j].symbol;

      const real_t N_i = element_counts.at(sym_i);
      const real_t N_j = element_counts.at(sym_j);

      const auto &H_ij = H_r.partials.at(key);

      g_r.partials[key].assign(settings.num_bins, 0.0);
      G_r.partials[key].assign(settings.num_bins, 0.0);
      J_r.partials[key].assign(settings.num_bins, 0.0);
      J_r.partials[inversekey].assign(settings.num_bins, 0.0);

      // Second Pass: Normalize the raw counts H(r) into target distribution
      // functions. g(r) normalization constant: V / (4 * pi * dr * N_i * N_j).
      // The r^2 term is applied per-bin inside the SIMD kernel.
      const real_t g_norm_constant = settings.volume / (correlation::math::four_pi * settings.bin_width * N_i * N_j);
      const real_t rho_j = N_j / settings.volume;
      const real_t inv_Ni_dr = static_cast<real_t>(1.0) / (N_i * settings.bin_width);
      const real_t inv_Nj_dr = static_cast<real_t>(1.0) / (N_j * settings.bin_width);
      const real_t pi4_rho_j = correlation::math::four_pi * rho_j;

      correlation::math::RDFNormalizationParams<real_t> params{
          .hist_data = H_ij.data(),
          .radial_bins = g_r.bins.data(),
          .g_norm = g_norm_constant,
          .inv_Ni_dr = inv_Ni_dr,
          .inv_Nj_dr = inv_Nj_dr,
          .pi4_rho_j = pi4_rho_j,
          .g_out = g_r.partials[key].data(),
          .G_out = G_r.partials[key].data(),
          .J_out = J_r.partials[key].data(),
          .Jinv_out = J_r.partials[inversekey].data(),
          .count = settings.num_bins,
      };
      correlation::math::normalize_rdf_bins(params);
    }
  }
}

struct RDFWeightingSettings {
  real_t rho_0;
  size_t num_bins;
};

void weightPartials(const correlation::core::Cell &cell, const std::map<std::string, real_t> &ashcroft_weights,
                    RDFWeightingSettings settings, correlation::analysis::Histogram &g_r,
                    correlation::analysis::Histogram &G_r) {
  const auto &elements = cell.elements();
  const size_t num_elements = elements.size();

  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      std::string const key = getPartialKey(cell, i, j);
      real_t const weight = ashcroft_weights.at(key);

      // 1. Weight g_r partial: g_ij_weighted(r) = w_ij * g_ij(r)
      auto &g_part = g_r.partials.at(key);
      for (size_t k = 0; k < settings.num_bins; ++k) {
        g_part[k] *= weight;
      }

      // 2. Weight G_r partial: G_ij_weighted(r) = w_ij * 4 * pi * rho_0 * r * (g_ij(r) - 1)
      // Since we already weighted g_part, we have g_part[k] = w_ij * g_ij(r).
      // Thus, G_ij_weighted(r) = 4 * pi * rho_0 * r * (g_part[k] - weight)
      auto &G_part = G_r.partials.at(key);
      for (size_t k = 0; k < settings.num_bins; ++k) {
        const real_t r_k = g_r.bins[k];
        if (r_k < 1e-9) {
          G_part[k] = 0.0;
        } else {
          G_part[k] = correlation::math::four_pi * settings.rho_0 * r_k * (g_part[k] - weight);
        }
      }
    }
  }
}
} // namespace

void RDFCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                   const correlation::analysis::AnalysisSettings &settings) const {
  auto results =
      calculate(dists.cell(), dists.neighbors(), dists.getAshcroftWeights(), settings.r_max, settings.r_bin_width);
  for (auto &[name, histogram] : results) {
    dists.addHistogram(name, std::move(histogram));
  }
}

std::map<std::string, correlation::analysis::Histogram>
RDFCalculator::calculate(const correlation::core::Cell &cell, const correlation::analysis::StructureAnalyzer *neighbors,
                         const std::map<std::string, real_t> &ashcroft_weights, real_t r_max, real_t r_bin_width) {
  if (r_bin_width <= 0) {
    throw std::invalid_argument("Bin width must be positive, got: " + std::to_string(r_bin_width));
  }
  if (r_max <= 0) {
    throw std::invalid_argument("Cutoff radius must be positive, got: " + std::to_string(r_max));
  }

  const real_t volume = cell.volume();
  if (volume <= std::numeric_limits<real_t>::epsilon()) {
    throw std::logic_error("Cell volume must be positive, got: " + std::to_string(volume));
  }

  const auto &elements = cell.elements();
  const size_t num_elements = elements.size();
  const auto num_atoms = static_cast<real_t>(cell.atomCount());
  if (num_atoms == 0.0) {
    return {};
  }

  std::map<std::string, real_t> element_counts;
  for (const auto &atom : cell.atoms()) {
    element_counts[atom.element().symbol]++;
  }

  const auto num_bins = static_cast<size_t>(std::floor(r_max / r_bin_width));
  const real_t Vol = static_cast<real_t>(cell.volume());
  const real_t d_r = r_bin_width;
  const real_t rho_0 = num_atoms / Vol;

  correlation::analysis::Histogram H_r;
  correlation::analysis::Histogram g_r;
  correlation::analysis::Histogram G_r;
  correlation::analysis::Histogram J_r;
  H_r.bins.resize(num_bins);
  g_r.bins.resize(num_bins);
  G_r.bins.resize(num_bins);
  J_r.bins.resize(num_bins);
  H_r.x_label = "r";
  H_r.title = "H(r) — Distance Histogram";
  H_r.y_label = "H(r)";
  H_r.x_unit = "Å";
  H_r.y_unit = "counts";
  H_r.description = "Radial Distribution Function";
  H_r.file_suffix = "_H";

  g_r.x_label = "r";
  g_r.title = "g(r) — Radial Distribution Function";
  g_r.y_label = "g(r)";
  g_r.x_unit = "Å";
  g_r.y_unit = "dimensionless";
  g_r.description = "Radial Distribution Function";
  g_r.file_suffix = "_g";

  G_r.x_label = "r";
  G_r.title = "G(r) — Reduced Pair Distribution Function";
  G_r.y_label = "G(r)";
  G_r.x_unit = "Å";
  G_r.y_unit = "Å⁻²";
  G_r.description = "Radial Distribution Function";
  G_r.file_suffix = "_G";

  J_r.x_label = "r";
  J_r.title = "J(r) — Reduced Pair Distribution";
  J_r.y_label = "J(r)";
  J_r.x_unit = "Å";
  J_r.y_unit = "Å⁻¹";
  J_r.description = "Radial Distribution Function";
  J_r.file_suffix = "_J";

  for (size_t i = 0; i < num_bins; ++i) {
    const real_t r_i = (static_cast<real_t>(i) + static_cast<real_t>(0.5)) * r_bin_width;
    H_r.bins[i] = r_i;
    g_r.bins[i] = r_i;
    G_r.bins[i] = r_i;
    J_r.bins[i] = r_i;
  }

  accumulateRawCounts(cell, neighbors,
                      {
                          .r_max = r_max,
                          .r_bin_width = r_bin_width,
                          .num_bins = num_bins,
                      },
                      H_r);

  normalizeDistributions(cell, element_counts,
                         {
                             .volume = Vol,
                             .bin_width = d_r,
                             .num_bins = num_bins,
                         },
                         H_r, g_r, G_r, J_r);

  auto &total_g = g_r.partials["Total"];
  total_g.assign(num_bins, 0.0);
  for (const auto &[key, g_partial] : g_r.partials) {
    if (key == "Total") {
      continue;
    }

    real_t const weight = ashcroft_weights.at(key);

    for (size_t k = 0; k < num_bins; ++k) {
      total_g[k] += g_partial[k] * weight;
    }
  }

  auto &total_J = J_r.partials["Total"];
  auto &total_G = G_r.partials["Total"];
  total_J.assign(num_bins, 0.0);
  total_G.assign(num_bins, 0.0);

  for (size_t k = 0; k < num_bins; ++k) {
    const real_t r_k = g_r.bins[k];
    if (r_k < 1e-9) {
      continue;
    }

    total_J[k] = correlation::math::four_pi * r_k * r_k * rho_0 * total_g[k];
    total_G[k] = correlation::math::four_pi * rho_0 * r_k * (total_g[k] - static_cast<real_t>(1.0));
  }

  weightPartials(cell, ashcroft_weights,
                 {
                     .rho_0 = rho_0,
                     .num_bins = num_bins,
                 },
                 g_r, G_r);

  std::map<std::string, correlation::analysis::Histogram> results;
  results["J_r"] = std::move(J_r);
  results["g_r"] = std::move(g_r);
  results["G_r"] = std::move(G_r);

  return results;
}

} // namespace correlation::calculators
