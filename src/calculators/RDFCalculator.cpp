/**
 * @file RDFCalculator.cpp
 * @brief Implementation of the radial distribution function calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "calculators/RDFCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include "math/SIMDUtils.hpp"

#include <cmath>
#include <stdexcept>

namespace {
std::string getPartialKey(const Cell &cell, int type1, int type2) {
  const auto &elements = cell.elements();
  if (type1 > type2)
    std::swap(type1, type2);
  return elements[type1].symbol + "-" + elements[type2].symbol;
}

std::string getInversePartialKey(const Cell &cell, int type1, int type2) {
  const auto &elements = cell.elements();
  if (type1 < type2)
    std::swap(type1, type2);
  return elements[type1].symbol + "-" + elements[type2].symbol;
}

// Static registration
bool registered = CalculatorFactory::instance().registerCalculator(
    std::make_unique<RDFCalculator>());
} // namespace

void RDFCalculator::calculateFrame(DistributionFunctions &df,
                                   const AnalysisSettings &settings) const {
  auto results = calculate(df.cell(), df.neighbors(), df.getAshcroftWeights(),
                           settings.r_max, settings.r_bin_width);
  for (auto &[name, histogram] : results) {
    df.addHistogram(name, std::move(histogram));
  }
}

std::map<std::string, Histogram>
RDFCalculator::calculate(const Cell &cell, const StructureAnalyzer *neighbors,
                         const std::map<std::string, double> &ashcroft_weights,
                         double r_max, double r_bin_width) {
  if (r_bin_width <= 0) {
    throw std::invalid_argument("Bin width must be positive, got: " +
                                std::to_string(r_bin_width));
  }
  if (r_max <= 0) {
    throw std::invalid_argument("Cutoff radius must be positive, got: " +
                                std::to_string(r_max));
  }

  const double volume = cell.volume();
  if (volume <= std::numeric_limits<double>::epsilon()) {
    throw std::logic_error("Cell volume must be positive, got: " +
                           std::to_string(volume));
  }

  const auto &elements = cell.elements();
  const size_t num_elements = elements.size();
  const double num_atoms = static_cast<double>(cell.atomCount());
  if (num_atoms == 0.0)
    return {};

  std::map<std::string, double> element_counts;
  for (const auto &atom : cell.atoms()) {
    element_counts[atom.element().symbol]++;
  }

  const size_t num_bins = static_cast<size_t>(std::floor(r_max / r_bin_width));
  const double V = cell.volume();
  const double dr = r_bin_width;
  const double rho_0 = num_atoms / V;

  Histogram H_r, g_r, G_r, J_r;
  H_r.bins.resize(num_bins);
  g_r.bins.resize(num_bins);
  G_r.bins.resize(num_bins);
  J_r.bins.resize(num_bins);
  H_r.x_label = "r (Å)";
  H_r.title = "H(r) — Distance Histogram";
  H_r.y_label = "H(r)";
  H_r.x_unit = "Å";
  H_r.y_unit = "counts";
  H_r.description = "Distance Histogram";
  H_r.file_suffix = "_H";

  g_r.x_label = "r (Å)";
  g_r.title = "g(r) — Pair Distribution";
  g_r.y_label = "g(r)";
  g_r.x_unit = "Å";
  g_r.y_unit = "Å^-1";
  g_r.description = "Pair Distribution Function";
  g_r.file_suffix = "_g";

  G_r.x_label = "r (Å)";
  G_r.title = "G(r) — Reduced Pair Distribution";
  G_r.y_label = "G(r)";
  G_r.x_unit = "Å";
  G_r.y_unit = "Å^-1";
  G_r.description = "Reduced Pair Distribution Function";
  G_r.file_suffix = "_G_reduced";

  J_r.x_label = "r (Å)";
  J_r.title = "J(r) — Reduced Pair Distribution";
  J_r.y_label = "J(r)";
  J_r.x_unit = "Å";
  J_r.y_unit = "Å^-1";
  J_r.description = "Radial Distribution Function";
  J_r.file_suffix = "_J";

  for (size_t i = 0; i < num_bins; ++i) {
    const double r = (i + 0.5) * r_bin_width;
    H_r.bins[i] = r;
    g_r.bins[i] = r;
    G_r.bins[i] = r;
    J_r.bins[i] = r;
  }

  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      std::string key = getPartialKey(cell, i, j);
      auto &partial_hist = H_r.partials[key];
      partial_hist.assign(num_bins, 0.0);

      // First Pass: Accumulate the raw distance counts H(r) into histogram
      // bins. Distances are pre-computed in
      // `StructureAnalyzer::distance_tensor_`.
      for (const auto &dist : neighbors->distances()[i][j]) {
        if (dist < r_max) {
          size_t bin = static_cast<size_t>(dist / r_bin_width);
          if (bin < num_bins) {
            partial_hist[bin] += 1.0;
          }
        }
      }

      // For self-pairs (A-A), each pair is counted once in the upper triangular
      // DistanceTensor. We multiply by 2 to account for both A_1 -> A_2 and A_2
      // -> A_1 interactions.
      if (i == j) {
        correlation::math::scale_bins(partial_hist.data(), 2.0, num_bins);
      }
    }
  }

  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      std::string key = getPartialKey(cell, i, j);
      std::string inversekey = getInversePartialKey(cell, i, j);

      const std::string &sym_i = elements[i].symbol;
      const std::string &sym_j = elements[j].symbol;

      const double Ni = element_counts.at(sym_i);
      const double Nj = element_counts.at(sym_j);

      const auto &H_ij = H_r.partials.at(key);

      g_r.partials[key].assign(num_bins, 0.0);
      G_r.partials[key].assign(num_bins, 0.0);
      J_r.partials[key].assign(num_bins, 0.0);
      J_r.partials[inversekey].assign(num_bins, 0.0);

      // Second Pass: Normalize the raw counts H(r) into target distribution
      // functions. g(r) normalization constant: V / (4 * pi * dr * N_i * N_j).
      // The r^2 term is applied per-bin inside the SIMD kernel.
      const double g_norm_constant =
          V / (correlation::math::four_pi * dr * Ni * Nj);
      const double rho_j = Nj / V;
      const double inv_Ni_dr = 1.0 / (Ni * dr);
      const double inv_Nj_dr = 1.0 / (Nj * dr);
      const double pi4_rho_j = correlation::math::four_pi * rho_j;

      correlation::math::normalize_rdf_bins(
          H_ij.data(), g_r.bins.data(), g_norm_constant, inv_Ni_dr, inv_Nj_dr,
          pi4_rho_j, g_r.partials[key].data(), G_r.partials[key].data(),
          J_r.partials[key].data(), J_r.partials[inversekey].data(), num_bins);
    }
  }

  auto &total_g = g_r.partials["Total"];
  total_g.assign(num_bins, 0.0);
  for (const auto &[key, g_partial] : g_r.partials) {
    if (key == "Total")
      continue;

    double weight = ashcroft_weights.at(key);

    for (size_t k = 0; k < num_bins; ++k) {
      total_g[k] += g_partial[k] * weight;
    }
  }

  auto &total_J = J_r.partials["Total"];
  auto &total_G = G_r.partials["Total"];
  total_J.assign(num_bins, 0.0);
  total_G.assign(num_bins, 0.0);

  for (size_t k = 0; k < num_bins; ++k) {
    const double r = g_r.bins[k];
    if (r < 1e-9)
      continue;

    total_J[k] = correlation::math::four_pi * r * r * rho_0 * total_g[k];
    total_G[k] = correlation::math::four_pi * rho_0 * r * (total_g[k] - 1.0);
  }

  std::map<std::string, Histogram> results;
  results["J_r"] = std::move(J_r);
  results["g_r"] = std::move(g_r);
  results["G_r"] = std::move(G_r);

  return results;
}
