// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/RDFCalculator.hpp"
#include "PhysicalData.hpp"
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
} // namespace

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
  H_r.bin_label = "r (Å)";
  g_r.bin_label = "r (Å)";
  G_r.bin_label = "r (Å)";
  J_r.bin_label = "r (Å)";

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

      for (const auto &dist : neighbors->distances()[i][j]) {
        if (dist < r_max) {
          size_t bin = static_cast<size_t>(dist / r_bin_width);
          if (bin < num_bins) {
            partial_hist[bin] += 1.0;
          }
        }
      }
      if (i == j) {
        for (size_t k = 0; k < num_bins; ++k) {
          partial_hist[k] *= 2.0;
        }
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

      const double g_norm_constant = (V) / (4.0 * constants::pi * dr * Ni * Nj);
      const double rho_j = Nj / V;

      for (size_t k = 0; k < num_bins; ++k) {
        const double r = g_r.bins[k];
        if (r < 1e-9)
          continue;

        const double H = H_ij[k];

        g_r.partials[key][k] = H * g_norm_constant / (r * r);
        J_r.partials[key][k] = H / (Ni * dr);
        J_r.partials[inversekey][k] = H / (Nj * dr);
        G_r.partials[key][k] =
            4.0 * constants::pi * rho_j * r * (g_r.partials[key][k] - 1.0);
      }
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

    total_J[k] = 4.0 * constants::pi * r * r * rho_0 * total_g[k];
    total_G[k] = 4.0 * constants::pi * rho_0 * r * (total_g[k] - 1.0);
  }

  std::map<std::string, Histogram> results;
  results["J(r)"] = std::move(J_r);
  results["g(r)"] = std::move(g_r);
  results["G(r)"] = std::move(G_r);

  return results;
}
