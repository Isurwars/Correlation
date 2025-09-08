// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/DistributionFunctions.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "../include/Constants.hpp"
#include "../include/Smoothing.hpp"

//---------------------------------------------------------------------------//
//------------------------------- Constructor -------------------------------//
//---------------------------------------------------------------------------//

DistributionFunctions::DistributionFunctions(const Cell &cell,
                                             const NeighborList &neighbors)
    : cell_(cell), neighbors_(neighbors) {}

//--------------------------------------------------------------------------//
//---------------------------- Helper Functions ----------------------------//
//--------------------------------------------------------------------------//

std::vector<std::string> DistributionFunctions::getAvailableHistograms() const {
  std::vector<std::string> keys;
  // Iterate through the map of histograms and extract the key for each entry.
  for (const auto &[key, val] : histograms_) {
    keys.push_back(key);
  }
  return keys;
}

std::string DistributionFunctions::getPartialKey(int type1, int type2) const {
  const auto &elements = cell_.elements();
  // Ensure consistent ordering for pairs (e.g., Si-O is same as O-Si)
  if (type1 > type2)
    std::swap(type1, type2);
  return elements[type1].symbol + "-" + elements[type2].symbol;
}

//---------------------------------------------------------------------------//
//--------------------------- Calculation Methods ---------------------------//
//---------------------------------------------------------------------------//

void DistributionFunctions::calculateCoordinationNumber() {
  // Modernized implementation would go here...
}

void DistributionFunctions::calculateRDF(double r_cut, double bin_width,
                                         bool normalize) {
  if (bin_width <= 0)
    throw std::invalid_argument("Bin width must be positive.");
  if (r_cut <= 0)
    throw std::invalid_argument("Cutoff radius must be positive.");
  if (cell_.volume() <= 1e-9)
    throw std::logic_error("Cell volume must be positive.");

  const auto &elements = cell_.elements();
  const size_t num_elements = elements.size();
  const size_t num_atoms = cell_.atomCount();
  if (num_atoms == 0)
    return;

  const size_t num_bins = static_cast<size_t>(std::floor(r_cut / bin_width));

  // Initialize histograms
  Histogram J_r, g_r, G_r;
  J_r.bins.resize(num_bins);
  g_r.bins.resize(num_bins);
  G_r.bins.resize(num_bins);

  for (size_t i = 0; i < num_bins; ++i) {
    const double r = (i + 0.5) * bin_width;
    J_r.bins[i] = r;
    g_r.bins[i] = r;
    G_r.bins[i] = r;
  }

  // Bin the distances
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      std::string key = getPartialKey(i, j);
      auto &partial_hist = J_r.partials[key];
      partial_hist.assign(num_bins, 0.0);

      for (const auto &dist : neighbors_.distances()[i][j]) {
        if (dist < r_cut) {
          size_t bin = static_cast<size_t>(dist / bin_width);
          if (bin < num_bins) {
            partial_hist[bin] += (i == j) ? 1.0 : 2.0;
          }
        }
      }
    }
  }

  // Calculate total J(r)
  auto &total_J = J_r.partials["Total"];
  total_J.assign(num_bins, 0.0);
  for (const auto &[key, partial] : J_r.partials) {
    if (key != "Total") {
      for (size_t i = 0; i < num_bins; ++i) {
        total_J[i] += partial[i];
      }
    }
  }

  // --- Normalization and Calculation of g(r) and G(r) ---
  const double total_rho = num_atoms / cell_.volume();

  // TODO: Implement Ashcroft-Waseda weighting if normalize is false

  for (auto const &[key, J_partial] : J_r.partials) {
    g_r.partials[key].assign(num_bins, 0.0);
    G_r.partials[key].assign(num_bins, 0.0);

    for (size_t i = 0; i < num_bins; ++i) {
      const double r = J_r.bins[i];
      if (r < 1e-9)
        continue;

      double norm_factor =
          4.0 * constants::pi * r * r * total_rho * bin_width * num_atoms;
      if (norm_factor > 1e-9) {
        g_r.partials[key][i] = J_partial[i] / norm_factor;
      }

      G_r.partials[key][i] =
          4.0 * constants::pi * r * total_rho * (g_r.partials[key][i] - 1.0);
    }
  }

  histograms_["J(r)"] = std::move(J_r);
  histograms_["g(r)"] = std::move(g_r);
  histograms_["G(r)"] = std::move(G_r);
}

void DistributionFunctions::calculatePAD(double theta_cut, double bin_width) {
  if (bin_width <= 0)
    throw std::invalid_argument("Bin width must be positive.");
  if (theta_cut <= 0)
    throw std::invalid_argument("Theta cutoff must be positive.");

  const auto &elements = cell_.elements();
  const size_t num_elements = elements.size();
  if (num_elements == 0)
    return;

  const size_t num_bins = static_cast<size_t>(theta_cut / bin_width);

  Histogram f_theta;
  f_theta.bins.resize(num_bins);
  for (size_t i = 0; i < num_bins; ++i) {
    f_theta.bins[i] = (i + 0.5) * bin_width;
  }

  // Bin the angles
  for (size_t i = 0; i < num_elements; ++i) { // central atom
    for (size_t j = 0; j < num_elements; ++j) {
      for (size_t k = j; k < num_elements; ++k) {
        std::string key = elements[j].symbol + "-" + elements[i].symbol + "-" +
                          elements[k].symbol;
        auto &partial_hist = f_theta.partials[key];
        partial_hist.assign(num_bins, 0.0);

        for (const auto &angle_rad : neighbors_.angles()[j][i][k]) {
          double angle_deg = angle_rad * constants::rad2deg;
          if (angle_deg < theta_cut) {
            size_t bin = static_cast<size_t>(angle_deg / bin_width);
            if (bin < num_bins) {
              partial_hist[bin]++;
            }
          }
        }
      }
    }
  }

  // Calculate total f(theta) and normalization factor
  auto &total_f = f_theta.partials["Total"];
  total_f.assign(num_bins, 0.0);
  double total_counts = 0;
  for (const auto &[key, partial] : f_theta.partials) {
    if (key != "Total") {
      for (double count : partial) {
        total_counts += count;
      }
    }
  }

  if (total_counts < 1) { // No angles found
    histograms_["f(theta)"] = std::move(f_theta);
    return;
  }

  // Normalize
  const double normalization_factor = 1.0 / (total_counts * bin_width);
  for (auto &[key, partial] : f_theta.partials) {
    for (size_t i = 0; i < num_bins; ++i) {
      partial[i] *= normalization_factor;
      if (key != "Total") {
        total_f[i] += partial[i];
      }
    }
  }
  histograms_["f(theta)"] = std::move(f_theta);
}

void DistributionFunctions::smooth(const std::string &name, double sigma) {
  if (histograms_.find(name) == histograms_.end()) {
    throw std::runtime_error("Histogram '" + name +
                             "' not found for smoothing.");
  }
  auto &hist = histograms_.at(name);
  hist.smoothed_partials.clear();

  for (const auto &[key, partial_values] : hist.partials) {
    hist.smoothed_partials[key] = KernelSmoothing(
        hist.bins, partial_values, sigma, 1); // Assuming kernel type 1
  }
}

void DistributionFunctions::smoothAll(double sigma) {
  // This method provides a convenient way to smooth all histograms
  // by iterating through the map and calling the single-histogram
  // 'smooth' method on each one.
  for (const auto &[name, histogram] : histograms_) {
    // The 'name' is the key from the map (e.g., "g(r)")
    smooth(name, sigma);
  }
}

// Placeholder for other functions - they would be modernized similarly
void DistributionFunctions::calculateSQ(double q_max, double q_bin_width,
                                        double r_integration_max) {
  // Modernized implementation would go here...
}
void DistributionFunctions::calculateXRD(double lambda, double theta_min,
                                         double theta_max, double bin_width) {
  // Modernized implementation would go here...
}
