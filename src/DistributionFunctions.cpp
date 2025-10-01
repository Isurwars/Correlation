// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/DistributionFunctions.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "../include/PhysicalData.hpp"
#include "../include/Smoothing.hpp"

//---------------------------------------------------------------------------//
//------------------------------- Constructor -------------------------------//
//---------------------------------------------------------------------------//

DistributionFunctions::DistributionFunctions(const Cell &cell, double cutoff,
                                             double bond_factor)
    : cell_(cell), neighbors_(nullptr), current_cutoff_(0.0),
      bond_factor_(bond_factor) {
  if (cutoff > 0.0) {
    ensureNeighborsComputed(cutoff);
  }
  calculateAshcroftWeights();
}

// Move Constructor
DistributionFunctions::DistributionFunctions(
    DistributionFunctions &&other) noexcept
    : cell_(other.cell_), neighbors_(std::move(other.neighbors_)),
      current_cutoff_(other.current_cutoff_), bond_factor_(other.bond_factor_),
      histograms_(std::move(other.histograms_)),
      ashcroft_weights_(std::move(other.ashcroft_weights_)) {

  other.current_cutoff_ = -1.0;
}

// Move Assignment Operator
DistributionFunctions &
DistributionFunctions::operator=(DistributionFunctions &&other) noexcept {
  if (this != &other) {
    neighbors_ = std::move(other.neighbors_);
    current_cutoff_ = other.current_cutoff_;
    bond_factor_ = other.bond_factor_;
    histograms_ = std::move(other.histograms_);
    ashcroft_weights_ = std::move(other.ashcroft_weights_);

    other.current_cutoff_ = -1.0;
  }
  return *this;
}

//--------------------------------------------------------------------------//
//---------------------------- Helper Functions ----------------------------//
//--------------------------------------------------------------------------//

const Histogram &
DistributionFunctions::getHistogram(const std::string &name) const {
  auto it = histograms_.find(name);
  if (it == histograms_.end()) {
    throw std::out_of_range("Histogram '" + name + "' not found.");
  }
  return it->second;
}

void DistributionFunctions::ensureNeighborsComputed(double r_max) {
  if (!neighbors_ || r_max > current_cutoff_) {
    neighbors_ =
        std::make_unique<StructureAnalyzer>(cell_, r_max, bond_factor_);
    current_cutoff_ = r_max;
  }
}

std::vector<std::string> DistributionFunctions::getAvailableHistograms() const {
  std::vector<std::string> keys;
  keys.reserve(histograms_.size());
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

void DistributionFunctions::calculateAshcroftWeights() {
  const auto &atoms = cell_.atoms();
  if (atoms.empty()) {
    throw std::invalid_argument("Cell must contain atoms");
  }
  const size_t num_atoms = atoms.size();

  // 1. Count the number of atoms for each element symbol.
  std::map<std::string, int> element_counts;
  for (const auto &atom : atoms) {
    element_counts[atom.element().symbol]++;
  }

  // 2. Get the list of unique elements present.
  const auto &elements = cell_.elements();

  // 3. Calculate w_ij = (N_i * N_j) / N_total^2 for all pairs.
  // This is the concentration-based weighting factor.
  for (size_t i = 0; i < elements.size(); ++i) {
    for (size_t j = i; j < elements.size(); ++j) {
      const auto &element_i = elements[i];
      const auto &element_j = elements[j];

      const double count_i =
          static_cast<double>(element_counts.at(element_i.symbol));
      const double count_j =
          static_cast<double>(element_counts.at(element_j.symbol));

      double weight = (count_i * count_j) / (num_atoms * num_atoms);
      if (i != j) {
        weight *= 2.0;
      }

      // Get the canonical key (e.g., "Si-O", not "O-Si")
      std::string key = getPartialKey(element_i.id.value, element_j.id.value);
      ashcroft_weights_[key] = weight;
    }
  }
}

//---------------------------------------------------------------------------//
//---------------------------- Smoothing Methods ----------------------------//
//---------------------------------------------------------------------------//
void DistributionFunctions::smooth(const std::string &name, double sigma,
                                   KernelType kernel) {
  if (histograms_.find(name) == histograms_.end()) {
    throw std::runtime_error("Histogram '" + name +
                             "' not found for smoothing.");
  }

  auto &hist = histograms_.at(name);
  hist.smoothed_partials.clear();

  for (const auto &[key, partial_values] : hist.partials) {
    hist.smoothed_partials[key] =
        KernelSmoothing(hist.bins, partial_values, sigma, kernel);
  }
}

void DistributionFunctions::smoothAll(double sigma, KernelType kernel) {
  for (const auto &[name, histogram] : histograms_) {
    smooth(name, sigma, kernel);
  }
}

//---------------------------------------------------------------------------//
//----------------------------- Calculation CN ------------------------------//
//---------------------------------------------------------------------------//

void DistributionFunctions::calculateCoordinationNumber() {
  if (!neighbors_) {
    throw std::logic_error(
        "Cannot calculate Coordination Number. Neighbor list has not been "
        "computed.");
  }

  const auto &atoms = cell_.atoms();
  const auto &all_neighbors = neighbors_->neighbors();

  // Stores the distribution of neighbors for each partial pair.
  // Key: "Central-Neighbor" (e.g. "Si-O"), Value: histogram vector.
  // The index of the vector is the CN, the value is the number of central atoms
  // with that CN. e.g., partial_dists["Si-O"][4] = count of Si atoms with 4 O
  // neighbors.
  std::map<std::string, std::vector<int>> partial_dists;
  size_t max_cn = 0; // Track the maximum CN found to size the final histogram.

  // 1. Iterate through each central atom to find its specific neighbor counts.
  for (size_t i = 0; i < atoms.size(); ++i) {
    const auto &central_atom = atoms[i];
    const std::string &central_symbol = central_atom.element().symbol;

    // Count neighbors of this specific atom, grouped by their element type.
    std::map<std::string, int> neighbor_counts_for_this_atom;
    for (const auto &neighbor : all_neighbors[i]) {
      const auto &neighbor_atom = atoms[neighbor.index.value];
      neighbor_counts_for_this_atom[neighbor_atom.element().symbol]++;
    }

    // 2. Tally these counts in the main distribution maps.
    for (const auto &[neighbor_symbol, count] : neighbor_counts_for_this_atom) {
      std::string key = central_symbol + "-" + neighbor_symbol;
      if (static_cast<size_t>(count) >= partial_dists[key].size()) {
        // Automatically resize the vector if we find a larger CN.
        partial_dists[key].resize(count + 1, 0);
      }
      partial_dists[key][count]++;
      if (static_cast<size_t>(count) > max_cn) {
        max_cn = count;
      }
    }
  }

  // 3. Assemble the final Histogram object.
  Histogram cn_histogram;
  const size_t num_bins = max_cn + 3;
  cn_histogram.bin_label = "# neighbors";

  // The bins are the integer coordination numbers [0, 1, 2, ...].
  cn_histogram.bins.resize(num_bins);
  std::iota(cn_histogram.bins.begin(), cn_histogram.bins.end(), 0.0);

  // 4. Populate the final partials, padding with zeros for a consistent size.
  for (auto &[key, dist_vector] : partial_dists) {
    dist_vector.resize(num_bins, 0); // Pad with zeros up to max_cn.
    cn_histogram.partials[key].assign(dist_vector.begin(), dist_vector.end());
  }

  // 5. Store the completed histogram in the main map.
  histograms_["CN"] = std::move(cn_histogram);
}

//---------------------------------------------------------------------------//
//----------------------------- Calculation RDF -----------------------------//
//---------------------------------------------------------------------------//

void DistributionFunctions::calculateRDF(double r_max, double r_bin_width) {
  if (r_bin_width <= 0) {
    throw std::invalid_argument("Bin width must be positive, got: " +
                                std::to_string(r_bin_width));
  }
  if (r_max <= 0) {
    throw std::invalid_argument("Cutoff radius must be positive, got: " +
                                std::to_string(r_max));
  }

  const double volume = cell_.volume();
  if (volume <= std::numeric_limits<double>::epsilon()) {
    throw std::logic_error("Cell volume must be positive, got: " +
                           std::to_string(volume));
  }

  const auto &elements = cell_.elements();
  const size_t num_elements = elements.size();
  const double num_atoms = static_cast<double>(cell_.atomCount());
  if (num_atoms == 0.0)
    return;

  // 1. Calculate and store element counts for normalization
  std::map<std::string, double> element_counts;
  for (const auto &atom : cell_.atoms()) {
    element_counts[atom.element().symbol]++;
  }

  const size_t num_bins = static_cast<size_t>(std::floor(r_max / r_bin_width));
  const double V = cell_.volume();
  const double dr = r_bin_width;
  const double rho_0 = num_atoms / V; // Total number density

  // Initialize histograms
  Histogram H_r, g_r, G_r, J_r; // H_r stores the raw counts H_ij(r)
  H_r.bins.resize(num_bins);
  g_r.bins.resize(num_bins);
  G_r.bins.resize(num_bins);
  J_r.bins.resize(num_bins);
  H_r.bin_label = "r";
  g_r.bin_label = "r";
  G_r.bin_label = "r";
  J_r.bin_label = "r";

  for (size_t i = 0; i < num_bins; ++i) {
    const double r = (i + 0.5) * r_bin_width;
    H_r.bins[i] = r;
    g_r.bins[i] = r;
    G_r.bins[i] = r;
    J_r.bins[i] = r;
  }

  // 2. Bin the distances (Raw Counts H_ij(r))
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      std::string key = getPartialKey(i, j);
      auto &partial_hist =
          H_r.partials[key]; // H_r.partials[key] stores H_ij(r)
      partial_hist.assign(num_bins, 0.0);

      for (const auto &dist : neighbors_->distances()[i][j]) {
        if (dist < r_max) {
          size_t bin = static_cast<size_t>(dist / r_bin_width);
          if (bin < num_bins) {
            // H_ij(r) is the raw pair count
            partial_hist[bin] += 1.0;
          }
        }
      }
      // Since StructureAnalyzer stores distances once per unique pair (i-j),
      // we must double the counts here to get the total pair count H_ij(r).
      if (i != j) {
        for (size_t k = 0; k < num_bins; ++k) {
          partial_hist[k] *= 2.0;
        }
      }
    }
  }

  // 3. Calculate Partials: g_ij(r), J_ij(r), G_ij(r) from H_ij(r)
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      std::string key = getPartialKey(i, j);

      const std::string &sym_i = elements[i].symbol;
      const std::string &sym_j = elements[j].symbol;

      const double Ni = element_counts.at(sym_i);
      const double Nj = element_counts.at(sym_j);

      const auto &H_ij = H_r.partials.at(key); // Raw pair count H_ij(r)

      g_r.partials[key].assign(num_bins, 0.0);
      G_r.partials[key].assign(num_bins, 0.0);
      J_r.partials[key].assign(num_bins, 0.0);

      // Normalization constant for g_ij(r): V / (4*pi*dr*N*N)
      const double g_norm_constant =
          (V) / (4.0 * constants::pi * dr * num_atoms * num_atoms);
      const double rho_j = Nj / V;
      const double w_ij = ashcroft_weights_[key];
      for (size_t k = 0; k < num_bins; ++k) {
        const double r = g_r.bins[k];
        if (r < 1e-9)
          continue;

        const double H = H_ij[k]; // Raw count

        // 3a. Calculate g_ij(r) = (H_ij / (r^2)) * g_norm_constant
        g_r.partials[key][k] = H * g_norm_constant / (r * r);

        // 3b. Calculate J_ij(r) (The RDF) = H_ij / (Ni * dr)
        // This is the number of j atoms around an i atom, per unit distance.
        J_r.partials[key][k] = H / (Ni * dr);

        // 3c. Calculate G_ij(r) = 4*pi*rho_j*r*(g_ij(r)-1)
        G_r.partials[key][k] =
            4.0 * constants::pi * rho_j * r * (g_r.partials[key][k] - w_ij);
      }
    }
  }

  // 4. Calculate Total Functions

  // 4a. Total g(r) = Sum of weighted partials g_ij(r)
  auto &total_g = g_r.partials["Total"];
  total_g.assign(num_bins, 0.0);
  for (const auto &[key, g_partial] : g_r.partials) {
    if (key == "Total")
      continue;

    for (size_t k = 0; k < num_bins; ++k) {
      total_g[k] += g_partial[k];
    }
  }

  // 4b. Total J(r) and G(r) from Total g(r)
  auto &total_J = J_r.partials["Total"];
  auto &total_G = G_r.partials["Total"];
  total_J.assign(num_bins, 0.0);
  total_G.assign(num_bins, 0.0);

  for (size_t k = 0; k < num_bins; ++k) {
    const double r = g_r.bins[k];
    if (r < 1e-9)
      continue;

    // J(r) = 4*pi*r^2*rho_0*g(r)
    total_J[k] = 4.0 * constants::pi * r * r * rho_0 * total_g[k];

    // G(r) = 4*pi*rho_0*r*(g(r)-1)
    total_G[k] = 4.0 * constants::pi * rho_0 * r * (total_g[k] - 1.0);
  }

  // 5. Store the final histograms
  // Note: H(r) is also stored, representing the raw pair counts.
  // histograms_["H(r)"] = std::move(H_r);
  histograms_["J(r)"] = std::move(J_r);
  histograms_["g(r)"] = std::move(g_r);
  histograms_["G(r)"] = std::move(G_r);
}

//---------------------------------------------------------------------------//
//----------------------------- Calculation PAD -----------------------------//
//---------------------------------------------------------------------------//

void DistributionFunctions::calculatePAD(double theta_cut, double bin_width) {
  if (bin_width <= 0) {
    throw std::invalid_argument("Bin width must be positive");
  }
  if (theta_cut <= 0) {
    throw std::invalid_argument("Theta cutoff must be positive");
  }
  if (theta_cut > 180.0) {
    throw std::invalid_argument("Theta cutoff cannot exceed 180 degrees");
  }

  const auto &elements = cell_.elements();
  const size_t num_elements = elements.size();
  if (num_elements == 0)
    return;

  const size_t num_bins = static_cast<size_t>(theta_cut / bin_width);

  Histogram f_theta;
  f_theta.bin_label = "theta";
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

        for (const auto &angle_rad : neighbors_->angles()[j][i][k]) {
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
      for (size_t i = 0; i < num_bins; ++i) {
        total_f[i] += partial[i];
        total_counts += partial[i];
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

//---------------------------------------------------------------------------//
//---------------------------- Calculation S(Q) -----------------------------//
//---------------------------------------------------------------------------//

void DistributionFunctions::calculateSQ(double q_max, double q_bin_width,
                                        double r_integration_max) {
  // 1. Input Validation & Dependency Check
  if (q_bin_width <= 0 || q_max <= 0) {
    throw std::invalid_argument("Q-space parameters must be positive.");
  }
  if (r_integration_max <= 0) {
    throw std::invalid_argument("Integration cutoff must be positive.");
  }
  if (histograms_.find("g(r)") == histograms_.end()) {
    throw std::logic_error(
        "Cannot calculate S(Q). Please calculate g(r) first by calling "
        "calculateRDF().");
  }

  const auto &g_r_hist = histograms_.at("g(r)");
  const auto &r_bins = g_r_hist.bins;
  if (r_bins.size() < 2) {
    throw std::logic_error("Insufficient r-bins for integration.");
  }

  // Determine r_max: Use the minimum of user's requested limit and the max
  // available in g(r).
  const double max_r_in_hist = r_bins.back();
  // Correct usage of min for safe integration range.
  const double r_max = std::min(r_integration_max, max_r_in_hist);

  const double dr = r_bins[1] - r_bins[0];
  const double total_rho = cell_.atomCount() / cell_.volume();

  // 2. Setup S(Q) Histogram
  Histogram s_q_hist;
  s_q_hist.bin_label = "Q";
  const double q_min = 1.0;
  const size_t num_q_bins =
      static_cast<size_t>(std::floor((q_max - q_min) / q_bin_width));

  if (num_q_bins == 0) {
    throw std::invalid_argument(
        "Q_max is too small for the given Q_bin_width.");
  }

  s_q_hist.bins.resize(num_q_bins);
  for (size_t i = 0; i < num_q_bins; ++i) {
    s_q_hist.bins[i] = q_min + (i + 0.5) * q_bin_width;
  }

  // 3. Find integration limit in r-bins
  size_t j_max = r_bins.size();
  for (size_t j = 0; j < r_bins.size(); ++j) {
    if (r_bins[j] > r_max) {
      j_max = j;
      break;
    }
  }

  // Initialize the vector for the total Number Structure Factor S_num(Q)
  std::vector<double> total_s_q_num(num_q_bins, 0.0);

  // 4. Iterate over all partials to calculate S_ij(Q)
  for (const auto &[key, g_r_partial] : g_r_hist.partials) {
    if (key == "Total") {
      continue;
    }
    // Ashcroft Weights
    double weight = ashcroft_weights_[key];
    // Faber-Ziman normalization
    double composition_sqroot_factor = std::sqrt(weight);

    auto &s_q_partial = s_q_hist.partials[key];
    s_q_partial.assign(num_q_bins, 0.0);

    // Precompute windowed integrand terms: r * (g(r) - 1) * Window(r) * dr
    std::vector<double> integrand_term(j_max);
    for (size_t j = 0; j < j_max; ++j) {
      const double r = r_bins[j];
      double window = 1.0;
      if (r > r_max * 0.8) {
        const double x = (r - 0.8 * r_max) / (0.2 * r_max);
        window = std::sin(constants::pi * x / 2.0);
      }

      integrand_term[j] = r * ((g_r_partial[j] / weight) - 1.0) * window * dr;
    }

    // 5. Calculate S_ij(Q) for each Q bin
    for (size_t i = 0; i < num_q_bins; ++i) {
      const double Q = s_q_hist.bins[i];
      double integral = 0.0;

      // Numerical integration for ∫ [r * (g(r)-1) * W(r) * sin(Qr)] dr
      for (size_t j = 0; j < j_max; ++j) {
        const double r = r_bins[j];
        // Integral: r * (g-1) * W * sin(Qr)
        integral += integrand_term[j] * std::sin(Q * r);
      }

      // Delta function: 1.0 for S_ii(Q), 0.0 for S_ij(Q) (i != j)
      double delta_ij = (key.find('-') == key.rfind('-')) ? 1.0 : 0.0;

      // S_ij(Q) = delta_ij + (4*pi*rho_total*sqrt(c_i*c_j)/Q) * Integral[...]
      s_q_partial[i] = delta_ij + (4.0 * constants::pi * total_rho *
                                   composition_sqroot_factor / Q) *
                                      integral;

      // 6a. Accumulate for the Total Number Structure Factor S_num(Q)
      total_s_q_num[i] += s_q_partial[i] * weight;
    }
  }

  // 6b. Finalize Total S(Q)
  s_q_hist.partials["Total"] = std::move(total_s_q_num);

  // Store the result
  histograms_["S(Q)"] = std::move(s_q_hist);
}

void DistributionFunctions::calculateXRD(double lambda, double theta_min,
                                         double theta_max, double bin_width) {
  // Modernized implementation would go here...
}
