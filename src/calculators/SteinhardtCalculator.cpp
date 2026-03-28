// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "math/LinearAlgebra.hpp"
#include "calculators/SteinhardtCalculator.hpp"
#include "math/SpecialFunctions.hpp"
#include "calculators/CalculatorFactory.hpp"
#include <cmath>
#include <stdexcept>

namespace {
bool registered = CalculatorFactory::instance().registerCalculator(
    std::make_unique<SteinhardtCalculator>());
} // namespace

std::complex<double> SteinhardtCalculator::sphericalHarmonic(int l, int m,
                                                             double theta,
                                                             double phi) {
  if (m >= 0) {
    double P_lm = correlation::math::special::sph_legendre(l, m, theta);
    return P_lm * std::polar(1.0, m * phi);
  } else {
    // For negative m: Y_l^{-m} = (-1)^m (Y_l^m)*
    int abs_m = -m;
    double P_lm = correlation::math::special::sph_legendre(l, abs_m, theta);
    std::complex<double> Y_l_m = P_lm * std::polar(1.0, abs_m * phi);
    std::complex<double> Y_l_minus_m = std::conj(Y_l_m);
    if (abs_m % 2 != 0) {
      Y_l_minus_m = -Y_l_minus_m;
    }
    return Y_l_minus_m;
  }
}

double SteinhardtCalculator::wigner3j(int j1, int j2, int j3, int m1, int m2,
                                      int m3) {
  if (m1 + m2 + m3 != 0)
    return 0.0;
  if (j3 < std::abs(j1 - j2) || j3 > j1 + j2)
    return 0.0;
  if (std::abs(m1) > j1 || std::abs(m2) > j2 || std::abs(m3) > j3)
    return 0.0;

  double delta = correlation::math::special::factorial(j1 + j2 - j3) *
                 correlation::math::special::factorial(j1 - j2 + j3) *
                 correlation::math::special::factorial(-j1 + j2 + j3) /
                 correlation::math::special::factorial(j1 + j2 + j3 + 1);
  delta = std::sqrt(delta);

  double comp =
      correlation::math::special::factorial(j1 - m1) * correlation::math::special::factorial(j1 + m1) *
      correlation::math::special::factorial(j2 - m2) * correlation::math::special::factorial(j2 + m2) *
      correlation::math::special::factorial(j3 - m3) * correlation::math::special::factorial(j3 + m3);
  comp = std::sqrt(comp);

  double phase1 = ((j1 - j2 - m3) % 2 != 0) ? -1.0 : 1.0;

  int k_min = std::max(0, std::max(j2 - j3 - m1, j1 - j3 + m2));
  int k_max = std::min(j1 + j2 - j3, std::min(j1 - m1, j2 + m2));

  double sum = 0.0;
  for (int k = k_min; k <= k_max; ++k) {
    double k_phase = (k % 2 != 0) ? -1.0 : 1.0;
    double denom = correlation::math::special::factorial(k) *
                   correlation::math::special::factorial(j1 + j2 - j3 - k) *
                   correlation::math::special::factorial(j1 - m1 - k) *
                   correlation::math::special::factorial(j2 + m2 - k) *
                   correlation::math::special::factorial(j3 - j2 + m1 + k) *
                   correlation::math::special::factorial(j3 - j1 - m2 + k);
    sum += k_phase / denom;
  }

  return phase1 * delta * comp * sum;
}

void SteinhardtCalculator::calculateFrame(
    DistributionFunctions &df, const AnalysisSettings &settings) const {
  auto histograms = calculate(df.cell(), df.neighbors());
  for (auto &[name, hist] : histograms) {
    df.addHistogram(name, std::move(hist));
  }
}

std::map<std::string, Histogram>
SteinhardtCalculator::calculate(const Cell &cell,
                                const StructureAnalyzer *neighbors) {
  if (!neighbors) {
    throw std::logic_error(
        "Cannot calculate Steinhardt Parameters. Neighbor list has not been "
        "computed.");
  }

  const auto &atoms = cell.atoms();
  const auto &neighbor_graph = neighbors->neighborGraph();
  size_t num_atoms = atoms.size();

  // We want to compute Q4, Q6, W6 histograms mapping Q values to probability.
  // First we calculate q_lm(i) for l=4, l=6.
  std::vector<std::vector<std::complex<double>>> q4m(
      num_atoms, std::vector<std::complex<double>>(9, 0.0));
  std::vector<std::vector<std::complex<double>>> q6m(
      num_atoms, std::vector<std::complex<double>>(13, 0.0));

  std::vector<double> Q4(num_atoms, 0.0);
  std::vector<double> Q6(num_atoms, 0.0);
  std::vector<double> W6(num_atoms, 0.0);
  std::vector<double> W6_hat(num_atoms, 0.0);

  double global_Q4_factor = std::sqrt(4.0 * M_PI / 9.0);
  double global_Q6_factor = std::sqrt(4.0 * M_PI / 13.0);

  for (size_t i = 0; i < num_atoms; ++i) {
    const auto &atom_neighbors = neighbor_graph.getNeighbors(i);
    size_t nb = atom_neighbors.size();
    if (nb < 2) {
      continue;
    }

    for (const auto &neighbor : atom_neighbors) {
      // r_ij vector
      correlation::math::linalg::Vector3<double> r_ij = neighbor.r_ij;
      double r = neighbor.distance;
      if (r == 0)
        continue; // Safety check

      double theta = std::acos(r_ij.z() / r);
      double phi = std::atan2(r_ij.y(), r_ij.x());

      for (int m = -4; m <= 4; ++m) {
        q4m[i][m + 4] += sphericalHarmonic(4, m, theta, phi);
      }
      for (int m = -6; m <= 6; ++m) {
        q6m[i][m + 6] += sphericalHarmonic(6, m, theta, phi);
      }
    }

    for (int m = -4; m <= 4; ++m) {
      q4m[i][m + 4] /= static_cast<double>(nb);
    }
    for (int m = -6; m <= 6; ++m) {
      q6m[i][m + 6] /= static_cast<double>(nb);
    }

    double sum_sq_4 = 0.0;
    for (int m = -4; m <= 4; ++m) {
      sum_sq_4 += std::norm(q4m[i][m + 4]);
    }
    Q4[i] = global_Q4_factor * std::sqrt(sum_sq_4);

    double sum_sq_6 = 0.0;
    for (int m = -6; m <= 6; ++m) {
      sum_sq_6 += std::norm(q6m[i][m + 6]);
    }
    Q6[i] = global_Q6_factor * std::sqrt(sum_sq_6);

    // Compute W6
    double w6 = 0.0;
    for (int m1 = -6; m1 <= 6; ++m1) {
      for (int m2 = -6; m2 <= 6; ++m2) {
        int m3 = -(m1 + m2);
        if (m3 < -6 || m3 > 6)
          continue;
        double w3j = wigner3j(6, 6, 6, m1, m2, m3);
        if (w3j == 0.0)
          continue;

        std::complex<double> prod =
            q6m[i][m1 + 6] * q6m[i][m2 + 6] * q6m[i][m3 + 6];
        w6 += w3j * prod.real(); // product should be real theoretically
      }
    }
    W6[i] = w6;
    if (sum_sq_6 > 1e-12) {
      W6_hat[i] = w6 / std::pow(sum_sq_6, 1.5);
    }
  }

  // Distribution settings
  size_t bins_Q = 100;
  double Q_max = 1.0;
  double dQ = Q_max / bins_Q;

  size_t bins_W = 100;
  double W_min = -0.2;
  double W_max = 0.2;
  double dW = (W_max - W_min) / bins_W;

  Histogram hist_Q4;
  hist_Q4.bin_label = "Q4";
  hist_Q4.bins.resize(bins_Q);
  for (size_t b = 0; b < bins_Q; ++b)
    hist_Q4.bins[b] = (b + 0.5) * dQ;

  Histogram hist_Q6;
  hist_Q6.bin_label = "Q6";
  hist_Q6.bins.resize(bins_Q);
  for (size_t b = 0; b < bins_Q; ++b)
    hist_Q6.bins[b] = (b + 0.5) * dQ;

  Histogram hist_W6;
  hist_W6.bin_label = "W6_hat";
  hist_W6.bins.resize(bins_W);
  for (size_t b = 0; b < bins_W; ++b)
    hist_W6.bins[b] = W_min + (b + 0.5) * dW;

  // Track partials by element symbol and total
  std::map<std::string, std::vector<double>> partials_Q4;
  std::map<std::string, std::vector<double>> partials_Q6;
  std::map<std::string, std::vector<double>> partials_W6;

  for (const auto &elem : cell.elements()) {
    partials_Q4[elem.symbol].assign(bins_Q, 0.0);
    partials_Q6[elem.symbol].assign(bins_Q, 0.0);
    partials_W6[elem.symbol].assign(bins_W, 0.0);
  }
  partials_Q4["Total"].assign(bins_Q, 0.0);
  partials_Q6["Total"].assign(bins_Q, 0.0);
  partials_W6["Total"].assign(bins_W, 0.0);

  double num_atoms_f = 0.0;
  for (size_t i = 0; i < num_atoms; ++i) {
    if (neighbor_graph.getNeighbors(i).size() < 2)
      continue;

    num_atoms_f += 1.0;
    const std::string &symbol = atoms[i].element().symbol;

    if (Q4[i] >= 0 && Q4[i] < Q_max) {
      size_t b = static_cast<size_t>(Q4[i] / dQ);
      partials_Q4[symbol][b] += 1.0;
      partials_Q4["Total"][b] += 1.0;
    }

    if (Q6[i] >= 0 && Q6[i] < Q_max) {
      size_t b = static_cast<size_t>(Q6[i] / dQ);
      partials_Q6[symbol][b] += 1.0;
      partials_Q6["Total"][b] += 1.0;
    }

    if (W6_hat[i] >= W_min && W6_hat[i] < W_max) {
      size_t b = static_cast<size_t>((W6_hat[i] - W_min) / dW);
      partials_W6[symbol][b] += 1.0;
      partials_W6["Total"][b] += 1.0;
    }
  }

  if (num_atoms_f > 0) {
    for (auto &[key, vec] : partials_Q4) {
      for (auto &val : vec)
        val /= (num_atoms_f * dQ);
      hist_Q4.partials[key] = vec;
    }
    for (auto &[key, vec] : partials_Q6) {
      for (auto &val : vec)
        val /= (num_atoms_f * dQ);
      hist_Q6.partials[key] = vec;
    }
    for (auto &[key, vec] : partials_W6) {
      for (auto &val : vec)
        val /= (num_atoms_f * dW);
      hist_W6.partials[key] = vec;
    }
  }

  std::map<std::string, Histogram> hists;
  hists["Q4"] = std::move(hist_Q4);
  hists["Q6"] = std::move(hist_Q6);
  hists["W6_hat"] = std::move(hist_W6);

  return hists;
}
