// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/XRDCalculator.hpp"
#include "PhysicalData.hpp"
#include <cmath>
#include <stdexcept>

Histogram
XRDCalculator::calculate(const Histogram &g_r_hist, const Cell &cell,
                         const std::map<std::string, double> &ashcroft_weights,
                         double lambda, double theta_min, double theta_max,
                         double bin_width) {
  if (bin_width <= 0) {
    throw std::invalid_argument("Bin width must be positive.");
  }

  const auto &r_bins = g_r_hist.bins;
  const double dr = r_bins[1] - r_bins[0];
  const double total_rho = cell.atomCount() / cell.volume();
  const double max_r = r_bins.back();

  size_t num_bins =
      static_cast<size_t>((theta_max - theta_min) / bin_width) + 1;
  Histogram xrd_hist;
  xrd_hist.bin_label = "2Theta (°)";
  xrd_hist.bins.resize(num_bins);
  std::vector<double> intensities(num_bins, 0.0);

  auto get_f_Q = [](const std::string &symbol, double Q) -> double {
    const auto &coeffs = AtomicFormFactors::get(symbol);
    double s = Q / (4.0 * constants::pi);
    double s2 = s * s;
    double f = coeffs[8];
    for (size_t i = 0; i < 4; ++i) {
      f += coeffs[2 * i] * std::exp(-coeffs[2 * i + 1] * s2);
    }
    return f;
  };

  std::map<std::string, std::vector<double>> partial_integrands;
  for (const auto &[key, g_partial] : g_r_hist.partials) {
    if (key == "Total")
      continue;
    partial_integrands[key].resize(g_partial.size());
    for (size_t k = 0; k < g_partial.size(); ++k) {
      partial_integrands[key][k] = r_bins[k] * (g_partial[k] - 1.0) * dr;
    }
  }

  std::map<std::string, double> concentrations;
  for (const auto &elem : cell.elements()) {
    double count = 0;
    for (const auto &atom : cell.atoms()) {
      if (atom.element().symbol == elem.symbol)
        count++;
    }
    concentrations[elem.symbol] = count / cell.atomCount();
  }

  for (size_t i = 0; i < num_bins; ++i) {
    double two_theta = theta_min + i * bin_width;
    xrd_hist.bins[i] = two_theta;

    double theta_rad = (two_theta / 2.0) * constants::deg2rad;
    double Q = 4.0 * constants::pi * std::sin(theta_rad) / lambda;

    if (Q < 1e-6) {
      intensities[i] = 0.0;
      continue;
    }

    double I_Q = 0.0;

    for (const auto &[sym, c] : concentrations) {
      double f = get_f_Q(sym, Q);
      I_Q += c * f * f;
    }

    for (const auto &[key, integrands] : partial_integrands) {
      double weight = ashcroft_weights.at(key);

      size_t dash_pos = key.find('-');
      std::string sym1 = key.substr(0, dash_pos);
      std::string sym2 = key.substr(dash_pos + 1);

      double f1 = get_f_Q(sym1, Q);
      double f2 = get_f_Q(sym2, Q);

      double integral = 0.0;
      for (size_t k = 0; k < integrands.size(); ++k) {
        integral += integrands[k] * std::sin(Q * r_bins[k]);
      }

      I_Q +=
          weight * f1 * f2 * (4.0 * constants::pi * total_rho / Q) * integral;
    }

    intensities[i] = I_Q;
  }

  xrd_hist.partials["Total"] = std::move(intensities);

  return xrd_hist;
}
