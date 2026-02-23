// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/SQCalculator.hpp"
#include "PhysicalData.hpp"
#include <cmath>
#include <stdexcept>

Histogram
SQCalculator::calculate(const Histogram &g_r_hist, const Cell &cell,
                        const std::map<std::string, double> &ashcroft_weights,
                        double q_max, double q_bin_width,
                        double r_integration_max) {
  if (q_bin_width <= 0 || q_max <= 0) {
    throw std::invalid_argument("Q-space parameters must be positive.");
  }
  if (r_integration_max <= 0) {
    throw std::invalid_argument("Integration cutoff must be positive.");
  }

  const auto &r_bins = g_r_hist.bins;
  if (r_bins.size() < 2) {
    throw std::logic_error("Insufficient r-bins for integration.");
  }

  const double max_r_in_hist = r_bins.back();
  const double r_max_val = std::min(r_integration_max, max_r_in_hist);

  const double dr = r_bins[1] - r_bins[0];
  const double total_rho = cell.atomCount() / cell.volume();

  Histogram s_q_hist;
  s_q_hist.bin_label = "Q (Å⁻¹)";
  const size_t num_q_bins =
      static_cast<size_t>(std::floor(q_max / q_bin_width));

  if (num_q_bins == 0) {
    throw std::invalid_argument(
        "Q_max is too small for the given Q_bin_width.");
  }

  s_q_hist.bins.resize(num_q_bins);
  for (size_t i = 0; i < num_q_bins; ++i) {
    s_q_hist.bins[i] = (i + 0.5) * q_bin_width;
  }

  size_t j_max = r_bins.size();
  for (size_t j = 0; j < r_bins.size(); ++j) {
    if (r_bins[j] > r_max_val) {
      j_max = j;
      break;
    }
  }

  std::vector<double> total_s_q_num(num_q_bins, 0.0);

  for (const auto &[key, g_r_partial] : g_r_hist.partials) {
    if (key == "Total") {
      continue;
    }
    double weight = ashcroft_weights.at(key);
    double composition_sqroot_factor;
    bool is_identical = false;

    size_t dash_pos = key.find('-');
    if (dash_pos != std::string::npos) {
      std::string sym1 = key.substr(0, dash_pos);
      std::string sym2 = key.substr(dash_pos + 1);
      is_identical = (sym1 == sym2);
    }

    if (is_identical) {
      composition_sqroot_factor = std::sqrt(weight);
    } else {
      composition_sqroot_factor = std::sqrt(weight / 2.0);
    }

    auto &s_q_partial = s_q_hist.partials[key];
    s_q_partial.assign(num_q_bins, 0.0);
    std::vector<double> integrand_term(j_max);
    for (size_t j = 0; j < j_max; ++j) {
      const double r = r_bins[j];
      double window = 1.0;
      if (r > r_max_val * 0.8) {
        const double x = (r - 0.8 * r_max_val) / (0.2 * r_max_val);
        window = std::cos(constants::pi * x / 2.0);
      }

      integrand_term[j] = r * (g_r_partial[j] - 1.0) * window * dr;
    }

    for (size_t i = 0; i < num_q_bins; ++i) {
      const double Q = s_q_hist.bins[i];
      double integral = 0.0;

      for (size_t j = 0; j < j_max; ++j) {
        const double r = r_bins[j];
        integral += integrand_term[j] * std::sin(Q * r);
      }

      double delta_ij = is_identical ? 1.0 : 0.0;

      s_q_partial[i] = delta_ij + (4.0 * constants::pi * total_rho *
                                   composition_sqroot_factor / Q) *
                                      integral;

      total_s_q_num[i] += s_q_partial[i] * weight;
    }
  }

  s_q_hist.partials["Total"] = std::move(total_s_q_num);

  return s_q_hist;
}
