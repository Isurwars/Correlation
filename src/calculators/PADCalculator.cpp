// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/PADCalculator.hpp"
#include "PhysicalData.hpp"
#include <stdexcept>

Histogram PADCalculator::calculate(const Cell &cell,
                                   const StructureAnalyzer *neighbors,
                                   double bin_width) {
  if (bin_width <= 0) {
    throw std::invalid_argument("Bin width must be positive");
  }

  const double theta_cut = 180.0;

  const auto &elements = cell.elements();
  const size_t num_elements = elements.size();
  if (num_elements == 0)
    return {};

  const size_t num_bins = static_cast<size_t>((theta_cut / bin_width) + 1);

  Histogram f_theta;
  f_theta.bin_label = "theta (°)";
  f_theta.bins.resize(num_bins);
  for (size_t i = 0; i < num_bins; ++i) {
    f_theta.bins[i] = (i + 0.5) * bin_width;
  }

  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
      for (size_t k = j; k < num_elements; ++k) {
        std::string key = elements[j].symbol + "-" + elements[i].symbol + "-" +
                          elements[k].symbol;
        auto &partial_hist = f_theta.partials[key];
        partial_hist.assign(num_bins, 0.0);

        for (const auto &angle_rad : neighbors->angles()[j][i][k]) {
          double angle_deg = angle_rad * constants::rad2deg;

          if (angle_deg <= theta_cut + 1e-5) {
            size_t bin = static_cast<size_t>(angle_deg / bin_width);

            if (bin == num_bins && bin > 0) {
              bin = num_bins - 1;
            }

            if (bin < num_bins) {
              partial_hist[bin]++;
            }
          }
        }
      }
    }
  }

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

  if (total_counts < 1) {
    return f_theta;
  }

  const double normalization_factor = 1.0 / (total_counts * bin_width);
  for (auto &[key, partial] : f_theta.partials) {
    for (size_t i = 0; i < num_bins; ++i) {
      partial[i] *= normalization_factor;
    }
  }
  return f_theta;
}
