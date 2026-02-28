// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/DADCalculator.hpp"
#include "PhysicalData.hpp"
#include <stdexcept>

Histogram DADCalculator::calculate(const Cell &cell,
                                   const StructureAnalyzer *neighbors,
                                   double bin_width) {
  if (bin_width <= 0) {
    throw std::invalid_argument("Bin width must be positive");
  }

  // Dihedral angles are from -180 to 180 degrees.
  const double theta_min = -180.0;
  const double theta_max = 180.0;
  const double theta_range = theta_max - theta_min;

  const auto &elements = cell.elements();
  const size_t num_elements = elements.size();
  if (num_elements == 0)
    return {};

  const size_t num_bins = static_cast<size_t>((theta_range / bin_width) + 1);

  Histogram f_dihedral;
  f_dihedral.bin_label = "Dihedral Angle (°)";
  f_dihedral.bins.resize(num_bins);
  for (size_t i = 0; i < num_bins; ++i) {
    f_dihedral.bins[i] = theta_min + (i + 0.5) * bin_width;
  }

  for (size_t a = 0; a < num_elements; ++a) {
    for (size_t b = 0; b < num_elements; ++b) {
      for (size_t c = 0; c < num_elements; ++c) {
        for (size_t d = 0; d < num_elements; ++d) {

          const auto &angles_rad = neighbors->dihedrals()[a][b][c][d];
          if (angles_rad.empty())
            continue;

          std::string key = elements[a].symbol + "-" + elements[b].symbol +
                            "-" + elements[c].symbol + "-" + elements[d].symbol;

          auto &partial_hist = f_dihedral.partials[key];
          if (partial_hist.empty()) {
            partial_hist.assign(num_bins, 0.0);
          }

          for (const auto &angle_rad : angles_rad) {
            double angle_deg = angle_rad * constants::rad2deg;

            // clamp angle into [-180, 180]
            while (angle_deg <= -180.0)
              angle_deg += 360.0;
            while (angle_deg > 180.0)
              angle_deg -= 360.0;

            if (angle_deg >= theta_min && angle_deg <= theta_max) {
              size_t bin =
                  static_cast<size_t>((angle_deg - theta_min) / bin_width);

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
  }

  auto &total_f = f_dihedral.partials["Total"];
  total_f.assign(num_bins, 0.0);
  double total_counts = 0;

  for (const auto &[key, partial] : f_dihedral.partials) {
    if (key != "Total") {
      for (size_t i = 0; i < num_bins; ++i) {
        total_f[i] += partial[i];
        total_counts += partial[i];
      }
    }
  }

  if (total_counts < 1) {
    return f_dihedral;
  }

  const double normalization_factor = 1.0 / (total_counts * bin_width);
  for (auto &[key, partial] : f_dihedral.partials) {
    for (size_t i = 0; i < num_bins; ++i) {
      partial[i] *= normalization_factor;
    }
  }
  return f_dihedral;
}
