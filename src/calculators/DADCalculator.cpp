/**
 * @file DADCalculator.cpp
 * @brief Implementation of the dihedral angle distribution calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/DADCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"

#include <stdexcept>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
// NOLINTNEXTLINE(cert-err58-cpp)
const bool registered = CalculatorFactory::registerTypeSafe<DADCalculator>("DADCalculator");

/**
 * @brief Helper to initialize the Dihedral-Angle Distribution histogram.
 */
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
correlation::analysis::Histogram initializeHistogram(size_t num_bins, real_t theta_min, real_t bin_width) {
  correlation::analysis::Histogram f_dihedral;
  f_dihedral.x_label = "φ";
  f_dihedral.title = "Dihedral-Angle Distribution";
  f_dihedral.y_label = "P(φ)";
  f_dihedral.x_unit = "°";
  f_dihedral.y_unit = "°⁻¹";
  f_dihedral.description = "Dihedral Angle Distribution";
  f_dihedral.file_suffix = "_DAD";
  f_dihedral.bins.resize(num_bins);
  for (size_t idx = 0; idx < num_bins; ++idx) {
    f_dihedral.bins[idx] = theta_min + (static_cast<real_t>(idx) + static_cast<real_t>(0.5)) * bin_width;
  }
  return f_dihedral;
}

/**
 * @brief Helper to process and bin a vector of dihedral angles.
 */
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void processDihedralAngles(const std::vector<real_t> &angles_rad, std::vector<real_t> &partial_hist, real_t bin_width,
                           size_t num_bins, real_t theta_min, real_t theta_max) {
  for (const auto &angle_rad : angles_rad) {
    real_t angle_deg = angle_rad * static_cast<real_t>(correlation::math::rad_to_deg);

    // clamp angle into [-180, 180]
    while (angle_deg <= static_cast<real_t>(-180.0)) {
      angle_deg += static_cast<real_t>(360.0);
    }
    while (angle_deg > static_cast<real_t>(180.0)) {
      angle_deg -= static_cast<real_t>(360.0);
    }

    if (angle_deg >= theta_min && angle_deg <= theta_max) {
      auto bin = static_cast<size_t>((angle_deg - theta_min) / bin_width);

      if (bin == num_bins && bin > 0) {
        bin = num_bins - 1;
      }

      if (bin < num_bins) {
        partial_hist[bin]++;
      }
    }
  }
}

/**
 * @brief Helper to sum partial histograms and normalize them.
 */
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void normalizeAndScale(correlation::analysis::Histogram &f_dihedral, size_t num_bins, real_t bin_width) {
  auto &total_f = f_dihedral.partials["Total"];
  total_f.assign(num_bins, static_cast<real_t>(0.0));
  real_t total_counts = static_cast<real_t>(0.0);

  for (const auto &[key, partial] : f_dihedral.partials) {
    if (key != "Total") {
      for (size_t idx = 0; idx < num_bins; ++idx) {
        total_f[idx] += partial[idx];
        total_counts += partial[idx];
      }
    }
  }

  if (total_counts >= static_cast<real_t>(1.0)) {
    const real_t normalization_factor = static_cast<real_t>(1.0) / (total_counts * bin_width);
    for (auto &[key, partial] : f_dihedral.partials) {
      for (auto &val : partial) {
        val *= normalization_factor;
      }
    }
  }
}
} // namespace

void DADCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                   const correlation::analysis::AnalysisSettings &settings) const {
  dists.addHistogram("DAD", calculate(dists.cell(), dists.neighbors(), settings.dihedral_bin_width));
}

correlation::analysis::Histogram DADCalculator::calculate(const correlation::core::Cell &cell,
                                                          const correlation::analysis::StructureAnalyzer *neighbors,
                                                          real_t bin_width) {
  if (bin_width <= 0) {
    throw std::invalid_argument("Bin width must be positive");
  }
  if (neighbors == nullptr) {
    throw std::logic_error("Cannot calculate DAD. Neighbor list has not been computed.");
  }

  // Dihedral angles are from -180 to 180 degrees.
  const real_t theta_min = static_cast<real_t>(-180.0);
  const real_t theta_max = static_cast<real_t>(180.0);
  const real_t theta_range = theta_max - theta_min;

  const auto &elements = cell.elements();
  const size_t num_elements = elements.size();
  if (num_elements == 0) {
    return {};
  }

  const auto num_bins = static_cast<size_t>((theta_range / bin_width) + 1);
  correlation::analysis::Histogram f_dihedral = initializeHistogram(num_bins, theta_min, bin_width);

  for (size_t idx_a = 0; idx_a < num_elements; ++idx_a) {
    for (size_t idx_b = 0; idx_b < num_elements; ++idx_b) {
      for (size_t idx_c = 0; idx_c < num_elements; ++idx_c) {
        for (size_t idx_d = 0; idx_d < num_elements; ++idx_d) {

          const auto &angles_rad = neighbors->dihedrals()[idx_a][idx_b][idx_c][idx_d];
          if (angles_rad.empty()) {
            continue;
          }

          std::string const key = elements[idx_a].symbol + "-" + elements[idx_b].symbol + "-" + elements[idx_c].symbol +
                                  "-" + elements[idx_d].symbol;

          auto &partial_hist = f_dihedral.partials[key];
          if (partial_hist.empty()) {
            partial_hist.assign(num_bins, 0.0);
          }

          processDihedralAngles(angles_rad, partial_hist, bin_width, num_bins, theta_min, theta_max);
        }
      }
    }
  }

  normalizeAndScale(f_dihedral, num_bins, bin_width);
  return f_dihedral;
}

} // namespace correlation::calculators
