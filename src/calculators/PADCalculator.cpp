/**
 * @file PADCalculator.cpp
 * @brief Implementation of the plane angle distribution calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "calculators/PADCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"
#include "math/SIMDUtils.hpp"

#include <stdexcept>

namespace correlation::calculators {

namespace {
bool registered = CalculatorFactory::instance().registerCalculator(
    std::make_unique<PADCalculator>());
} // namespace

void PADCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &df,
    const correlation::analysis::AnalysisSettings &settings) const {
  df.addHistogram(
      "BAD", calculate(df.cell(), df.neighbors(), settings.angle_bin_width));
}

correlation::analysis::Histogram PADCalculator::calculate(
    const correlation::core::Cell &cell,
    const correlation::analysis::StructureAnalyzer *neighbors,
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

  correlation::analysis::Histogram f_theta;
  f_theta.x_label = "θ";
  f_theta.title = "Plane-Angle Distribution";
  f_theta.y_label = "P(θ)";
  f_theta.x_unit = "°";
  f_theta.y_unit = "°⁻¹";
  f_theta.description = "Bond Angle Distribution";
  f_theta.file_suffix = "_PAD";
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
          double angle_deg = angle_rad * correlation::math::rad_to_deg;

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
    correlation::math::scale_bins(partial.data(), normalization_factor,
                                  num_bins);
  }
  return f_theta;
}

} // namespace correlation::calculators
