/**
 * @file HyperuniformityCalculator.cpp
 * @brief Implementation of the hyperuniformity (local number variance) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/HyperuniformityCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/LinearAlgebra.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
// NOLINTNEXTLINE(cert-err58-cpp)
const bool registered = CalculatorFactory::registerTypeSafe<HyperuniformityCalculator>("HyperuniformityCalculator");
} // namespace

void HyperuniformityCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                               const correlation::analysis::AnalysisSettings &settings) const {
  auto results =
      calculate(dists.cell(), {.num_samples = settings.hyperuniformity_samples, .r_bin_width = settings.r_bin_width});
  for (auto &[name, histogram] : results) {
    dists.addHistogram(name, std::move(histogram));
  }
}

std::map<std::string, correlation::analysis::Histogram>
HyperuniformityCalculator::calculate(const correlation::core::Cell &cell, const HyperuniformityParams &params) {
  const auto num_samples = params.num_samples;
  const auto r_bin_width = params.r_bin_width;

  if (r_bin_width <= 0) {
    throw std::invalid_argument("Bin width must be positive, got: " + std::to_string(r_bin_width));
  }
  if (num_samples == 0) {
    throw std::invalid_argument("Number of samples must be positive.");
  }

  const auto &lattice_params = cell.lattice_parameters();
  const double l_x = lattice_params[0];
  const double l_y = lattice_params[1];
  const double l_z = lattice_params[2];

  if (l_x <= 0 || l_y <= 0 || l_z <= 0) {
    return {};
  }

  // Half the minimum box length is the maximum window radius
  const double l_min = std::min({l_x, l_y, l_z});
  const double r_max = l_min / 2.0;
  constexpr double r_min = 2.0; // Minimum window radius in Angstroms

  if (r_max <= r_min) {
    return {}; // Box too small for meaningful hyperuniformity analysis
  }

  const auto num_bins = static_cast<size_t>(std::ceil((r_max - r_min) / r_bin_width));
  if (num_bins == 0) {
    return {};
  }

  // Pre-compute bin radii
  std::vector<double> radii(num_bins);
  for (size_t k = 0; k < num_bins; ++k) {
    radii[k] = r_min + (static_cast<double>(k) + 0.5) * r_bin_width;
  }

  // Squared radii for distance comparison
  std::vector<double> radii_sq(num_bins);
  for (size_t k = 0; k < num_bins; ++k) {
    double r_edge = r_min + static_cast<double>(k + 1) * r_bin_width;
    radii_sq[k] = r_edge * r_edge;
  }

  // Collect atom positions
  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();
  if (num_atoms == 0) {
    return {};
  }

  std::vector<math::Vector3<double>> positions;
  positions.reserve(num_atoms);
  for (const auto &atom : atoms) {
    positions.push_back(atom.position());
  }

  // Accumulators: sum of N(R) and sum of N²(R) for each bin
  std::vector<double> sum_N(num_bins, 0.0);
  std::vector<double> sum_N2(num_bins, 0.0);

  // Seed generated once per process from std::random_device — same across frames
  // within a run, but different between program executions.
  static const auto seed = std::random_device{}();
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  const auto &lattice = cell.latticeVectors();

  for (size_t sample_idx = 0; sample_idx < num_samples; ++sample_idx) {
    // Generate a random point in fractional coordinates, then convert to Cartesian
    const double f_x = dist(rng);
    const double f_y = dist(rng);
    const double f_z = dist(rng);

    math::Vector3<double> const frac(f_x, f_y, f_z);
    math::Vector3<double> const sample_point = lattice * frac;

    // Count atoms within each window radius
    // First compute all minimum-image squared distances
    std::vector<double> dsq_all(num_atoms);
    for (size_t i = 0; i < num_atoms; ++i) {
      math::Vector3<double> const diff = positions[i] - sample_point;
      math::Vector3<double> const mic = cell.minimumImage(diff);
      dsq_all[i] = math::norm_sq(mic);
    }

    // Sort squared distances for efficient counting
    std::ranges::sort(dsq_all);

    // Count atoms in windows of increasing radius using sorted distances
    size_t count = 0;
    size_t atom_idx = 0;
    for (size_t k = 0; k < num_bins; ++k) {
      while (atom_idx < num_atoms && dsq_all[atom_idx] <= radii_sq[k]) {
        ++count;
        ++atom_idx;
      }
      auto const n_k = static_cast<double>(count);
      sum_N[k] += n_k;
      sum_N2[k] += n_k * n_k;
    }
  }

  // Compute statistics
  const auto n_samples = static_cast<double>(num_samples);

  // σ²_N(R) histogram
  correlation::analysis::Histogram sigma2_hist;
  sigma2_hist.bins = radii;
  sigma2_hist.x_label = "R";
  sigma2_hist.title = "σ²_N(R) — Local Number Variance";
  sigma2_hist.y_label = "σ²_N(R)";
  sigma2_hist.x_unit = "Å";
  sigma2_hist.y_unit = "";
  sigma2_hist.description = "Local Number Variance";
  sigma2_hist.file_suffix = "_sigma2_N";

  auto &sigma2_total = sigma2_hist.partials["Total"];
  sigma2_total.resize(num_bins);

  // χ_H(R) histogram
  correlation::analysis::Histogram chi_hist;
  chi_hist.bins = radii;
  chi_hist.x_label = "R";
  chi_hist.title = "χ_H(R) — Hyperuniformity Index";
  chi_hist.y_label = "χ_H(R)";
  chi_hist.x_unit = "Å";
  chi_hist.y_unit = "";
  chi_hist.description = "Hyperuniformity Index";
  chi_hist.file_suffix = "_chi_H";

  auto &chi_total = chi_hist.partials["Total"];
  chi_total.resize(num_bins);

  for (size_t k = 0; k < num_bins; ++k) {
    const double mean_N = sum_N[k] / n_samples;
    const double mean_N2 = sum_N2[k] / n_samples;
    const double variance = mean_N2 - mean_N * mean_N;

    sigma2_total[k] = variance;
    chi_total[k] = (mean_N > 0.0) ? variance / mean_N : 0.0;
  }

  std::map<std::string, correlation::analysis::Histogram> results;
  results["sigma2_N"] = std::move(sigma2_hist);
  results["chi_H"] = std::move(chi_hist);

  return results;
}

} // namespace correlation::calculators
