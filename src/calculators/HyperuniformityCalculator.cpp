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
#include "math/Precision.hpp"

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

/**
 * @brief Portable 53-bit uniform generator in [0, 1) to guarantee bit-for-bit reproducible random sampling across compilers and platforms.
 */
[[nodiscard]] real_t generate_canonical_portable(std::mt19937_64 &rng) noexcept {
  return static_cast<real_t>(static_cast<double>(rng() >> 11) * (1.0 / 9007199254740992.0));
}
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
  const real_t l_x = lattice_params[0];
  const real_t l_y = lattice_params[1];
  const real_t l_z = lattice_params[2];

  if (l_x <= 0 || l_y <= 0 || l_z <= 0) {
    return {};
  }

  // Half the minimum box length is the maximum window radius
  const real_t l_min = std::min({l_x, l_y, l_z});
  const auto r_max = static_cast<real_t>(l_min / static_cast<real_t>(2.0));
  constexpr real_t r_min = static_cast<real_t>(2.0); // Minimum window radius in Angstroms

  if (r_max <= r_min) {
    return {}; // Box too small for meaningful hyperuniformity analysis
  }

  const auto num_bins = static_cast<size_t>(std::ceil((r_max - r_min) / r_bin_width));
  if (num_bins == 0) {
    return {};
  }

  // Pre-compute bin radii
  std::vector<real_t> radii(num_bins);
  for (size_t k = 0; k < num_bins; ++k) {
    radii[k] = static_cast<real_t>(r_min + (static_cast<real_t>(k) + static_cast<real_t>(0.5)) * r_bin_width);
  }

  // Squared radii for distance comparison
  std::vector<real_t> radii_sq(num_bins);
  for (size_t k = 0; k < num_bins; ++k) {
    const real_t r_edge = r_min + static_cast<real_t>(k + 1) * r_bin_width;
    radii_sq[k] = r_edge * r_edge;
  }

  // Collect atom positions
  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();
  if (num_atoms == 0) {
    return {};
  }

  std::vector<math::Vector3<real_t>> positions;
  positions.reserve(num_atoms);
  for (const auto &atom : atoms) {
    positions.push_back(atom.position());
  }

  // Accumulators: sum of N(R) and sum of N²(R) for each bin
  std::vector<real_t> sum_N(num_bins, static_cast<real_t>(0.0));
  std::vector<real_t> sum_N2(num_bins, static_cast<real_t>(0.0));

  // Deterministic seed for reproducible calculations across runs
  std::mt19937_64 rng(12345); // NOLINT(cert-msc51-cpp, cert-msc32-c, bugprone-random-generator-seed)

  const auto &lattice = cell.latticeVectors();

  std::vector<real_t> dsq_all(num_atoms);

  for (size_t sample_idx = 0; sample_idx < num_samples; ++sample_idx) {
    // Generate a random point in fractional coordinates, then convert to Cartesian
    const auto f_x = generate_canonical_portable(rng);
    const auto f_y = generate_canonical_portable(rng);
    const auto f_z = generate_canonical_portable(rng);

    math::Vector3<real_t> const frac(f_x, f_y, f_z);
    math::Vector3<real_t> const sample_point = lattice * frac;

    // Count atoms within each window radius
    // First compute all minimum-image squared distances
    for (size_t i = 0; i < num_atoms; ++i) {
      math::Vector3<real_t> const diff = positions[i] - sample_point;
      math::Vector3<real_t> const mic = cell.minimumImage(diff);
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
      auto const n_k = static_cast<real_t>(count);
      sum_N[k] += n_k;
      sum_N2[k] += n_k * n_k;
    }
  }

  // Compute statistics
  const auto n_samples = static_cast<real_t>(num_samples);

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
    const real_t mean_N = sum_N[k] / n_samples;
    const real_t mean_N2 = sum_N2[k] / n_samples;
    const real_t raw_var = mean_N2 - mean_N * mean_N;
    const auto variance =
        (raw_var > static_cast<real_t>(1e-12)) ? raw_var : static_cast<real_t>(0.0);

    sigma2_total[k] = variance;
    chi_total[k] = (mean_N > static_cast<real_t>(0.0)) ? (variance / mean_N) : static_cast<real_t>(0.0);
  }

  std::map<std::string, correlation::analysis::Histogram> results;
  results["sigma2_N"] = std::move(sigma2_hist);
  results["chi_H"] = std::move(chi_hist);

  return results;
}

} // namespace correlation::calculators
