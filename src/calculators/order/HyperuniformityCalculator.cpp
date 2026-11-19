/**
 * @file HyperuniformityCalculator.cpp
 * @brief Implementation of the hyperuniformity (local number variance) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/order/HyperuniformityCalculator.hpp"
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
const bool REGISTERED =
    CalculatorFactory::registerTypeSafe<HyperuniformityCalculator>("HyperuniformityCalculator");

/**
 * @brief Portable 53-bit uniform generator in [0, 1) to guarantee bit-for-bit reproducible random
 * sampling across compilers and platforms.
 */
[[nodiscard]] real_t generateCanonicalPortable(std::mt19937_64 &rng) noexcept {
  return static_cast<real_t>(rng() >> 11U) * static_cast<real_t>(1.0 / 9007199254740992.0);
}
} // namespace

void HyperuniformityCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &dists,
    const correlation::analysis::AnalysisSettings &settings) const {
  auto results = calculate(dists.cell(), {
                                             .num_samples = settings.hyperuniformity_samples,
                                             .r_bin_width = settings.r_bin_width,
                                         });
  for (auto &[name, histogram] : results) {
    dists.addHistogram(name, std::move(histogram));
  }
}

std::map<std::string, correlation::analysis::Histogram>
HyperuniformityCalculator::calculate(const correlation::core::Cell &cell,
                                     const HyperuniformityParams &params) {
  const auto num_samples = params.num_samples;
  const auto r_bin_width = params.r_bin_width;

  if (r_bin_width <= 0) {
    throw std::invalid_argument("Bin width must be positive, got: " + std::to_string(r_bin_width));
  }
  if (num_samples == 0) {
    throw std::invalid_argument("Number of samples must be positive.");
  }

  const auto &lattice_params = cell.latticeParameters();
  const real_t l_x = lattice_params[0];
  const real_t l_y = lattice_params[1];
  const real_t l_z = lattice_params[2];

  if (l_x <= 0 || l_y <= 0 || l_z <= 0) {
    return {};
  }

  // Half the minimum box length is the maximum window radius
  const real_t l_min = std::min({l_x, l_y, l_z});
  const auto r_max = (l_min / static_cast<real_t>(2.0));
  constexpr auto R_MIN = static_cast<real_t>(2.0); // Minimum window radius in Angstroms

  if (r_max <= R_MIN) {
    return {}; // Box too small for meaningful hyperuniformity analysis
  }

  const auto num_bins = static_cast<size_t>(std::ceil((r_max - R_MIN) / r_bin_width));
  if (num_bins == 0) {
    return {};
  }

  // Pre-compute bin radii
  std::vector<real_t> radii(num_bins);
  for (size_t k = 0; k < num_bins; ++k) {
    radii[k] = (R_MIN + (static_cast<real_t>(k) + static_cast<real_t>(0.5)) * r_bin_width);
  }

  // Squared radii for distance comparison
  std::vector<real_t> radii_sq(num_bins);
  for (size_t k = 0; k < num_bins; ++k) {
    const real_t r_edge = R_MIN + static_cast<real_t>(k + 1) * r_bin_width;
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
  std::vector<real_t> sum_n(num_bins, static_cast<real_t>(0.0));
  std::vector<real_t> sum_n2(num_bins, static_cast<real_t>(0.0));

  // Deterministic seed for reproducible calculations across runs
  std::seed_seq seed{12345};
  std::mt19937_64 rng(seed);

  const auto &lattice = cell.latticeVectors();

  std::vector<real_t> dsq_all(num_atoms);

  for (size_t sample_idx = 0; sample_idx < num_samples; ++sample_idx) {
    // Generate a random point in fractional coordinates, then convert to Cartesian
    const auto f_x = generateCanonicalPortable(rng);
    const auto f_y = generateCanonicalPortable(rng);
    const auto f_z = generateCanonicalPortable(rng);

    math::Vector3<real_t> const frac(f_x, f_y, f_z);
    math::Vector3<real_t> const sample_point = lattice * frac;

    // Count atoms within each window radius
    // First compute all minimum-image squared distances
    for (size_t i = 0; i < num_atoms; ++i) {
      math::Vector3<real_t> const diff = positions[i] - sample_point;
      math::Vector3<real_t> const mic = cell.minimumImage(diff);
      dsq_all[i] = math::normSq(mic);
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
      sum_n[k] += n_k;
      sum_n2[k] += n_k * n_k;
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
    const real_t mean_n = sum_n[k] / n_samples;
    const real_t mean_n2 = sum_n2[k] / n_samples;
    const real_t raw_var = mean_n2 - mean_n * mean_n;
    const auto variance =
        (raw_var > static_cast<real_t>(1e-12)) ? raw_var : static_cast<real_t>(0.0);

    sigma2_total[k] = variance;
    chi_total[k] =
        (mean_n > static_cast<real_t>(0.0)) ? (variance / mean_n) : static_cast<real_t>(0.0);
  }

  std::map<std::string, correlation::analysis::Histogram> results;
  results["sigma2_N"] = std::move(sigma2_hist);
  results["chi_H"] = std::move(chi_hist);

  return results;
}

namespace {

[[nodiscard]] double fitLogLogSlope(std::span<const double> x_vals,
                                    std::span<const double> y_vals) noexcept {
  const size_t count = std::min(x_vals.size(), y_vals.size());
  if (count < 2U) {
    return 0.0;
  }

  double sum_x = 0.0;
  double sum_y = 0.0;
  for (size_t i = 0; i < count; ++i) {
    sum_x += x_vals[i];
    sum_y += y_vals[i];
  }
  const double mean_x = sum_x / static_cast<double>(count);
  const double mean_y = sum_y / static_cast<double>(count);

  double cov_xy = 0.0;
  double var_x = 0.0;
  for (size_t i = 0; i < count; ++i) {
    const double delta_x = x_vals[i] - mean_x;
    cov_xy += delta_x * (y_vals[i] - mean_y);
    var_x += delta_x * delta_x;
  }
  return (var_x > 1e-18) ? (cov_xy / var_x) : 0.0;
}

[[nodiscard]] double extrapolateZeroIntercept(std::span<const double> x_vals,
                                              std::span<const double> y_vals) noexcept {
  const size_t count = std::min(x_vals.size(), y_vals.size());
  if (count == 0U) {
    return 0.0;
  }
  if (count == 1U) {
    return std::max(0.0, y_vals[0]);
  }

  double sum_x = 0.0;
  double sum_y = 0.0;
  for (size_t i = 0; i < count; ++i) {
    sum_x += x_vals[i];
    sum_y += y_vals[i];
  }
  const double mean_x = sum_x / static_cast<double>(count);
  const double mean_y = sum_y / static_cast<double>(count);

  double cov_xy = 0.0;
  double var_x = 0.0;
  for (size_t i = 0; i < count; ++i) {
    const double delta_x = x_vals[i] - mean_x;
    cov_xy += delta_x * (y_vals[i] - mean_y);
    var_x += delta_x * delta_x;
  }
  const double slope = (var_x > 1e-18) ? (cov_xy / var_x) : 0.0;
  const double intercept = mean_y - slope * mean_x;
  if (intercept > 0.0) {
    return intercept;
  }
  return std::max(0.0, y_vals[0]);
}

} // namespace

correlation::analysis::HyperuniformityDiagnostics
HyperuniformityCalculator::computeRealSpaceDiagnostics(std::span<const real_t> radii,
                                                       std::span<const real_t> sigma2_n,
                                                       std::span<const real_t> chi_h) noexcept {
  correlation::analysis::HyperuniformityDiagnostics diag{};
  const size_t bin_count = std::min(radii.size(), sigma2_n.size());
  if (bin_count < 3U) {
    return diag;
  }

  // Skip the smallest 30% of radii when enough bins exist to avoid short-range discreteness noise.
  const size_t fit_start = (bin_count >= 6U) ? (bin_count * 3U / 10U) : 0U;
  std::vector<double> log_r;
  std::vector<double> log_var;
  log_r.reserve(bin_count - fit_start);
  log_var.reserve(bin_count - fit_start);

  for (size_t i = fit_start; i < bin_count; ++i) {
    const auto r_val = static_cast<double>(radii[i]);
    const auto v_val = static_cast<double>(sigma2_n[i]);
    if (r_val > 0.0 && v_val > 1e-12) {
      log_r.push_back(std::log(r_val));
      log_var.push_back(std::log(v_val));
    }
  }

  if (log_r.size() >= 3U) {
    diag.variance_exponent_beta = static_cast<real_t>(fitLogLogSlope(log_r, log_var));
    diag.has_real_space = true;
  }

  const size_t chi_n = std::min(radii.size(), chi_h.size());
  if (chi_n > 0U) {
    const size_t tail_start = (chi_n >= 4U) ? (chi_n * 3U / 4U) : 0U;
    double chi_sum = 0.0;
    size_t chi_count = 0U;
    for (size_t i = tail_start; i < chi_n; ++i) {
      const auto c_val = static_cast<double>(chi_h[i]);
      if (std::isfinite(c_val) && c_val >= 0.0) {
        chi_sum += c_val;
        ++chi_count;
      }
    }
    if (chi_count > 0U) {
      diag.asymptotic_chi_h = static_cast<real_t>(chi_sum / static_cast<double>(chi_count));
      diag.has_real_space = true;
    }
  }

  return diag;
}

correlation::analysis::HyperuniformityDiagnostics
HyperuniformityCalculator::computeReciprocalDiagnostics(std::span<const real_t> q_bins,
                                                        std::span<const real_t> sq_total) noexcept {
  correlation::analysis::HyperuniformityDiagnostics diag{};
  const size_t bin_count = std::min(q_bins.size(), sq_total.size());
  if (bin_count < 4U) {
    return diag;
  }

  // Locate primary diffraction peak S(q_peak) for positive q.
  double s_peak = 0.0;
  double q_peak = 0.0;
  size_t peak_idx = 0U;
  for (size_t i = 0; i < bin_count; ++i) {
    const auto q_val = static_cast<double>(q_bins[i]);
    const auto s_val = static_cast<double>(sq_total[i]);
    if (q_val > 0.0 && std::isfinite(s_val) && s_val > s_peak) {
      s_peak = s_val;
      q_peak = q_val;
      peak_idx = i;
    }
  }
  if (s_peak <= 1e-12 || peak_idx == 0U) {
    return diag;
  }

  // Gather positive low-q points below 0.5 * q_peak (with fallback to lowest pre-peak points).
  const double q_cutoff = 0.5 * q_peak;
  std::vector<double> low_q;
  std::vector<double> low_s;
  std::vector<double> log_q;
  std::vector<double> log_s;
  for (size_t i = 0; i < peak_idx; ++i) {
    const auto q_val = static_cast<double>(q_bins[i]);
    const auto s_val = static_cast<double>(sq_total[i]);
    if (q_val <= 0.0 || !std::isfinite(s_val) || s_val <= 0.0) {
      continue;
    }
    if (q_val <= q_cutoff || low_q.size() < 4U) {
      low_q.push_back(q_val);
      low_s.push_back(s_val);
      log_q.push_back(std::log(q_val));
      log_s.push_back(std::log(s_val));
    }
  }

  if (low_q.empty()) {
    return diag;
  }

  const size_t extrap_n = std::min<size_t>(low_q.size(), 5U);
  const double s_zero = extrapolateZeroIntercept(std::span<const double>(low_q).first(extrap_n),
                                                 std::span<const double>(low_s).first(extrap_n));
  diag.torquato_h_index = static_cast<real_t>(s_zero / s_peak);
  if (log_q.size() >= 3U) {
    diag.low_q_exponent_alpha = static_cast<real_t>(fitLogLogSlope(log_q, log_s));
  }
  diag.has_reciprocal_space = true;
  return diag;
}

HyperuniformityVerdict HyperuniformityCalculator::classifyHyperuniformity(
    const correlation::analysis::HyperuniformityDiagnostics &diag) noexcept {
  if (!diag.has_real_space && !diag.has_reciprocal_space) {
    return {.label = "", .severity_level = -1};
  }

  const bool strong_h =
      diag.has_reciprocal_space && diag.torquato_h_index <= static_cast<real_t>(1e-4);
  if (strong_h) {
    return {.label = "Strongly Hyperuniform", .severity_level = 1};
  }

  const bool eff_hu_h =
      diag.has_reciprocal_space && diag.torquato_h_index <= static_cast<real_t>(1e-3);
  const bool eff_hu_beta = diag.has_real_space &&
                           diag.variance_exponent_beta > static_cast<real_t>(0.0) &&
                           diag.variance_exponent_beta <= static_cast<real_t>(2.35);
  if (eff_hu_h || eff_hu_beta) {
    return {.label = "Effectively Hyperuniform", .severity_level = 1};
  }

  const bool near_hu_h =
      diag.has_reciprocal_space && diag.torquato_h_index <= static_cast<real_t>(1e-2);
  const bool near_hu_beta = diag.has_real_space &&
                            diag.variance_exponent_beta > static_cast<real_t>(0.0) &&
                            diag.variance_exponent_beta <= static_cast<real_t>(2.70);
  if (near_hu_h || near_hu_beta) {
    return {.label = "Near-Hyperuniform", .severity_level = 2};
  }

  return {.label = "Non-Hyperuniform Liquid", .severity_level = 0};
}

} // namespace correlation::calculators
