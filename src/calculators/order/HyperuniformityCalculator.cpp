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
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
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

struct ThreadLocalWelford {
  std::vector<double> mean;
  std::vector<double> m2;
  std::vector<size_t> bin_counts;
  size_t count{0};
};

struct HyperuniformityGrid {
  real_t r_min;
  real_t r_bin_width;
  real_t inv_dr;
  real_t max_edge_sq;
  size_t num_bins;
};

[[nodiscard]] std::vector<math::Vector3<real_t>>
generateSamplePoints(const correlation::core::Cell &cell, size_t num_samples) {
  std::seed_seq seed{12345};
  std::mt19937_64 rng(seed);
  const auto &lattice = cell.latticeVectors();
  std::vector<math::Vector3<real_t>> sample_points;
  sample_points.reserve(num_samples);
  for (size_t sample_idx = 0; sample_idx < num_samples; ++sample_idx) {
    const auto f_x = generateCanonicalPortable(rng);
    const auto f_y = generateCanonicalPortable(rng);
    const auto f_z = generateCanonicalPortable(rng);
    sample_points.push_back(lattice * math::Vector3<real_t>(f_x, f_y, f_z));
  }
  return sample_points;
}

void accumulateSampleWelford(const correlation::core::Cell &cell,
                             std::span<const math::Vector3<real_t>> positions,
                             const math::Vector3<real_t> &sample_point,
                             const HyperuniformityGrid &grid, ThreadLocalWelford &local) {
  std::ranges::fill(local.bin_counts, 0U);

  // Direct O(N_atoms) binning instead of O(N log N) sort
  for (const auto &pos : positions) {
    const math::Vector3<real_t> mic = cell.minimumImage(pos - sample_point);
    const real_t d_sq = math::normSq(mic);
    if (d_sq > grid.max_edge_sq) {
      continue;
    }
    const real_t dist = std::sqrt(d_sq);
    size_t bin = 0;
    if (dist > grid.r_min) {
      bin = static_cast<size_t>((dist - grid.r_min) * grid.inv_dr);
      if (bin > 0 && dist <= grid.r_min + static_cast<real_t>(bin) * grid.r_bin_width) {
        --bin;
      }
    }
    if (bin < grid.num_bins) {
      ++local.bin_counts[bin];
    }
  }

  ++local.count;
  const auto count_d = static_cast<double>(local.count);
  size_t running_atoms = 0;
  for (size_t k = 0; k < grid.num_bins; ++k) {
    running_atoms += local.bin_counts[k];
    const auto val = static_cast<double>(running_atoms);
    const double delta = val - local.mean[k];
    local.mean[k] += delta / count_d;
    const double delta2 = val - local.mean[k];
    local.m2[k] += delta * delta2;
  }
}

struct MergedWelford {
  std::vector<double> mean;
  std::vector<double> m2;
  size_t count{0};
};

[[nodiscard]] MergedWelford
mergeWelfordStatistics(const tbb::enumerable_thread_specific<ThreadLocalWelford> &ets,
                       size_t num_bins) {
  MergedWelford merged{
      .mean = std::vector<double>(num_bins, 0.0),
      .m2 = std::vector<double>(num_bins, 0.0),
      .count = 0U,
  };
  for (const auto &local : ets) {
    if (local.count == 0) {
      continue;
    }
    const size_t combined_count = merged.count + local.count;
    const auto n_a = static_cast<double>(merged.count);
    const auto n_b = static_cast<double>(local.count);
    const auto n_ab = static_cast<double>(combined_count);
    for (size_t k = 0; k < num_bins; ++k) {
      const double delta = local.mean[k] - merged.mean[k];
      merged.mean[k] += delta * (n_b / n_ab);
      merged.m2[k] += local.m2[k] + delta * delta * (n_a * n_b / n_ab);
    }
    merged.count = combined_count;
  }
  return merged;
}

[[nodiscard]] std::map<std::string, correlation::analysis::Histogram>
buildHyperuniformityHistograms(const std::vector<real_t> &radii, const MergedWelford &merged) {
  const size_t num_bins = radii.size();

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

  const double inv_total_count =
      (merged.count > 0) ? (1.0 / static_cast<double>(merged.count)) : 0.0;
  for (size_t k = 0; k < num_bins; ++k) {
    const double mean_n = merged.mean[k];
    const double raw_var = merged.m2[k] * inv_total_count;
    const auto variance =
        (raw_var > 1e-12) ? static_cast<real_t>(raw_var) : static_cast<real_t>(0.0);

    sigma2_total[k] = variance;
    chi_total[k] = (mean_n > 0.0) ? static_cast<real_t>(static_cast<double>(variance) / mean_n)
                                  : static_cast<real_t>(0.0);
  }

  std::map<std::string, correlation::analysis::Histogram> results;
  results["sigma2_N"] = std::move(sigma2_hist);
  results["chi_H"] = std::move(chi_hist);
  return results;
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

  // Half the minimum perpendicular cell width prevents sphere self-overlap in triclinic boxes
  const auto widths = cell.perpendicularWidths();
  const real_t l_min = std::min({widths[0], widths[1], widths[2]});
  if (l_min <= static_cast<real_t>(0.0)) {
    return {};
  }
  const auto r_max = (l_min / static_cast<real_t>(2.0));
  constexpr auto R_MIN = static_cast<real_t>(2.0); // Minimum window radius in Angstroms

  if (r_max <= R_MIN) {
    return {}; // Box too small for meaningful hyperuniformity analysis
  }

  const auto num_bins = static_cast<size_t>(std::ceil((r_max - R_MIN) / r_bin_width));
  if (num_bins == 0) {
    return {};
  }

  const auto &atoms = cell.atoms();
  if (atoms.empty()) {
    return {};
  }

  std::vector<real_t> radii(num_bins);
  for (size_t k = 0; k < num_bins; ++k) {
    radii[k] = (R_MIN + (static_cast<real_t>(k) + static_cast<real_t>(0.5)) * r_bin_width);
  }
  const real_t max_edge = R_MIN + static_cast<real_t>(num_bins) * r_bin_width;
  const HyperuniformityGrid grid{
      .r_min = R_MIN,
      .r_bin_width = r_bin_width,
      .inv_dr = static_cast<real_t>(1.0) / r_bin_width,
      .max_edge_sq = max_edge * max_edge,
      .num_bins = num_bins,
  };

  std::vector<math::Vector3<real_t>> positions;
  positions.reserve(atoms.size());
  for (const auto &atom : atoms) {
    positions.push_back(atom.position());
  }

  const auto sample_points = generateSamplePoints(cell, num_samples);

  tbb::enumerable_thread_specific<ThreadLocalWelford> ets([num_bins] {
    return ThreadLocalWelford{
        .mean = std::vector<double>(num_bins, 0.0),
        .m2 = std::vector<double>(num_bins, 0.0),
        .bin_counts = std::vector<size_t>(num_bins, 0U),
        .count = 0U,
    };
  });

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, num_samples), [&](const tbb::blocked_range<size_t> &range) {
        auto &local = ets.local();
        for (size_t sample_idx = range.begin(); sample_idx != range.end(); ++sample_idx) {
          accumulateSampleWelford(cell, positions, sample_points[sample_idx], grid, local);
        }
      });

  const auto merged = mergeWelfordStatistics(ets, num_bins);
  return buildHyperuniformityHistograms(radii, merged);
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
