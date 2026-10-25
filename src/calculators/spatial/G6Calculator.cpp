/**
 * @file G6Calculator.cpp
 * @brief Implementation of radial bond-orientational correlation function g6(r).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/spatial/G6Calculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "calculators/spatial/DistanceCalculator.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <tbb/enumerable_thread_specific.h>

namespace correlation::calculators {

namespace {

const bool REGISTERED = CalculatorFactory::registerTypeSafe<G6Calculator>("G6Calculator");

struct RegressionResult {
  real_t slope{0.0};
  real_t intercept{0.0};
  real_t r_squared{0.0};
  bool valid{false};
};

[[nodiscard]] RegressionResult computeLinearRegression(std::span<const real_t> x_coords,
                                                       std::span<const real_t> y_coords) noexcept {
  RegressionResult result;
  const size_t num_pts = x_coords.size();
  if (num_pts != y_coords.size() || num_pts < 2) {
    return result;
  }

  real_t sum_x = 0.0;
  real_t sum_y = 0.0;
  real_t sum_xx = 0.0;
  real_t sum_xy = 0.0;
  for (size_t pt_idx = 0; pt_idx < num_pts; ++pt_idx) {
    sum_x += x_coords[pt_idx];
    sum_y += y_coords[pt_idx];
    sum_xx += x_coords[pt_idx] * x_coords[pt_idx];
    sum_xy += x_coords[pt_idx] * y_coords[pt_idx];
  }

  const auto num_pts_f = static_cast<real_t>(num_pts);
  const real_t denom = num_pts_f * sum_xx - sum_x * sum_x;
  if (std::abs(denom) < static_cast<real_t>(1e-12)) {
    return result;
  }

  const real_t slope = (num_pts_f * sum_xy - sum_x * sum_y) / denom;
  const real_t mean_x = sum_x / num_pts_f;
  const real_t mean_y = sum_y / num_pts_f;
  const real_t intercept = mean_y - slope * mean_x;

  real_t ss_res = 0.0;
  real_t ss_tot = 0.0;
  for (size_t pt_idx = 0; pt_idx < num_pts; ++pt_idx) {
    const real_t pred = slope * x_coords[pt_idx] + intercept;
    ss_res += (y_coords[pt_idx] - pred) * (y_coords[pt_idx] - pred);
    ss_tot += (y_coords[pt_idx] - mean_y) * (y_coords[pt_idx] - mean_y);
  }

  result.slope = slope;
  result.intercept = intercept;
  result.r_squared = (ss_tot > static_cast<real_t>(1e-12))
                         ? (static_cast<real_t>(1.0) - ss_res / ss_tot)
                         : static_cast<real_t>(0.0);
  result.valid = true;
  return result;
}

struct FitDataPoints {
  std::vector<real_t> r_pts;
  std::vector<real_t> ln_g6_pts;
};

[[nodiscard]] FitDataPoints extractEnvelopePoints(std::span<const real_t> r_bins,
                                                  std::span<const real_t> g6_values) noexcept {
  const size_t num_bins = std::min(r_bins.size(), g6_values.size());
  FitDataPoints pts;

  // 1. First pass: attempt peak detection
  for (size_t pt_idx = 1; pt_idx + 1 < num_bins; ++pt_idx) {
    const real_t g_curr = g6_values[pt_idx];
    if (g_curr > g6_values[pt_idx - 1] && g_curr >= g6_values[pt_idx + 1] &&
        g_curr > static_cast<real_t>(1e-4)) {
      pts.r_pts.push_back(r_bins[pt_idx]);
      pts.ln_g6_pts.push_back(std::log(g_curr));
    }
  }

  // 2. If fewer than 2 peaks detected, use all strictly positive non-noise values
  if (pts.r_pts.size() < 2) {
    pts.r_pts.clear();
    pts.ln_g6_pts.clear();

    real_t max_val = 0.0;
    for (size_t pt_idx = 0; pt_idx < num_bins; ++pt_idx) {
      max_val = std::max(max_val, g6_values[pt_idx]);
    }
    const real_t threshold =
        std::max(static_cast<real_t>(1e-4), max_val * static_cast<real_t>(0.05));

    for (size_t pt_idx = 0; pt_idx < num_bins; ++pt_idx) {
      if (g6_values[pt_idx] > threshold) {
        pts.r_pts.push_back(r_bins[pt_idx]);
        pts.ln_g6_pts.push_back(std::log(g6_values[pt_idx]));
      }
    }
  }

  return pts;
}

[[nodiscard]] inline real_t computePairCoherence(const std::complex<real_t> *q6m_i,
                                                 const std::complex<real_t> *q6m_j) noexcept {
  auto chi = static_cast<real_t>(0.0);
  for (size_t mi = 0; mi < 13; ++mi) {
    chi += q6m_i[mi].real() * q6m_j[mi].real() + q6m_i[mi].imag() * q6m_j[mi].imag();
  }
  return chi;
}

struct ThreadLocalG6 {
  std::vector<real_t> sum_chi;
  std::vector<size_t> pair_counts;

  explicit ThreadLocalG6(size_t num_bins) : sum_chi(num_bins, 0.0), pair_counts(num_bins, 0) {}
};

[[nodiscard]] const std::vector<std::complex<real_t>> *
resolveQ6mVectors(const correlation::core::Cell &cell,
                  const correlation::analysis::StructureAnalyzer *neighbors,
                  const SteinhardtCalculator::SteinhardtParams *steinhardt_params,
                  std::optional<SteinhardtCalculator::SteinhardtParams> &local_steinhardt) {
  if (steinhardt_params != nullptr && steinhardt_params->q6m_normalized.has_value()) {
    return &(*steinhardt_params->q6m_normalized);
  }
  local_steinhardt = SteinhardtCalculator::calculatePerAtom(cell, neighbors, {.keep_q6m = true});
  if (local_steinhardt->q6m_normalized.has_value()) {
    return &(*local_steinhardt->q6m_normalized);
  }
  return nullptr;
}

struct G6GridParams {
  real_t r_max{0.0};
  real_t r_bin_width{0.0};
  size_t num_bins{0};
};

struct AccumulatedG6Data {
  std::vector<real_t> g6_values;
  std::vector<real_t> pair_counts;
};

[[nodiscard]] AccumulatedG6Data
accumulateG6(const correlation::core::Cell &cell,
             const correlation::analysis::StructureAnalyzer *neighbors,
             std::span<const std::complex<real_t>> q6m_norm, const G6GridParams &grid_params) {
  const auto &neighbor_graph = neighbors->neighborGraph();
  const size_t num_atoms = cell.atomCount();

  tbb::enumerable_thread_specific<ThreadLocalG6> ets(
      [num_bins = grid_params.num_bins] { return ThreadLocalG6(num_bins); });

  DistanceCalculator::forEachPairWithin(
      cell, grid_params.r_max, true,
      [&](size_t orig_idx, size_t j_idx, real_t d_sq,
          const correlation::math::Vector3<real_t> & /*r_ij*/) {
        if (orig_idx >= num_atoms || j_idx >= num_atoms) {
          return;
        }

        if (neighbor_graph.getNeighbors(orig_idx).size() < 2 ||
            neighbor_graph.getNeighbors(j_idx).size() < 2) {
          return;
        }

        const real_t dist = std::sqrt(d_sq);
        const auto bin_idx = static_cast<size_t>(dist / grid_params.r_bin_width);
        if (bin_idx >= grid_params.num_bins) {
          return;
        }

        const real_t chi = computePairCoherence(&q6m_norm[orig_idx * 13], &q6m_norm[j_idx * 13]);

        auto &local = ets.local();
        local.sum_chi[bin_idx] += chi;
        local.pair_counts[bin_idx] += 1;
      });

  AccumulatedG6Data data{
      .g6_values = std::vector<real_t>(grid_params.num_bins, 0.0),
      .pair_counts = std::vector<real_t>(grid_params.num_bins, 0.0),
  };

  for (const auto &local : ets) {
    for (size_t bin_idx = 0; bin_idx < grid_params.num_bins; ++bin_idx) {
      data.g6_values[bin_idx] += local.sum_chi[bin_idx];
      data.pair_counts[bin_idx] += static_cast<real_t>(local.pair_counts[bin_idx]);
    }
  }

  for (size_t bin_idx = 0; bin_idx < grid_params.num_bins; ++bin_idx) {
    if (data.pair_counts[bin_idx] > 0.0) {
      data.g6_values[bin_idx] /= data.pair_counts[bin_idx];
    }
  }

  return data;
}

} // namespace

CorrelationLengthFit
G6Calculator::fitCorrelationLength(std::span<const real_t> r_bins,
                                   std::span<const real_t> g6_values) noexcept {
  const auto pts = extractEnvelopePoints(r_bins, g6_values);
  if (pts.r_pts.size() < 2) {
    return CorrelationLengthFit{};
  }

  const auto reg = computeLinearRegression(pts.r_pts, pts.ln_g6_pts);
  if (!reg.valid) {
    return CorrelationLengthFit{};
  }

  if (reg.slope < static_cast<real_t>(-1e-4)) {
    const auto xi_decay = -static_cast<real_t>(1.0) / reg.slope;
    return CorrelationLengthFit{
        .xi_6 = xi_decay,
        .r_squared = std::clamp(reg.r_squared, static_cast<real_t>(0.0), static_cast<real_t>(1.0)),
        .is_plateau = false,
    };
  }

  return CorrelationLengthFit{
      .xi_6 = std::numeric_limits<real_t>::infinity(),
      .r_squared = std::clamp(reg.r_squared, static_cast<real_t>(0.0), static_cast<real_t>(1.0)),
      .is_plateau = true,
  };
}

void G6Calculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                  const correlation::analysis::AnalysisSettings &settings) const {
  if (dists.neighbors() == nullptr) {
    return;
  }
  const real_t r_max = settings.r_max > 0.0 ? settings.r_max : static_cast<real_t>(10.0);
  const real_t r_bin_width =
      settings.r_bin_width > 0.0 ? settings.r_bin_width : static_cast<real_t>(0.05);

  auto hist = calculate(dists.cell(), dists.neighbors(), r_max, r_bin_width);
  dists.addHistogram("G6_r", std::move(hist));
}

correlation::analysis::Histogram
G6Calculator::calculate(const correlation::core::Cell &cell,
                        const correlation::analysis::StructureAnalyzer *neighbors, real_t r_max,
                        real_t r_bin_width,
                        const SteinhardtCalculator::SteinhardtParams *steinhardt_params) {
  if (neighbors == nullptr) {
    throw std::logic_error("Cannot calculate G6. Neighbor list has not been computed.");
  }
  if (r_max <= static_cast<real_t>(0.0)) {
    throw std::invalid_argument("r_max must be strictly positive.");
  }
  if (r_bin_width <= static_cast<real_t>(0.0)) {
    throw std::invalid_argument("r_bin_width must be strictly positive.");
  }

  const size_t num_atoms = cell.atomCount();
  const size_t num_bins = std::max(size_t{1}, static_cast<size_t>(std::ceil(r_max / r_bin_width)));

  std::vector<real_t> bins(num_bins);
  for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
    bins[bin_idx] = (static_cast<real_t>(bin_idx) + static_cast<real_t>(0.5)) * r_bin_width;
  }

  correlation::analysis::Histogram hist;
  hist.title = "g₆(r) — Bond-Orientational Correlation Function";
  hist.x_label = "r";
  hist.y_label = "g₆(r)";
  hist.bins = bins;

  if (num_atoms < 2) {
    hist.partials["Total"] = std::vector<real_t>(num_bins, 0.0);
    hist.partials["Counts"] = std::vector<real_t>(num_bins, 0.0);
    return hist;
  }

  std::optional<SteinhardtCalculator::SteinhardtParams> local_steinhardt;
  const auto *q6m_norm_ptr =
      resolveQ6mVectors(cell, neighbors, steinhardt_params, local_steinhardt);

  if (q6m_norm_ptr == nullptr || q6m_norm_ptr->size() < num_atoms * 13) {
    hist.partials["Total"] = std::vector<real_t>(num_bins, 0.0);
    hist.partials["Counts"] = std::vector<real_t>(num_bins, 0.0);
    return hist;
  }

  const G6GridParams grid_params{
      .r_max = r_max,
      .r_bin_width = r_bin_width,
      .num_bins = num_bins,
  };
  auto data = accumulateG6(cell, neighbors, *q6m_norm_ptr, grid_params);

  const auto fit = fitCorrelationLength(bins, data.g6_values);
  if (fit.is_plateau) {
    hist.description = "Crystalline plateau (xi_6 -> inf)";
  } else if (fit.xi_6 > 0.0) {
    hist.description = std::format("Decay xi_6 = {:.3f} A (R^2 = {:.3f})", fit.xi_6, fit.r_squared);
  } else {
    hist.description = "No correlation fit";
  }

  hist.file_suffix = "g6_r";
  hist.x_unit = "Å";
  hist.y_unit = "";
  hist.partials["Total"] = std::move(data.g6_values);
  hist.partials["Counts"] = std::move(data.pair_counts);

  return hist;
}

} // namespace correlation::calculators
