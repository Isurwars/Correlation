/**
 * @file PeakIntegrationService.cpp
 * @brief Implementation of numerical integration, peak detection, and coordination numbers.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/PeakIntegrationService.hpp"
#include "math/Constants.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace correlation::analysis {

namespace {

constexpr real_t MIN_PEAK_THRESHOLD_FRAC = 0.05;
using correlation::math::FOUR_PI;

[[nodiscard]] std::pair<std::size_t, std::size_t>
findBoundingIndices(std::span<const real_t> x_bins, real_t r_min, real_t r_max) noexcept {
  const auto n = x_bins.size();
  std::size_t min_idx = 0;
  while (min_idx < n && x_bins[min_idx] < r_min) {
    ++min_idx;
  }
  std::size_t max_idx = min_idx;
  while (max_idx < n && x_bins[max_idx] <= r_max) {
    ++max_idx;
  }
  if (max_idx > 0 && max_idx == n) {
    max_idx = n - 1;
  }
  return {min_idx, max_idx};
}

[[nodiscard]] std::size_t locatePeakIndex(std::span<const real_t> y_vals, std::size_t min_idx,
                                          std::size_t max_idx) noexcept {
  std::size_t peak_idx = min_idx;
  real_t max_y = std::numeric_limits<real_t>::lowest();
  for (std::size_t i = min_idx; i <= max_idx && i < y_vals.size(); ++i) {
    if (y_vals[i] > max_y) {
      max_y = y_vals[i];
      peak_idx = i;
    }
  }
  return peak_idx;
}

[[nodiscard]] real_t interpolateHalfCross(real_t x1, real_t y1, real_t x2, real_t y2,
                                          real_t target_y) noexcept {
  const real_t dy = y2 - y1;
  if (std::abs(dy) < static_cast<real_t>(1e-12)) {
    return (x1 + x2) * static_cast<real_t>(0.5);
  }
  const real_t frac = (target_y - y1) / dy;
  return x1 + frac * (x2 - x1);
}

} // namespace

real_t PeakIntegrationService::calculateFWHM(std::span<const real_t> x_bins,
                                             std::span<const real_t> y_vals, std::size_t peak_idx,
                                             std::size_t min_idx, std::size_t max_idx) noexcept {
  if (x_bins.empty() || y_vals.empty() || peak_idx >= y_vals.size()) {
    return 0.0;
  }
  const real_t half_max = y_vals[peak_idx] * static_cast<real_t>(0.5);
  if (half_max <= static_cast<real_t>(0.0)) {
    return 0.0;
  }

  real_t left_x = x_bins[min_idx];
  for (std::size_t i = peak_idx; i > min_idx; --i) {
    if (y_vals[i - 1] <= half_max && y_vals[i] >= half_max) {
      left_x = interpolateHalfCross(x_bins[i - 1], y_vals[i - 1], x_bins[i], y_vals[i], half_max);
      break;
    }
  }

  real_t right_x = x_bins[max_idx];
  for (std::size_t i = peak_idx; i < max_idx && (i + 1) < y_vals.size(); ++i) {
    if (y_vals[i] >= half_max && y_vals[i + 1] <= half_max) {
      right_x = interpolateHalfCross(x_bins[i], y_vals[i], x_bins[i + 1], y_vals[i + 1], half_max);
      break;
    }
  }

  return std::max(static_cast<real_t>(0.0), right_x - left_x);
}

PeakIntegralResult PeakIntegrationService::integrateRange(std::span<const real_t> x_bins,
                                                          std::span<const real_t> y_vals,
                                                          real_t r_min, real_t r_max,
                                                          std::optional<real_t> number_density,
                                                          bool is_radial) {
  if (r_min > r_max) {
    std::swap(r_min, r_max);
  }
  const auto num_points = std::min(x_bins.size(), y_vals.size());
  if (num_points < 2) {
    return {};
  }

  const auto [min_idx, max_idx] = findBoundingIndices(x_bins, r_min, r_max);
  if (min_idx >= max_idx) {
    return {};
  }

  const std::size_t peak_idx = locatePeakIndex(y_vals, min_idx, max_idx);
  const real_t peak_x = x_bins[peak_idx];
  const real_t peak_y = y_vals[peak_idx];

  real_t area = 0.0;
  real_t moment = 0.0;
  real_t cn_sum = 0.0;

  for (std::size_t i = min_idx; i < max_idx && (i + 1) < num_points; ++i) {
    const real_t x0 = x_bins[i];
    const real_t x1 = x_bins[i + 1];
    const real_t y0 = y_vals[i];
    const real_t y1 = y_vals[i + 1];
    const real_t dx = x1 - x0;

    area += static_cast<real_t>(0.5) * (y0 + y1) * dx;
    moment += static_cast<real_t>(0.5) * (x0 * y0 + x1 * y1) * dx;

    if (is_radial && number_density.has_value()) {
      const real_t shell0 = x0 * x0 * y0;
      const real_t shell1 = x1 * x1 * y1;
      cn_sum += static_cast<real_t>(0.5) * (shell0 + shell1) * dx;
    }
  }

  const real_t centroid = (area > static_cast<real_t>(1e-12)) ? (moment / area) : peak_x;
  const real_t fwhm = calculateFWHM(x_bins, y_vals, peak_idx, min_idx, max_idx);

  std::optional<real_t> cn;
  if (is_radial && number_density.has_value() && *number_density > static_cast<real_t>(0.0)) {
    cn = FOUR_PI * (*number_density) * cn_sum;
  }

  return PeakIntegralResult{
      .r_peak = peak_x,
      .r_min = x_bins[min_idx],
      .r_max = x_bins[max_idx],
      .peak_height = peak_y,
      .area = area,
      .coordination_number = cn,
      .fwhm = fwhm,
      .centroid = centroid,
  };
}

std::optional<PeakIntegralResult>
PeakIntegrationService::autoDetectPeakNear(std::span<const real_t> x_bins,
                                           std::span<const real_t> y_vals, real_t target_x,
                                           std::optional<real_t> number_density, bool is_radial) {
  const auto num_points = std::min(x_bins.size(), y_vals.size());
  if (num_points < 3) {
    return std::nullopt;
  }

  // 1. Locate closest data point to target_x
  std::size_t closest_idx = 0;
  real_t min_dist = std::numeric_limits<real_t>::max();
  for (std::size_t i = 0; i < num_points; ++i) {
    const real_t d = std::abs(x_bins[i] - target_x);
    if (d < min_dist) {
      min_dist = d;
      closest_idx = i;
    }
  }

  // 2. Search local neighborhood (within 20 bins or until bounds) for local maximum
  const std::size_t search_radius = 25;
  const std::size_t left_bound = (closest_idx > search_radius) ? (closest_idx - search_radius) : 0;
  const std::size_t right_bound = std::min(num_points - 1, closest_idx + search_radius);

  std::size_t peak_idx = closest_idx;
  real_t max_val = y_vals[closest_idx];
  for (std::size_t i = left_bound; i <= right_bound; ++i) {
    if (y_vals[i] > max_val) {
      max_val = y_vals[i];
      peak_idx = i;
    }
  }

  if (max_val <= static_cast<real_t>(1e-5)) {
    return std::nullopt;
  }

  // 3. Find left boundary (onset / minimum)
  const real_t threshold = max_val * MIN_PEAK_THRESHOLD_FRAC;
  std::size_t onset_idx = peak_idx;
  for (std::size_t i = peak_idx; i > 0; --i) {
    if (y_vals[i] <= threshold) {
      onset_idx = i;
      break;
    }
    if (i + 1 < num_points && y_vals[i - 1] > y_vals[i]) {
      onset_idx = i;
      break;
    }
    onset_idx = i - 1;
  }

  // 4. Find right boundary (first minimum / trough)
  std::size_t trough_idx = peak_idx;
  for (std::size_t i = peak_idx; (i + 1) < num_points; ++i) {
    if (y_vals[i] <= threshold) {
      trough_idx = i;
      break;
    }
    if (i > 0 && y_vals[i + 1] > y_vals[i]) {
      trough_idx = i;
      break;
    }
    trough_idx = i + 1;
  }

  if (onset_idx >= trough_idx) {
    return std::nullopt;
  }

  return integrateRange(x_bins, y_vals, x_bins[onset_idx], x_bins[trough_idx], number_density,
                        is_radial);
}

} // namespace correlation::analysis
