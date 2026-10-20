/**
 * @file FirstMinimumDetector.cpp
 * @brief Implementation of FirstMinimumDetector.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/FirstMinimumDetector.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace correlation::analysis {

namespace {

constexpr size_t MIN_DATA_POINTS = 5;
constexpr real_t MIN_PEAK_INTENSITY = 0.05;

[[nodiscard]] std::vector<real_t> smoothSeries(std::span<const real_t> input) {
  std::vector<real_t> output(input.size());
  if (input.size() < 3) {
    std::ranges::copy(input, output.begin());
    return output;
  }

  output.front() = input.front();
  for (size_t i = 1; i + 1 < input.size(); ++i) {
    output[i] = static_cast<real_t>(0.25) * input[i - 1] + static_cast<real_t>(0.5) * input[i] +
                static_cast<real_t>(0.25) * input[i + 1];
  }
  output.back() = input.back();
  return output;
}

[[nodiscard]] std::expected<size_t, DetectionError>
findPrimaryPeak(std::span<const real_t> r_bins, std::span<const real_t> curve,
                const DetectionOptions &opts) noexcept {
  size_t peak_idx = 0;
  real_t max_val = -std::numeric_limits<real_t>::infinity();
  bool found_any = false;

  for (size_t i = 0; i < r_bins.size(); ++i) {
    const real_t r = r_bins[i];
    if (r < opts.min_search_radius) {
      continue;
    }
    if (r > opts.max_search_radius) {
      break;
    }

    if (curve[i] > max_val) {
      max_val = curve[i];
      peak_idx = i;
      found_any = true;
    }
  }

  if (!found_any || max_val < MIN_PEAK_INTENSITY) {
    return std::unexpected(DetectionError::NoPeakFound);
  }
  return peak_idx;
}

[[nodiscard]] real_t findPeakOnset(std::span<const real_t> r_bins, std::span<const real_t> curve,
                                   size_t peak_idx, const DetectionOptions &opts) noexcept {
  const real_t peak_val = curve[peak_idx];
  const real_t threshold = opts.onset_threshold_frac * peak_val;

  real_t r_onset = r_bins[peak_idx];
  for (size_t i = peak_idx; i > 0; --i) {
    if (curve[i] <= threshold) {
      // Linear interpolation between i and i + 1
      const real_t delta_y = curve[i + 1] - curve[i];
      if (std::abs(delta_y) > static_cast<real_t>(1e-6)) {
        const real_t frac = (threshold - curve[i]) / delta_y;
        r_onset = r_bins[i] + frac * (r_bins[i + 1] - r_bins[i]);
      } else {
        r_onset = r_bins[i];
      }
      break;
    }
    r_onset = r_bins[i];
  }

  return std::max(r_onset, opts.core_floor_radius);
}

[[nodiscard]] std::expected<std::pair<real_t, real_t>, DetectionError>
findFirstTrough(std::span<const real_t> r_bins, std::span<const real_t> curve, size_t peak_idx,
                const DetectionOptions &opts) noexcept {
  const real_t r_peak = r_bins[peak_idx];
  const real_t r_limit = r_peak + opts.search_window;

  size_t trough_idx = peak_idx;
  real_t min_val = curve[peak_idx];
  std::optional<size_t> zero_crossing_idx;

  real_t prev_deriv = 0.0;
  bool prev_deriv_valid = false;

  for (size_t i = peak_idx + 1; i + 1 < r_bins.size(); ++i) {
    if (r_bins[i] > r_limit) {
      break;
    }

    if (curve[i] < min_val) {
      min_val = curve[i];
      trough_idx = i;
    }

    const real_t delta_r = r_bins[i + 1] - r_bins[i - 1];
    if (delta_r > static_cast<real_t>(1e-6)) {
      const real_t deriv = (curve[i + 1] - curve[i - 1]) / delta_r;
      if (prev_deriv_valid && prev_deriv <= static_cast<real_t>(0.0) &&
          deriv > static_cast<real_t>(0.0)) {
        zero_crossing_idx = i;
        break;
      }
      prev_deriv = deriv;
      prev_deriv_valid = true;
    }
  }

  if (zero_crossing_idx.has_value()) {
    const size_t z_idx = *zero_crossing_idx;
    return std::make_pair(r_bins[z_idx], curve[z_idx]);
  }

  if (trough_idx > peak_idx && min_val < curve[peak_idx]) {
    return std::make_pair(r_bins[trough_idx], min_val);
  }

  return std::unexpected(DetectionError::NoMinimumFound);
}

} // namespace

std::expected<CutoffResult, DetectionError>
FirstMinimumDetector::detectCutoffs(std::span<const real_t> r_bins, std::span<const real_t> g_r,
                                    const DetectionOptions &opts) {
  if (r_bins.size() != g_r.size()) {
    return std::unexpected(DetectionError::MismatchedSizes);
  }
  if (r_bins.size() < MIN_DATA_POINTS || g_r.size() < MIN_DATA_POINTS) {
    return std::unexpected(DetectionError::InsufficientData);
  }
  if (opts.min_search_radius >= opts.max_search_radius || opts.search_window <= 0.0) {
    return std::unexpected(DetectionError::InvalidRange);
  }

  std::vector<real_t> smoothed_storage;
  std::span<const real_t> active_curve = g_r;
  if (opts.apply_smoothing) {
    smoothed_storage = smoothSeries(g_r);
    active_curve = smoothed_storage;
  }

  const auto peak_res = findPrimaryPeak(r_bins, active_curve, opts);
  if (!peak_res.has_value()) {
    return std::unexpected(peak_res.error());
  }

  const size_t peak_idx = *peak_res;
  const real_t r_peak = r_bins[peak_idx];
  const real_t peak_height = active_curve[peak_idx];

  const real_t r_min = findPeakOnset(r_bins, active_curve, peak_idx, opts);

  const auto trough_res = findFirstTrough(r_bins, active_curve, peak_idx, opts);
  if (!trough_res.has_value()) {
    return std::unexpected(trough_res.error());
  }

  const auto [r_max, min_height] = *trough_res;

  return CutoffResult{
      .r_min = r_min,
      .r_max = r_max,
      .r_peak = r_peak,
      .peak_height = peak_height,
      .min_height = min_height,
  };
}

std::expected<CutoffResult, DetectionError>
FirstMinimumDetector::detectCutoffs(const Histogram &hist, std::string_view partial_key,
                                    const DetectionOptions &opts) {
  const std::string key(partial_key);
  const auto iter_partial = hist.partials.find(key);
  if (iter_partial == hist.partials.end()) {
    // Check if total requested and available under alternative casing
    if (key == "Total" && !hist.partials.empty()) {
      return detectCutoffs(hist.bins, hist.partials.begin()->second, opts);
    }
    return std::unexpected(DetectionError::PartialNotFound);
  }

  // If pre-smoothed partial exists and smoothing is disabled or requested, use smoothed directly
  if (const auto iter_smoothed = hist.smoothed_partials.find(key);
      iter_smoothed != hist.smoothed_partials.end() && !iter_smoothed->second.empty()) {
    DetectionOptions direct_opts = opts;
    direct_opts.apply_smoothing = false;
    return detectCutoffs(hist.bins, iter_smoothed->second, direct_opts);
  }

  return detectCutoffs(hist.bins, iter_partial->second, opts);
}

std::string_view FirstMinimumDetector::formatError(DetectionError error) noexcept {
  switch (error) {
  case DetectionError::InsufficientData:
    return "Insufficient data points to perform derivative analysis.";
  case DetectionError::MismatchedSizes:
    return "Bin array and curve intensity array dimensions do not match.";
  case DetectionError::NoPeakFound:
    return "No primary coordination peak detected within radial bounds.";
  case DetectionError::NoMinimumFound:
    return "No coordination trough or derivative zero-crossing detected after primary peak.";
  case DetectionError::InvalidRange:
    return "Invalid radial search range or search window specified.";
  case DetectionError::PartialNotFound:
    return "Target element partial curve was not found in dataset.";
  }
  return "Unknown detection error.";
}

} // namespace correlation::analysis
