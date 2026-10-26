/**
 * @file DistributionPostProcessor.cpp
 * @brief Implementation of histogram smoothing and post-processing service.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/DistributionPostProcessor.hpp"
#include "analysis/DistributionFunctions.hpp"

#include <algorithm>
#include <stdexcept>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <vector>

namespace correlation::analysis {

namespace {

[[nodiscard]] constexpr bool isAngularDistribution(std::string_view name) noexcept {
  return name == "BAD" || name == "PAD" || name == "PAD_raw" || name == "DAD" || name == "DAD_raw";
}

[[nodiscard]] real_t computeEffectiveSigma(std::string_view name, real_t sigma,
                                           real_t bin_dx) noexcept {
  real_t effective_sigma = sigma;
  if (isAngularDistribution(name)) {
    constexpr auto ANGULAR_DOMAIN_SCALE = static_cast<real_t>(15.0);
    effective_sigma = sigma * ANGULAR_DOMAIN_SCALE;
  }
  return std::max(bin_dx, effective_sigma);
}

} // namespace

void DistributionPostProcessor::smooth(Histogram &hist, std::string_view name, real_t sigma,
                                       correlation::math::KernelType kernel) {
  hist.smoothed_partials.clear();

  if (hist.bins.size() < 2) {
    return;
  }

  const real_t bin_dx = hist.bins[1] - hist.bins[0];
  if (bin_dx <= static_cast<real_t>(0.0)) {
    return;
  }

  const real_t min_sigma = computeEffectiveSigma(name, sigma, bin_dx);

  using PartialEntry = std::pair<const std::string *, const std::vector<real_t> *>;
  std::vector<PartialEntry> entries;
  entries.reserve(hist.partials.size());
  for (const auto &[key, vals] : hist.partials) {
    entries.emplace_back(&key, &vals);
  }

  std::vector<std::vector<real_t>> results(entries.size());
  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, entries.size()), [&](const tbb::blocked_range<size_t> &range) {
        for (size_t i = range.begin(); i != range.end(); ++i) {
          results[i] =
              correlation::math::KernelSmoothing(*entries[i].second, {
                                                                         .bin_width = bin_dx,
                                                                         .sigma = min_sigma,
                                                                         .type = kernel,
                                                                     });
        }
      });

  for (size_t i = 0; i < entries.size(); ++i) {
    hist.smoothed_partials[*entries[i].first] = std::move(results[i]);
  }
}

void DistributionPostProcessor::smoothAll(std::map<std::string, Histogram> &histograms,
                                          real_t sigma, correlation::math::KernelType kernel) {
  for (auto &[name, histogram] : histograms) {
    smooth(histogram, name, sigma, kernel);
  }
}

void DistributionPostProcessor::smooth(DistributionFunctions &dists, const std::string &name,
                                       real_t sigma, correlation::math::KernelType kernel) {
  if (!dists.hasHistogram(name)) {
    throw std::runtime_error("Histogram '" + name + "' not found for smoothing.");
  }
  smooth(dists.getMutableHistogram(name), name, sigma, kernel);
}

void DistributionPostProcessor::smoothAll(DistributionFunctions &dists, real_t sigma,
                                          correlation::math::KernelType kernel) {
  for (auto &[name, histogram] : dists.getAllHistograms()) {
    smooth(histogram, name, sigma, kernel);
  }
}

} // namespace correlation::analysis
