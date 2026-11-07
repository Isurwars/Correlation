/**
 * @file ChargeCorrelationCalculator.cpp
 * @brief Implementation of charge-weighted pair distribution g_qq(r) and screening calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/ml/ChargeCorrelationCalculator.hpp"
#include "math/Constants.hpp"

#include <cmath>
#include <stdexcept>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace correlation::calculators {

namespace {

struct alignas(64) LocalChargeHistograms {
  std::vector<real_t> charge_prod;
  std::vector<uint64_t> counts;

  explicit LocalChargeHistograms(size_t num_bins)
      : charge_prod(num_bins, static_cast<real_t>(0.0)), counts(num_bins, 0) {}
};

[[nodiscard]] size_t computeBinIndex(real_t distance, real_t r_bin_width,
                                     size_t num_bins) noexcept {
  if (distance <= static_cast<real_t>(0.0) || r_bin_width <= static_cast<real_t>(0.0)) {
    return num_bins;
  }
  const auto bin = static_cast<size_t>(distance / r_bin_width);
  return (bin < num_bins) ? bin : num_bins;
}

void processChargePair(const correlation::math::Vector3<real_t> &pos_i,
                       const correlation::math::Vector3<real_t> &pos_j, real_t q_i, real_t q_j,
                       const correlation::core::Cell &cell, bool is_periodic, real_t r_max,
                       real_t r_bin_width, size_t num_bins,
                       LocalChargeHistograms &local_hist) noexcept {
  correlation::math::Vector3<real_t> delta = pos_j - pos_i;
  if (is_periodic) {
    delta = cell.minimumImage(delta);
  }

  const real_t dist = correlation::math::norm(delta);
  if (dist <= static_cast<real_t>(0.0) || dist >= r_max) {
    return;
  }

  const size_t bin = computeBinIndex(dist, r_bin_width, num_bins);
  if (bin >= num_bins) {
    return;
  }

  local_hist.charge_prod[bin] += (q_i * q_j);
  local_hist.counts[bin] += 1;
}

void accumulateAllChargePairs(const correlation::core::Cell &cell, std::span<const real_t> charges,
                              real_t r_max, real_t r_bin_width, size_t num_bins,
                              ChargeCorrelationResult &result) {
  const size_t n_atoms = cell.atoms().size();
  const bool is_periodic = (cell.volume() > static_cast<real_t>(1e-6));

  tbb::enumerable_thread_specific<LocalChargeHistograms> thread_histograms(
      [num_bins]() { return LocalChargeHistograms(num_bins); });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, n_atoms),
                    [&](const tbb::blocked_range<size_t> &range) {
                      auto &local_hist = thread_histograms.local();
                      for (size_t i = range.begin(); i < range.end(); ++i) {
                        const auto &pos_i = cell.atoms()[i].position();
                        const real_t q_i = charges[i];
                        for (size_t j = i + 1; j < n_atoms; ++j) {
                          const auto &pos_j = cell.atoms()[j].position();
                          const real_t q_j = charges[j];
                          processChargePair(pos_i, pos_j, q_i, q_j, cell, is_periodic, r_max,
                                            r_bin_width, num_bins, local_hist);
                        }
                      }
                    });

  for (const auto &local_hist : thread_histograms) {
    for (size_t b = 0; b < num_bins; ++b) {
      result.g_qq[b] += local_hist.charge_prod[b];
      result.pair_counts[b] += local_hist.counts[b];
    }
  }
}

real_t computeSumSquaredCharges(std::span<const real_t> charges) noexcept {
  real_t sum_q2 = static_cast<real_t>(0.0);
  for (const real_t q : charges) {
    sum_q2 += (q * q);
  }
  return sum_q2;
}

void finalizeChargeNormalization(const correlation::core::Cell &cell,
                                 std::span<const real_t> charges, size_t num_bins,
                                 real_t r_bin_width, ChargeCorrelationResult &result) noexcept {
  const real_t sum_q2 = computeSumSquaredCharges(charges);
  real_t volume = cell.volume();
  if (volume <= static_cast<real_t>(1e-6)) {
    volume = static_cast<real_t>(1000.0);
  }

  real_t running_integral = static_cast<real_t>(0.0);
  constexpr real_t four_pi = correlation::math::FOUR_PI;

  for (size_t b = 0; b < num_bins; ++b) {
    const real_t r = (static_cast<real_t>(b) + static_cast<real_t>(0.5)) * r_bin_width;
    result.r_bins[b] = r;

    // Shell volume: dV = 4 * pi * r^2 * dr
    const real_t dV = four_pi * r * r * r_bin_width;

    if (dV > static_cast<real_t>(0.0) && sum_q2 > static_cast<real_t>(1e-9)) {
      // 2.0 factor because pair sum runs over i < j
      const real_t raw_prod_sum = static_cast<real_t>(2.0) * result.g_qq[b];
      result.g_qq[b] = (volume * raw_prod_sum) / (dV * sum_q2);
      running_integral += raw_prod_sum;
      result.charge_integral[b] = running_integral;
    } else {
      result.g_qq[b] = static_cast<real_t>(0.0);
      result.charge_integral[b] = static_cast<real_t>(0.0);
    }
  }
}

} // namespace

void ChargeCorrelationCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions & /*dists*/,
    const correlation::analysis::AnalysisSettings & /*settings*/) const {}

ChargeCorrelationResult ChargeCorrelationCalculator::calculate(const correlation::core::Cell &cell,
                                                               std::span<const real_t> charges,
                                                               real_t r_max, real_t r_bin_width) {
  if (r_max <= static_cast<real_t>(0.0) || r_bin_width <= static_cast<real_t>(0.0)) {
    throw std::invalid_argument(
        "[ChargeCorrelationCalculator] Cutoff radius and bin width must be positive.");
  }

  const size_t n_atoms = cell.atoms().size();
  if (charges.size() < n_atoms) {
    throw std::invalid_argument(
        "[ChargeCorrelationCalculator] Mismatch between atom count and charge span.");
  }

  const auto num_bins = static_cast<size_t>(std::ceil(r_max / r_bin_width));
  ChargeCorrelationResult result{
      .r_bins = std::vector<real_t>(num_bins, static_cast<real_t>(0.0)),
      .g_qq = std::vector<real_t>(num_bins, static_cast<real_t>(0.0)),
      .charge_integral = std::vector<real_t>(num_bins, static_cast<real_t>(0.0)),
      .pair_counts = std::vector<uint64_t>(num_bins, 0),
  };

  if (n_atoms < 2) {
    finalizeChargeNormalization(cell, charges, num_bins, r_bin_width, result);
    return result;
  }

  accumulateAllChargePairs(cell, charges, r_max, r_bin_width, num_bins, result);
  finalizeChargeNormalization(cell, charges, num_bins, r_bin_width, result);
  return result;
}

ChargeCorrelationResult
ChargeCorrelationCalculator::calculateFromMLIP(const correlation::core::Cell &cell,
                                               const correlation::mlip::MLIPOutput &mlip_out,
                                               real_t r_max, real_t r_bin_width) {
  return calculate(cell, mlip_out.charges, r_max, r_bin_width);
}

} // namespace correlation::calculators
