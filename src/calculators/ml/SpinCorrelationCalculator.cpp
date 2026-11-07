/**
 * @file SpinCorrelationCalculator.cpp
 * @brief Implementation of spatial spin-spin and magnetic correlation calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/ml/SpinCorrelationCalculator.hpp"
#include "math/Constants.hpp"

#include <cmath>
#include <stdexcept>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace correlation::calculators {

namespace {

struct alignas(64) LocalSpinHistograms {
  std::vector<real_t> spin_dot;
  std::vector<real_t> mag_dot;
  std::vector<uint64_t> counts;

  explicit LocalSpinHistograms(size_t num_bins)
      : spin_dot(num_bins, static_cast<real_t>(0.0)), mag_dot(num_bins, static_cast<real_t>(0.0)),
        counts(num_bins, 0) {}
};

[[nodiscard]] size_t computeBinIndex(real_t distance, real_t r_bin_width,
                                     size_t num_bins) noexcept {
  if (distance <= static_cast<real_t>(0.0) || r_bin_width <= static_cast<real_t>(0.0)) {
    return num_bins;
  }
  const auto bin = static_cast<size_t>(distance / r_bin_width);
  return (bin < num_bins) ? bin : num_bins;
}

void processAtomPair(const correlation::math::Vector3<real_t> &pos_i,
                     const correlation::math::Vector3<real_t> &pos_j,
                     const correlation::math::Vector3<real_t> &m_i,
                     const correlation::math::Vector3<real_t> &m_j,
                     const correlation::core::Cell &cell, bool is_periodic, real_t r_max,
                     real_t r_bin_width, size_t num_bins,
                     LocalSpinHistograms &local_hist) noexcept {
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

  const real_t mag_i_len = correlation::math::norm(m_i);
  const real_t mag_j_len = correlation::math::norm(m_j);
  const real_t dot_m = m_i.x() * m_j.x() + m_i.y() * m_j.y() + m_i.z() * m_j.z();

  real_t dot_s = static_cast<real_t>(0.0);
  if (mag_i_len > static_cast<real_t>(1e-6) && mag_j_len > static_cast<real_t>(1e-6)) {
    dot_s = dot_m / (mag_i_len * mag_j_len);
    if (dot_s > static_cast<real_t>(1.0))
      dot_s = static_cast<real_t>(1.0);
    if (dot_s < static_cast<real_t>(-1.0))
      dot_s = static_cast<real_t>(-1.0);
  }

  local_hist.spin_dot[bin] += dot_s;
  local_hist.mag_dot[bin] += dot_m;
  local_hist.counts[bin] += 1;
}

void accumulateAllPairs(const correlation::core::Cell &cell,
                        std::span<const correlation::math::Vector3<real_t>> magmoms, real_t r_max,
                        real_t r_bin_width, size_t num_bins, SpinCorrelationResult &result) {
  const size_t n_atoms = cell.atoms().size();
  const bool is_periodic = (cell.volume() > static_cast<real_t>(1e-6));

  tbb::enumerable_thread_specific<LocalSpinHistograms> thread_histograms(
      [num_bins]() { return LocalSpinHistograms(num_bins); });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, n_atoms),
                    [&](const tbb::blocked_range<size_t> &range) {
                      auto &local_hist = thread_histograms.local();
                      for (size_t i = range.begin(); i < range.end(); ++i) {
                        const auto &pos_i = cell.atoms()[i].position();
                        const auto &m_i = magmoms[i];
                        for (size_t j = i + 1; j < n_atoms; ++j) {
                          const auto &pos_j = cell.atoms()[j].position();
                          const auto &m_j = magmoms[j];
                          processAtomPair(pos_i, pos_j, m_i, m_j, cell, is_periodic, r_max,
                                          r_bin_width, num_bins, local_hist);
                        }
                      }
                    });

  for (const auto &local_hist : thread_histograms) {
    for (size_t b = 0; b < num_bins; ++b) {
      result.spin_correlation[b] += local_hist.spin_dot[b];
      result.magnetic_moment_corr[b] += local_hist.mag_dot[b];
      result.pair_counts[b] += local_hist.counts[b];
    }
  }
}

void finalizeAverages(size_t num_bins, real_t r_bin_width, SpinCorrelationResult &result) noexcept {
  for (size_t b = 0; b < num_bins; ++b) {
    result.r_bins[b] = (static_cast<real_t>(b) + static_cast<real_t>(0.5)) * r_bin_width;
    if (result.pair_counts[b] > 0) {
      result.spin_correlation[b] /= static_cast<real_t>(result.pair_counts[b]);
      result.magnetic_moment_corr[b] /= static_cast<real_t>(result.pair_counts[b]);
    } else {
      result.spin_correlation[b] = static_cast<real_t>(0.0);
      result.magnetic_moment_corr[b] = static_cast<real_t>(0.0);
    }
  }
}

} // namespace

void SpinCorrelationCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions & /*dists*/,
    const correlation::analysis::AnalysisSettings & /*settings*/) const {}

SpinCorrelationResult
SpinCorrelationCalculator::calculate(const correlation::core::Cell &cell,
                                     std::span<const correlation::math::Vector3<real_t>> magmoms,
                                     real_t r_max, real_t r_bin_width) {
  if (r_max <= static_cast<real_t>(0.0) || r_bin_width <= static_cast<real_t>(0.0)) {
    throw std::invalid_argument(
        "[SpinCorrelationCalculator] Cutoff radius and bin width must be positive.");
  }

  const size_t n_atoms = cell.atoms().size();
  if (magmoms.size() < n_atoms) {
    throw std::invalid_argument(
        "[SpinCorrelationCalculator] Mismatch between atom count and magnetic moment span.");
  }

  const auto num_bins = static_cast<size_t>(std::ceil(r_max / r_bin_width));
  SpinCorrelationResult result{
      .r_bins = std::vector<real_t>(num_bins, static_cast<real_t>(0.0)),
      .spin_correlation = std::vector<real_t>(num_bins, static_cast<real_t>(0.0)),
      .magnetic_moment_corr = std::vector<real_t>(num_bins, static_cast<real_t>(0.0)),
      .pair_counts = std::vector<uint64_t>(num_bins, 0),
  };

  if (n_atoms < 2) {
    finalizeAverages(num_bins, r_bin_width, result);
    return result;
  }

  accumulateAllPairs(cell, magmoms, r_max, r_bin_width, num_bins, result);
  finalizeAverages(num_bins, r_bin_width, result);
  return result;
}

SpinCorrelationResult
SpinCorrelationCalculator::calculateFromMLIP(const correlation::core::Cell &cell,
                                             const correlation::mlip::MLIPOutput &mlip_out,
                                             real_t r_max, real_t r_bin_width) {
  return calculate(cell, mlip_out.magmoms, r_max, r_bin_width);
}

} // namespace correlation::calculators
