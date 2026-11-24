/**
 * @file LDOSPostProcessor.cpp
 * @brief Implementation of electronic LDOS/TDOS post-processing utilities.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/scattering/LDOSPostProcessor.hpp"
#include "math/Constants.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace correlation::calculators {

namespace {

constexpr double K_BOLTZMANN_EV_PER_K = 8.617333262e-5;
constexpr double K_MIN_POSITIVE = 1e-12;

[[nodiscard]] double computeTrapezoidalIntegral(std::span<const real_t> energies,
                                                std::span<const double> values) noexcept {
  if (energies.size() < 2 || energies.size() != values.size()) {
    return 0.0;
  }
  double area = 0.0;
  for (size_t i = 0; i + 1 < energies.size(); ++i) {
    const double d_e = static_cast<double>(energies[i + 1]) - static_cast<double>(energies[i]);
    area += 0.5 * (values[i] + values[i + 1]) * d_e;
  }
  return area;
}

[[nodiscard]] double computeBinWidth(std::span<const real_t> energies, size_t idx) noexcept {
  const size_t count = energies.size();
  if (count < 2) {
    return 1.0;
  }
  if (idx == 0) {
    return static_cast<double>(energies[1]) - static_cast<double>(energies[0]);
  }
  if (idx + 1 == count) {
    return static_cast<double>(energies[count - 1]) - static_cast<double>(energies[count - 2]);
  }
  return 0.5 * (static_cast<double>(energies[idx + 1]) - static_cast<double>(energies[idx - 1]));
}

[[nodiscard]] double
evaluateKernelWeight(double delta_e, double sigma,
                     correlation::analysis::LDOSBroadeningKernel kernel) noexcept {
  if (kernel == correlation::analysis::LDOSBroadeningKernel::Gaussian) {
    const double ratio = delta_e / sigma;
    return std::exp(-0.5 * ratio * ratio) / (sigma * std::sqrt(correlation::math::TWO_PI));
  }
  const double denom = delta_e * delta_e + sigma * sigma;
  return (sigma / correlation::math::PI) / denom;
}

[[nodiscard]] double interpolateLinearAt(std::span<const real_t> energies,
                                         std::span<const real_t> values, double target_e) noexcept {
  const size_t count = energies.size();
  if (count == 0 || values.size() != count) {
    return 0.0;
  }
  if (target_e <= static_cast<double>(energies.front())) {
    return static_cast<double>(values.front());
  }
  if (target_e >= static_cast<double>(energies.back())) {
    return static_cast<double>(values.back());
  }
  const auto iter = std::ranges::lower_bound(energies, static_cast<real_t>(target_e));
  const auto idx = static_cast<size_t>(std::distance(energies.begin(), iter));
  if (idx == 0) {
    return static_cast<double>(values[0]);
  }
  const auto energy_lo = static_cast<double>(energies[idx - 1]);
  const auto energy_hi = static_cast<double>(energies[idx]);
  const double span = energy_hi - energy_lo;
  if (span <= K_MIN_POSITIVE) {
    return static_cast<double>(values[idx]);
  }
  const double frac = (target_e - energy_lo) / span;
  return static_cast<double>(values[idx - 1]) * (1.0 - frac) +
         static_cast<double>(values[idx]) * frac;
}

[[nodiscard]] double evaluateFermiDiracDerivativeWeight(double energy_ev, double kbt_ev) noexcept {
  const double arg = 0.5 * energy_ev / kbt_ev;
  if (std::abs(arg) > 35.0) {
    return 0.0;
  }
  const double cosh_val = std::cosh(arg);
  return 1.0 / (4.0 * kbt_ev * cosh_val * cosh_val);
}

} // namespace

std::vector<real_t>
LDOSPostProcessor::broadenSpectrum(std::span<const real_t> energies, std::span<const real_t> dos,
                                   correlation::analysis::LDOSBroadeningKernel kernel,
                                   real_t sigma) {
  const size_t count = energies.size();
  if (count == 0 || dos.size() != count) {
    return {};
  }
  const auto sigma_d = static_cast<double>(sigma);
  if (kernel == correlation::analysis::LDOSBroadeningKernel::None || sigma_d <= K_MIN_POSITIVE) {
    return {dos.begin(), dos.end()};
  }

  std::vector<double> raw_d(count);
  std::vector<double> widths(count);
  for (size_t i = 0; i < count; ++i) {
    raw_d[i] = static_cast<double>(dos[i]);
    widths[i] = std::max(0.0, computeBinWidth(energies, i));
  }
  const double orig_area = computeTrapezoidalIntegral(energies, raw_d);

  std::vector<double> convolved(count, 0.0);
  for (size_t i = 0; i < count; ++i) {
    const auto e_i = static_cast<double>(energies[i]);
    double acc = 0.0;
    double norm = 0.0;
    for (size_t j = 0; j < count; ++j) {
      const double weight =
          evaluateKernelWeight(e_i - static_cast<double>(energies[j]), sigma_d, kernel) * widths[j];
      acc += raw_d[j] * weight;
      norm += weight;
    }
    convolved[i] = (norm > K_MIN_POSITIVE) ? (acc / norm) : raw_d[i];
  }

  const double conv_area = computeTrapezoidalIntegral(energies, convolved);
  const double scale =
      (orig_area > K_MIN_POSITIVE && conv_area > K_MIN_POSITIVE) ? (orig_area / conv_area) : 1.0;

  std::vector<real_t> result(count);
  for (size_t i = 0; i < count; ++i) {
    result[i] = static_cast<real_t>(std::max(0.0, convolved[i] * scale));
  }
  return result;
}

std::vector<real_t> LDOSPostProcessor::rebinSpectrum(std::span<const real_t> src_energies,
                                                     std::span<const real_t> src_dos,
                                                     std::span<const real_t> dst_energies) {
  if (src_energies.empty() || src_dos.size() != src_energies.size() || dst_energies.empty()) {
    return {};
  }
  std::vector<real_t> result(dst_energies.size(), static_cast<real_t>(0.0));
  const auto src_min = static_cast<double>(src_energies.front());
  const auto src_max = static_cast<double>(src_energies.back());
  for (size_t i = 0; i < dst_energies.size(); ++i) {
    const auto e_val = static_cast<double>(dst_energies[i]);
    if (e_val < src_min || e_val > src_max) {
      result[i] = static_cast<real_t>(0.0);
    } else {
      result[i] =
          static_cast<real_t>(std::max(0.0, interpolateLinearAt(src_energies, src_dos, e_val)));
    }
  }
  return result;
}

real_t LDOSPostProcessor::evaluateFermiDos(std::span<const real_t> energies,
                                           std::span<const real_t> dos,
                                           real_t electronic_temp_k) noexcept {
  const size_t count = energies.size();
  if (count == 0 || dos.size() != count) {
    return static_cast<real_t>(0.0);
  }
  if (static_cast<double>(energies.front()) > 0.0 || static_cast<double>(energies.back()) < 0.0) {
    return static_cast<real_t>(0.0);
  }

  const double kbt_ev = static_cast<double>(electronic_temp_k) * K_BOLTZMANN_EV_PER_K;
  if (kbt_ev <= 1e-6 || count < 2) {
    return static_cast<real_t>(std::max(0.0, interpolateLinearAt(energies, dos, 0.0)));
  }

  double weighted_sum = 0.0;
  double weight_norm = 0.0;
  for (size_t i = 0; i < count; ++i) {
    const double d_e = std::max(0.0, computeBinWidth(energies, i));
    const double weight =
        evaluateFermiDiracDerivativeWeight(static_cast<double>(energies[i]), kbt_ev) * d_e;
    weighted_sum += static_cast<double>(dos[i]) * weight;
    weight_norm += weight;
  }

  if (weight_norm <= K_MIN_POSITIVE) {
    return static_cast<real_t>(std::max(0.0, interpolateLinearAt(energies, dos, 0.0)));
  }
  return static_cast<real_t>(std::max(0.0, weighted_sum / weight_norm));
}

BandGapInfo LDOSPostProcessor::detectBandGap(std::span<const real_t> energies,
                                             std::span<const real_t> dos,
                                             real_t threshold) noexcept {
  BandGapInfo info{};
  const size_t count = energies.size();
  if (count < 2 || dos.size() != count) {
    return info;
  }

  const double thresh = std::max(0.0, static_cast<double>(threshold));
  const double dos_at_ef = interpolateLinearAt(energies, dos, 0.0);
  if (dos_at_ef > thresh) {
    info.is_metallic = true;
    return info;
  }

  bool found_vbm = false;
  bool found_cbm = false;
  double vbm = 0.0;
  double cbm = 0.0;

  for (size_t i = 0; i < count; ++i) {
    const auto e_val = static_cast<double>(energies[i]);
    const auto d_val = static_cast<double>(dos[i]);
    if (d_val <= thresh) {
      continue;
    }
    if (e_val <= 0.0) {
      vbm = e_val;
      found_vbm = true;
    } else if (!found_cbm) {
      cbm = e_val;
      found_cbm = true;
    }
  }

  if (found_vbm && found_cbm && cbm > vbm) {
    info.vbm_energy = static_cast<real_t>(vbm);
    info.cbm_energy = static_cast<real_t>(cbm);
    info.band_gap = static_cast<real_t>(cbm - vbm);
  }
  return info;
}

std::vector<real_t> LDOSPostProcessor::computeIntegratedDos(std::span<const real_t> energies,
                                                            std::span<const real_t> dos) {
  const size_t count = energies.size();
  if (count == 0 || dos.size() != count) {
    return {};
  }
  std::vector<real_t> idos(count, static_cast<real_t>(0.0));
  double cumulative = 0.0;
  for (size_t i = 1; i < count; ++i) {
    const double d_e =
        std::max(0.0, static_cast<double>(energies[i]) - static_cast<double>(energies[i - 1]));
    const double avg_dos = 0.5 * (std::max(0.0, static_cast<double>(dos[i - 1])) +
                                  std::max(0.0, static_cast<double>(dos[i])));
    cumulative += avg_dos * d_e;
    idos[i] = static_cast<real_t>(cumulative);
  }
  return idos;
}

} // namespace correlation::calculators
