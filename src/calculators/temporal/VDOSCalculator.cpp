/**
 * @file VDOSCalculator.cpp
 * @brief Implementation of vibrational density of states calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/temporal/VDOSCalculator.hpp"
#include "analysis/DynamicsAnalyzer.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "core/Cell.hpp"
#include "math/Constants.hpp"
#include "physics/PhysicalData.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
const bool REGISTERED = CalculatorFactory::registerTypeSafe<VDOSCalculator>("VDOSCalculator");

[[nodiscard]] double computeTrapezoidalArea(const std::vector<real_t> &freqs,
                                            const std::vector<real_t> &vals) noexcept {
  if (freqs.size() < 2 || freqs.size() != vals.size()) {
    return 0.0;
  }
  double area = 0.0;
  for (size_t i = 0; i + 1 < freqs.size(); ++i) {
    const double d_nu = static_cast<double>(freqs[i + 1]) - static_cast<double>(freqs[i]);
    area += 0.5 * (static_cast<double>(vals[i]) + static_cast<double>(vals[i + 1])) * d_nu;
  }
  return area;
}

[[nodiscard]] std::map<std::string, real_t>
computeElementMassWeights(const correlation::core::Cell *cell,
                          const correlation::analysis::Histogram &vacf_hist) {
  std::map<std::string, real_t> weights;
  if (cell == nullptr || cell->atoms().empty()) {
    return weights;
  }
  std::map<std::string, size_t> counts;
  for (const auto &atom : cell->atoms()) {
    counts[atom.element().symbol]++;
  }
  double total_mass = 0.0;
  std::map<std::string, double> raw_mass;
  for (const auto &[sym, cnt] : counts) {
    if (!vacf_hist.partials.contains(sym)) {
      continue;
    }
    double mass = 1.0;
    if (const auto m_opt = correlation::physics::findAtomicMass(sym);
        m_opt.has_value() && *m_opt > static_cast<real_t>(0.0)) {
      mass = static_cast<double>(*m_opt);
    }
    const double species_mass = static_cast<double>(cnt) * mass;
    raw_mass[sym] = species_mass;
    total_mass += species_mass;
  }
  if (total_mass > 0.0) {
    for (const auto &[sym, species_mass] : raw_mass) {
      weights[sym] = static_cast<real_t>(species_mass / total_mass);
    }
  }
  return weights;
}

void normalizeVdosSpectra(const std::vector<real_t> &freqs, std::vector<real_t> &total,
                          std::map<std::string, std::vector<real_t>> &partials) noexcept {
  const double total_area = computeTrapezoidalArea(freqs, total);
  if (total_area <= 0.0) {
    return;
  }
  const auto inv_area = static_cast<real_t>(1.0 / total_area);
  for (real_t &val : total) {
    val *= inv_area;
  }
  for (auto &[sym, partial_dos] : partials) {
    for (real_t &val : partial_dos) {
      val *= inv_area;
    }
  }
}

/// Solves the 2PT hard-sphere fluidicity equation F(f) = 0 on f ∈ [0, 1] via bisection.
[[nodiscard]] double solveTwoPhaseFluidicity(double delta) noexcept {
  if (delta <= 1e-12) {
    return 0.0;
  }
  double low = 0.0;
  double high = 1.0;
  const double inv_delta_32 = std::pow(delta, -1.5);
  for (int iter = 0; iter < 60; ++iter) {
    const double mid = 0.5 * (low + high);
    const double z_val = std::min(1e6, inv_delta_32 * std::pow(mid, 2.5));
    const double f_val =
        2.0 * z_val * z_val * z_val - 6.0 * z_val * z_val + (6.0 - mid) * z_val + 2.0 * mid - 2.0;
    if (f_val > 0.0) {
      high = mid;
    } else {
      low = mid;
    }
  }
  return 0.5 * (low + high);
}

/// Estimates 2PT fluidicity f ∈ [0, 1] from (cell, C(0), g(0)) or spectral half-width fallback.
[[nodiscard]] double computeFluidicity(double g0_omega, double vacf_zero,
                                       const std::vector<real_t> &full_omega,
                                       const std::vector<real_t> &full_total,
                                       const correlation::core::Cell *cell) noexcept {
  if (g0_omega <= 1e-12) {
    return 0.0;
  }
  if (cell != nullptr && cell->volume() > static_cast<real_t>(1e-9) && !cell->atoms().empty() &&
      vacf_zero > 0.0) {
    const double rho =
        static_cast<double>(cell->atoms().size()) / static_cast<double>(cell->volume());
    const double v_therm = std::sqrt(correlation::math::PI * vacf_zero / 3.0) * 1000.0;
    const double hs_factor = std::cbrt(rho) * std::pow(6.0 / correlation::math::PI, 2.0 / 3.0);
    const double delta = (4.0 * correlation::math::PI / 3.0) * g0_omega * v_therm * hs_factor;
    return solveTwoPhaseFluidicity(delta);
  }
  // Fallback without a Cell: match Lorentzian half-width at half-maximum g(ω_1/2) = g(0)/2.
  const double half_g0 = 0.5 * g0_omega;
  for (size_t i = 1; i < full_omega.size(); ++i) {
    if (static_cast<double>(full_total[i]) <= half_g0) {
      const auto omega_half = static_cast<double>(full_omega[i]);
      return std::clamp(0.5 * correlation::math::PI * g0_omega * omega_half, 1e-6, 1.0);
    }
  }
  return 1.0;
}

/// Splits a normalized spectrum into Lorentzian translational and vibrational-only components.
void splitSingleSpectrum(const std::vector<real_t> &omega, const std::vector<real_t> &total,
                         double fluidicity, std::vector<real_t> &vib_out,
                         std::vector<real_t> &trans_out) {
  const size_t num_pts = omega.size();
  vib_out.assign(num_pts, static_cast<real_t>(0.0));
  trans_out.assign(num_pts, static_cast<real_t>(0.0));
  if (num_pts == 0) {
    return;
  }
  const auto g0_val = static_cast<double>(total[0]);
  if (fluidicity <= 1e-12 || g0_val <= 1e-12) {
    vib_out = total;
    return;
  }
  const double coeff = (correlation::math::PI * g0_val) / (2.0 * fluidicity);
  for (size_t i = 0; i < num_pts; ++i) {
    const double arg = coeff * static_cast<double>(omega[i]);
    const double g_trans = g0_val / (1.0 + arg * arg);
    const double g_vib = std::max(0.0, static_cast<double>(total[i]) - g_trans);
    trans_out[i] = static_cast<real_t>(g_trans);
    vib_out[i] = static_cast<real_t>(g_vib);
  }
  vib_out[0] = static_cast<real_t>(0.0);
}

struct FullRawSpectra {
  std::vector<real_t> freqs_thz;
  std::vector<real_t> omega_rad_ps;
  std::vector<real_t> total;
  std::map<std::string, std::vector<real_t>> partials;
};

[[nodiscard]] FullRawSpectra
computeRawFullSpectra(const correlation::analysis::Histogram &vacf_hist, const VDOSParams &params,
                      const correlation::core::Cell *cell) {
  FullRawSpectra raw;
  const real_t time_step = vacf_hist.bins[1] - vacf_hist.bins[0];
  const auto elem_weights = computeElementMassWeights(cell, vacf_hist);

  if (!elem_weights.empty()) {
    for (const auto &[sym, weight] : elem_weights) {
      auto [freqs, partial_dos] = correlation::analysis::DynamicsAnalyzer::calculateVDOS(
          vacf_hist.partials.at(sym), time_step, params.window_type, params.zero_pad_factor);
      if (raw.freqs_thz.empty()) {
        raw.freqs_thz = std::move(freqs);
        raw.total.assign(raw.freqs_thz.size(), static_cast<real_t>(0.0));
      }
      for (size_t i = 0; i < partial_dos.size(); ++i) {
        partial_dos[i] *= weight;
        raw.total[i] += partial_dos[i];
      }
      raw.partials[sym] = std::move(partial_dos);
    }
  } else {
    auto [freqs, total_dos] = correlation::analysis::DynamicsAnalyzer::calculateVDOS(
        vacf_hist.partials.at("Total"), time_step, params.window_type, params.zero_pad_factor);
    raw.freqs_thz = std::move(freqs);
    raw.total = std::move(total_dos);
  }

  raw.omega_rad_ps.resize(raw.freqs_thz.size());
  for (size_t i = 0; i < raw.freqs_thz.size(); ++i) {
    raw.omega_rad_ps[i] = raw.freqs_thz[i] * correlation::math::TWO_PI;
  }
  normalizeVdosSpectra(raw.omega_rad_ps, raw.total, raw.partials);
  return raw;
}

struct SplitSpectra {
  std::vector<real_t> vib_total;
  std::vector<real_t> trans_total;
  std::map<std::string, std::vector<real_t>> vib_partials;
  std::map<std::string, std::vector<real_t>> trans_partials;
};

[[nodiscard]] SplitSpectra decomposeSpectra(const FullRawSpectra &raw, double fluidicity,
                                            bool normalize) {
  SplitSpectra split;
  if (!raw.partials.empty()) {
    split.vib_total.assign(raw.omega_rad_ps.size(), static_cast<real_t>(0.0));
    split.trans_total.assign(raw.omega_rad_ps.size(), static_cast<real_t>(0.0));
    for (const auto &[sym, partial_dos] : raw.partials) {
      std::vector<real_t> vib_p;
      std::vector<real_t> trans_p;
      splitSingleSpectrum(raw.omega_rad_ps, partial_dos, fluidicity, vib_p, trans_p);
      for (size_t i = 0; i < raw.omega_rad_ps.size(); ++i) {
        split.vib_total[i] += vib_p[i];
        split.trans_total[i] += trans_p[i];
      }
      split.vib_partials[sym] = std::move(vib_p);
      split.trans_partials[sym] = std::move(trans_p);
    }
  } else {
    splitSingleSpectrum(raw.omega_rad_ps, raw.total, fluidicity, split.vib_total,
                        split.trans_total);
  }
  if (normalize) {
    normalizeVdosSpectra(raw.freqs_thz, split.vib_total, split.vib_partials);
    normalizeVdosSpectra(raw.freqs_thz, split.trans_total, split.trans_partials);
  }
  return split;
}

[[nodiscard]] correlation::analysis::Histogram
buildVdosHistogram(std::string title, std::string y_label, std::string file_suffix, bool normalized,
                   const std::vector<real_t> &out_thz, const std::vector<real_t> &out_cm_inv,
                   const std::vector<real_t> &out_mev, const std::vector<real_t> &full_total,
                   const std::map<std::string, std::vector<real_t>> &full_partials,
                   size_t keep_count) {
  correlation::analysis::Histogram hist;
  hist.x_label = "ω";
  hist.title = std::move(title);
  hist.description = hist.title;
  hist.file_suffix = std::move(file_suffix);
  hist.y_label = std::move(y_label);
  hist.x_unit = "THz";
  hist.y_unit = normalized ? "THz⁻¹" : "arbitrary units";
  hist.bins = out_thz;
  hist.partials["Total"].assign(full_total.begin(),
                                full_total.begin() + static_cast<std::ptrdiff_t>(keep_count));
  for (const auto &[sym, partial_dos] : full_partials) {
    hist.partials[sym].assign(partial_dos.begin(),
                              partial_dos.begin() + static_cast<std::ptrdiff_t>(keep_count));
  }
  hist.partials["Frequency_cm_1"] = out_cm_inv;
  hist.partials["Frequency_meV"] = out_mev;
  return hist;
}

} // namespace

real_t VDOSCalculator::computeRecommendedMaxFrequency(const correlation::core::Cell &cell,
                                                      real_t time_step) noexcept {
  real_t min_mass = std::numeric_limits<real_t>::max();
  for (const auto &elem : cell.elements()) {
    if (const auto mass = correlation::physics::findAtomicMass(elem.symbol);
        mass.has_value() && *mass > static_cast<real_t>(0.0) && *mass < min_mass) {
      min_mass = *mass;
    }
  }

  // Continuous mass scaling: ν_max = max(10.0, 50.0 / sqrt(m_min))
  // For H (m ~ 1.008): ~49.8 THz
  // For C (m ~ 12.01): ~14.4 THz
  // For Si (m ~ 28.09): ~9.4 THz -> clamped to 10.0 THz minimum
  auto max_freq = static_cast<real_t>(10.0);
  if (min_mass < std::numeric_limits<real_t>::max() && min_mass > static_cast<real_t>(0.0)) {
    const real_t scaled_freq = static_cast<real_t>(50.0) / std::sqrt(min_mass);
    max_freq = std::max(static_cast<real_t>(10.0), scaled_freq);
  }

  // Bounded by Nyquist frequency: ν_nyquist = 1000 / (2 * time_step) THz (time_step in fs)
  if (time_step > static_cast<real_t>(0.0)) {
    const real_t nyquist_thz = static_cast<real_t>(1000.0) / (static_cast<real_t>(2.0) * time_step);
    max_freq = std::min(max_freq, nyquist_thz);
  }

  return max_freq;
}

void VDOSCalculator::calculateTrajectory(
    correlation::analysis::DistributionFunctions &dists, const correlation::core::Trajectory &traj,
    const correlation::analysis::AnalysisSettings &settings) const {
  const auto &all = dists.getAllHistograms();
  if (!all.contains("VACF")) {
    return; // VACF must be computed first
  }
  const auto &vacf_hist = dists.getHistogram("VACF");
  auto const it_total = vacf_hist.partials.find("Total");
  if (it_total == vacf_hist.partials.end() || it_total->second.size() < 2) {
    return; // VACF data is too short for VDOS calculation (e.g. single-frame calculation)
  }

  real_t const time_step =
      vacf_hist.bins.size() >= 2 ? (vacf_hist.bins[1] - vacf_hist.bins[0]) : traj.getTimeStep();

  VDOSParams params;
  params.window_type = settings.vdos_window;
  params.zero_pad_factor = settings.vdos_zero_pad;
  const correlation::core::Cell *ref_cell = nullptr;
  if (traj.getFrameCount() > 0) {
    ref_cell = &traj.firstFrame();
    params.max_real_freq = computeRecommendedMaxFrequency(*ref_cell, time_step);
  }

  auto decomposed = calculateDecomposed(vacf_hist, params, ref_cell);
  dists.addHistogram("VDOS", std::move(decomposed.vibrational));
  dists.addHistogram("VDOS_Translational", std::move(decomposed.translational));
}

VDOSDecomposition
VDOSCalculator::calculateDecomposed(const correlation::analysis::Histogram &vacf_hist,
                                    const VDOSParams &params, const correlation::core::Cell *cell) {
  const auto &vacf_total = vacf_hist.partials.at("Total");
  if (vacf_total.size() < 2) {
    throw std::logic_error("VACF data is too short for VDOS calculation.");
  }

  const FullRawSpectra raw = computeRawFullSpectra(vacf_hist, params, cell);
  const double g0_omega = raw.total.empty() ? 0.0 : static_cast<double>(raw.total[0]);
  const double fluidicity = computeFluidicity(g0_omega, static_cast<double>(vacf_total[0]),
                                              raw.omega_rad_ps, raw.total, cell);
  const SplitSpectra split = decomposeSpectra(raw, fluidicity, params.normalize);

  size_t keep_count = 0;
  while (keep_count < raw.freqs_thz.size() && raw.freqs_thz[keep_count] <= params.max_real_freq) {
    ++keep_count;
  }

  const std::vector<real_t> out_thz(
      raw.freqs_thz.begin(), raw.freqs_thz.begin() + static_cast<std::ptrdiff_t>(keep_count));
  std::vector<real_t> out_cm_inv(keep_count);
  std::vector<real_t> out_mev(keep_count);
  for (size_t i = 0; i < keep_count; ++i) {
    out_cm_inv[i] = out_thz[i] * static_cast<real_t>(correlation::math::THZ_TO_CMINV);
    out_mev[i] = out_thz[i] * static_cast<real_t>(correlation::math::THZ_TO_MEV);
  }

  return VDOSDecomposition{
      .vibrational = buildVdosHistogram("Vibrational Density of States", "g_vib(ω)", "_VDOS",
                                        params.normalize, out_thz, out_cm_inv, out_mev,
                                        split.vib_total, split.vib_partials, keep_count),
      .translational = buildVdosHistogram(
          "Translational Density of States", "g_trans(ω)", "_VDOS_trans", params.normalize, out_thz,
          out_cm_inv, out_mev, split.trans_total, split.trans_partials, keep_count),
      .fluidicity = static_cast<real_t>(fluidicity),
  };
}

correlation::analysis::Histogram
VDOSCalculator::calculate(const correlation::analysis::Histogram &vacf_hist,
                          const VDOSParams &params, const correlation::core::Cell *cell) {
  return calculateDecomposed(vacf_hist, params, cell).vibrational;
}

} // namespace correlation::calculators
