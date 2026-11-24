/**
 * @file EliashbergCalculator.cpp
 * @brief Implementation of Hybrid MLIP-VDOS Eliashberg α²F(ω) and Allen-Dynes T_c estimator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/ml/EliashbergCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "calculators/scattering/LDOSPostProcessor.hpp"
#include "calculators/temporal/VDOSCalculator.hpp"
#include "mlip/OrbV3GNNModel.hpp"
#include "physics/PhysicalData.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace correlation::calculators {

namespace {

const bool REGISTERED =
    CalculatorFactory::registerTypeSafe<EliashbergCalculator>("EliashbergCalculator");

constexpr double K_THZ_TO_MEV = 4.135667696;
constexpr double K_MEV_TO_KELVIN = 11.60451812;
constexpr double K_HBAR2_OVER_2AMU_EV_A2 = 2.09008e-3;
constexpr double K_MIN_POSITIVE = 1e-12;

struct TrajectorySamples {
  std::vector<std::vector<double>> fermi_frames;
  std::vector<std::vector<math::Vector3<double>>> pos_frames;
};

[[nodiscard]] std::vector<real_t> buildOutputEnergyGrid(const mlip::MLIPOutput &output) {
  const size_t n_bins = output.ldos_bins;
  if (n_bins == 0) {
    return {};
  }
  const bool has_bounds = output.ldos_e_min < output.ldos_e_max;
  const real_t e_min = has_bounds ? output.ldos_e_min : static_cast<real_t>(-30.0);
  const real_t e_max = has_bounds ? output.ldos_e_max : static_cast<real_t>(15.0);
  std::vector<real_t> bins(n_bins);
  const real_t delta = (e_max - e_min) / static_cast<real_t>(n_bins);
  for (size_t i = 0; i < n_bins; ++i) {
    bins[i] = e_min + (static_cast<real_t>(i) + static_cast<real_t>(0.5)) * delta;
  }
  return bins;
}

[[nodiscard]] std::vector<double> extractFrameFermiDos(const core::Cell &cell,
                                                       const mlip::MLIPOutput &output,
                                                       std::span<const real_t> energies,
                                                       real_t electronic_temp_k) {
  const size_t n_atoms = cell.atoms().size();
  std::vector<double> fermi(n_atoms, 0.0);
  if (energies.empty()) {
    return fermi;
  }
  for (size_t i = 0; i < n_atoms; ++i) {
    const auto atom_ldos = output.getLdos(i);
    if (atom_ldos.size() < energies.size()) {
      continue;
    }
    fermi[i] = static_cast<double>(LDOSPostProcessor::evaluateFermiDos(
        energies, atom_ldos.subspan(0, energies.size()), electronic_temp_k));
  }
  return fermi;
}

[[nodiscard]] std::vector<size_t> selectSampleFrameIndices(size_t total_frames, size_t max_frames) {
  if (total_frames == 0 || max_frames == 0) {
    return {};
  }
  const size_t count = std::min(total_frames, max_frames);
  std::vector<size_t> indices;
  indices.reserve(count);
  if (count == 1) {
    indices.push_back(0);
    return indices;
  }
  for (size_t i = 0; i < count; ++i) {
    indices.push_back((i * (total_frames - 1)) / (count - 1));
  }
  return indices;
}

[[nodiscard]] TrajectorySamples collectTrajectorySamples(const core::Trajectory &traj,
                                                         const mlip::MLIPInterface &model,
                                                         const EliashbergParams &params,
                                                         size_t n_atoms,
                                                         const std::atomic<bool> *cancel_flag) {
  TrajectorySamples samples;
  const auto indices = selectSampleFrameIndices(traj.getFrameCount(), params.max_coupling_frames);
  samples.fermi_frames.reserve(indices.size());
  samples.pos_frames.reserve(indices.size());

  std::vector<real_t> energies;
  for (const size_t f_idx : indices) {
    if (cancel_flag != nullptr && cancel_flag->load()) {
      return {};
    }
    const auto &cell = traj.getFrame(f_idx);
    if (cell.atoms().size() != n_atoms) {
      continue;
    }
    const auto out = model.evaluate(cell);
    if (out.ldos_bins == 0 || out.ldos.empty()) {
      continue;
    }
    if (energies.empty()) {
      energies = buildOutputEnergyGrid(out);
    }
    samples.fermi_frames.push_back(
        extractFrameFermiDos(cell, out, energies, params.electronic_temp_k));

    std::vector<math::Vector3<double>> frame_pos(n_atoms);
    for (size_t i = 0; i < n_atoms; ++i) {
      const auto &pos = cell.atoms()[i].position();
      frame_pos[i] = {static_cast<double>(pos.x()), static_cast<double>(pos.y()),
                      static_cast<double>(pos.z())};
    }
    samples.pos_frames.push_back(std::move(frame_pos));
  }
  return samples;
}

[[nodiscard]] double computeAtomSlopeSquared(const TrajectorySamples &samples, size_t atom_idx,
                                             double &mean_n_out) noexcept {
  const size_t n_samples = samples.fermi_frames.size();
  const double inv_samples = 1.0 / static_cast<double>(n_samples);
  double mean_n = 0.0;
  math::Vector3<double> mean_r{0.0, 0.0, 0.0};
  for (size_t sample_idx = 0; sample_idx < n_samples; ++sample_idx) {
    mean_n += samples.fermi_frames[sample_idx][atom_idx] * inv_samples;
    mean_r.x() += samples.pos_frames[sample_idx][atom_idx].x() * inv_samples;
    mean_r.y() += samples.pos_frames[sample_idx][atom_idx].y() * inv_samples;
    mean_r.z() += samples.pos_frames[sample_idx][atom_idx].z() * inv_samples;
  }
  mean_n_out = mean_n;
  if (n_samples < 2) {
    return -1.0;
  }

  double sum_dn2 = 0.0;
  double sum_dr2 = 0.0;
  for (size_t sample_idx = 0; sample_idx < n_samples; ++sample_idx) {
    const double delta_n = samples.fermi_frames[sample_idx][atom_idx] - mean_n;
    const double delta_x = samples.pos_frames[sample_idx][atom_idx].x() - mean_r.x();
    const double delta_y = samples.pos_frames[sample_idx][atom_idx].y() - mean_r.y();
    const double delta_z = samples.pos_frames[sample_idx][atom_idx].z() - mean_r.z();
    sum_dn2 += delta_n * delta_n;
    sum_dr2 += delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
  }
  return (sum_dr2 > 1e-10) ? (sum_dn2 / sum_dr2) : -1.0;
}

[[nodiscard]] double evaluateAtomFiniteGradientSq(const core::Cell &ref_cell, size_t atom_idx,
                                                  const mlip::MLIPInterface &model,
                                                  std::span<const real_t> energies, double delta,
                                                  real_t electronic_temp_k) {
  const double inv_2delta = 1.0 / (2.0 * delta);
  double grad_sq = 0.0;
  for (int axis = 0; axis < 3; ++axis) {
    core::Cell cell_plus = ref_cell;
    core::Cell cell_minus = ref_cell;
    auto pos_p = cell_plus.atoms()[atom_idx].position();
    auto pos_m = cell_minus.atoms()[atom_idx].position();
    const auto d_real = static_cast<real_t>(delta);
    if (axis == 0) {
      pos_p.x() += d_real;
      pos_m.x() -= d_real;
    } else if (axis == 1) {
      pos_p.y() += d_real;
      pos_m.y() -= d_real;
    } else {
      pos_p.z() += d_real;
      pos_m.z() -= d_real;
    }
    cell_plus.atoms()[atom_idx].setPosition(pos_p);
    cell_minus.atoms()[atom_idx].setPosition(pos_m);

    const auto out_p = model.evaluate(cell_plus);
    const auto out_m = model.evaluate(cell_minus);
    const auto fermi_plus = extractFrameFermiDos(cell_plus, out_p, energies, electronic_temp_k);
    const auto fermi_minus = extractFrameFermiDos(cell_minus, out_m, energies, electronic_temp_k);
    const double deriv = (fermi_plus[atom_idx] - fermi_minus[atom_idx]) * inv_2delta;
    grad_sq += deriv * deriv;
  }
  return grad_sq;
}

[[nodiscard]] double computeMeanFrequencyMev(std::span<const double> omega_mev,
                                             std::span<const double> spectrum) noexcept {
  if (omega_mev.size() < 2 || omega_mev.size() != spectrum.size()) {
    return 10.0;
  }
  double num = 0.0;
  double den = 0.0;
  for (size_t i = 0; i + 1 < omega_mev.size(); ++i) {
    const double d_w = omega_mev[i + 1] - omega_mev[i];
    const double dos_lo = std::max(0.0, spectrum[i]);
    const double dos_hi = std::max(0.0, spectrum[i + 1]);
    num += 0.5 * (omega_mev[i] * dos_lo + omega_mev[i + 1] * dos_hi) * d_w;
    den += 0.5 * (dos_lo + dos_hi) * d_w;
  }
  return (den > K_MIN_POSITIVE) ? std::max(1.0, num / den) : 10.0;
}

[[nodiscard]] double lookupSpeciesMass(const std::string &symbol) noexcept {
  if (const auto m_opt = physics::findAtomicMass(symbol);
      m_opt.has_value() && *m_opt > static_cast<real_t>(0.0)) {
    return static_cast<double>(*m_opt);
  }
  return 12.0;
}

[[nodiscard]] double computeSpeciesPrefactor(const std::string &symbol, double eta_alpha,
                                             double mean_omega_mev) noexcept {
  const double mass_amu = lookupSpeciesMass(symbol);
  const double mean_omega_ev = std::max(1e-3, mean_omega_mev * 1e-3);
  return std::max(0.0, (K_HBAR2_OVER_2AMU_EV_A2 * eta_alpha) / (mass_amu * mean_omega_ev));
}

[[nodiscard]] std::vector<double> convertVdosToMevDensity(std::span<const real_t> g_vals,
                                                          double unit_to_mev) {
  std::vector<double> g_mev(g_vals.size());
  for (size_t i = 0; i < g_vals.size(); ++i) {
    g_mev[i] = std::max(0.0, static_cast<double>(g_vals[i])) / unit_to_mev;
  }
  return g_mev;
}

void synthesizeSpeciesSpectra(const analysis::Histogram &vdos_vib,
                              std::span<const double> omega_mev_d, double unit_to_mev,
                              const SpeciesCouplingWeights &coupling,
                              std::vector<double> &total_a2f,
                              std::map<std::string, std::vector<real_t>> &species_a2f) {
  const size_t n_pts = omega_mev_d.size();
  for (const auto &[sym, g_vals] : vdos_vib.partials) {
    if (sym == "Total" || sym == "total" || sym.starts_with("Frequency_") ||
        g_vals.size() != n_pts) {
      continue;
    }
    const auto eta_it = coupling.eta_by_species.find(sym);
    const double eta_alpha =
        (eta_it != coupling.eta_by_species.end()) ? static_cast<double>(eta_it->second) : 0.0;
    const auto g_mev = convertVdosToMevDensity(g_vals, unit_to_mev);
    const double prefactor =
        computeSpeciesPrefactor(sym, eta_alpha, computeMeanFrequencyMev(omega_mev_d, g_mev));

    std::vector<real_t> partial_a2f(n_pts, static_cast<real_t>(0.0));
    for (size_t i = 0; i < n_pts; ++i) {
      const double val = prefactor * g_mev[i] * omega_mev_d[i];
      partial_a2f[i] = static_cast<real_t>(val);
      total_a2f[i] += val;
    }
    species_a2f[sym] = std::move(partial_a2f);
  }
}

void synthesizeFallbackTotalSpectrum(const analysis::Histogram &vdos_vib,
                                     std::span<const double> omega_mev_d, double unit_to_mev,
                                     const SpeciesCouplingWeights &coupling,
                                     std::vector<double> &total_a2f) {
  const auto it_tot = vdos_vib.partials.contains("Total") ? vdos_vib.partials.find("Total")
                                                          : vdos_vib.partials.find("total");
  if (it_tot == vdos_vib.partials.end() || it_tot->second.size() != omega_mev_d.size()) {
    return;
  }
  double avg_eta = 0.0;
  double avg_mass = 0.0;
  for (const auto &[sym, eta_val] : coupling.eta_by_species) {
    avg_eta += static_cast<double>(eta_val);
    avg_mass += lookupSpeciesMass(sym);
  }
  const size_t n_spec = std::max<size_t>(1, coupling.eta_by_species.size());
  avg_eta /= static_cast<double>(n_spec);
  avg_mass = (avg_mass > 0.0) ? (avg_mass / static_cast<double>(n_spec)) : 12.0;

  const auto g_mev = convertVdosToMevDensity(it_tot->second, unit_to_mev);
  const double mean_w_ev = std::max(1e-3, computeMeanFrequencyMev(omega_mev_d, g_mev) * 1e-3);
  const double prefactor = (K_HBAR2_OVER_2AMU_EV_A2 * avg_eta) / (avg_mass * mean_w_ev);
  for (size_t i = 0; i < omega_mev_d.size(); ++i) {
    total_a2f[i] = std::max(0.0, prefactor * g_mev[i] * omega_mev_d[i]);
  }
}

void computeEliashbergMoments(std::span<const real_t> omega_mev, std::span<const real_t> a2f,
                              double omega_min_mev, double &lambda_out, double &omega_log_out,
                              double &omega_2_out) noexcept {
  lambda_out = 0.0;
  omega_log_out = 0.0;
  omega_2_out = 0.0;
  if (omega_mev.size() < 2 || omega_mev.size() != a2f.size()) {
    return;
  }

  double int_inv_w = 0.0;
  double int_ln_w = 0.0;
  double int_w = 0.0;
  const double w_cut = std::max(1e-4, omega_min_mev);

  for (size_t i = 0; i + 1 < omega_mev.size(); ++i) {
    const auto omega_lo = static_cast<double>(omega_mev[i]);
    const auto omega_hi = static_cast<double>(omega_mev[i + 1]);
    if (omega_hi <= w_cut || omega_hi <= omega_lo) {
      continue;
    }
    const double spec_lo = (omega_lo >= w_cut) ? std::max(0.0, static_cast<double>(a2f[i])) : 0.0;
    const double spec_hi = std::max(0.0, static_cast<double>(a2f[i + 1]));
    const double delta_omega = omega_hi - omega_lo;
    const double eff_w0 = std::max(omega_lo, w_cut);

    int_inv_w += 0.5 * ((spec_lo / eff_w0) + (spec_hi / omega_hi)) * delta_omega;
    int_ln_w +=
        0.5 * ((spec_lo * std::log(eff_w0) / eff_w0) + (spec_hi * std::log(omega_hi) / omega_hi)) *
        delta_omega;
    int_w += 0.5 * (spec_lo * eff_w0 + spec_hi * omega_hi) * delta_omega;
  }

  lambda_out = 2.0 * int_inv_w;
  if (lambda_out > K_MIN_POSITIVE) {
    omega_log_out = std::exp((2.0 / lambda_out) * int_ln_w);
    omega_2_out = std::sqrt(std::max(0.0, (2.0 / lambda_out) * int_w));
  }
}

void populateCriticalTemperatures(double lambda, double omega_log_mev, double omega_2_mev,
                                  double mu_star,
                                  analysis::SuperconductingDiagnostics &diag) noexcept {
  const double omega_log_k = omega_log_mev * K_MEV_TO_KELVIN;
  diag.lambda = static_cast<real_t>(lambda);
  diag.omega_log_mev = static_cast<real_t>(omega_log_mev);
  diag.omega_log_k = static_cast<real_t>(omega_log_k);
  diag.omega_2_mev = static_cast<real_t>(omega_2_mev);
  diag.mu_star = static_cast<real_t>(mu_star);
  diag.has_eliashberg = true;

  const double denom = lambda - mu_star * (1.0 + 0.62 * lambda);
  if (lambda <= K_MIN_POSITIVE || omega_log_k <= K_MIN_POSITIVE || denom <= 1e-6) {
    diag.tc_mcmillan_k = static_cast<real_t>(0.0);
    diag.tc_allen_dynes_k = static_cast<real_t>(0.0);
    return;
  }

  const double exponent = -1.04 * (1.0 + lambda) / denom;
  const double tc_mcmillan = (omega_log_k / 1.20) * std::exp(exponent);

  const double w2_ratio = (omega_log_mev > K_MIN_POSITIVE) ? (omega_2_mev / omega_log_mev) : 1.0;
  const double lambda1 = 2.46 * (1.0 + 3.8 * mu_star);
  const double lambda2 = 1.82 * (1.0 + 6.3 * mu_star) * w2_ratio;
  const double strong_factor = std::cbrt(1.0 + std::pow(lambda / lambda1, 1.5));
  const double shape_factor =
      1.0 + ((w2_ratio - 1.0) * lambda * lambda) / (lambda * lambda + lambda2 * lambda2);

  diag.tc_mcmillan_k = static_cast<real_t>(std::max(0.0, tc_mcmillan));
  diag.tc_allen_dynes_k =
      static_cast<real_t>(std::max(0.0, strong_factor * shape_factor * tc_mcmillan));
}

void buildSpectralHistogramMetadata(analysis::Histogram &hist, std::vector<real_t> omega_mev,
                                    int frame_count) {
  hist.bins = std::move(omega_mev);
  hist.title = "Eliashberg Spectral Function α²F(ω)";
  hist.x_label = "Phonon Energy ħω";
  hist.y_label = "α²F(ω)";
  hist.x_unit = "meV";
  hist.y_unit = "dimensionless";
  hist.description = "Hybrid MLIP-VDOS Eliashberg Spectral Function α²F(ω) estimate";
  hist.file_suffix = "_a2f";
  hist.compute_count = frame_count;
}

} // namespace

SpeciesCouplingWeights TrajectoryFluctuationCouplingEstimator::estimate(
    const core::Trajectory &traj, const mlip::MLIPInterface &model, const EliashbergParams &params,
    const std::atomic<bool> *cancel_flag) const {
  SpeciesCouplingWeights weights;
  if (traj.getFrameCount() == 0) {
    return weights;
  }

  const auto &ref_cell = traj.getFrame(0);
  const size_t n_atoms = ref_cell.atoms().size();
  if (n_atoms == 0) {
    return weights;
  }

  const auto samples = collectTrajectorySamples(traj, model, params, n_atoms, cancel_flag);
  if (samples.fermi_frames.empty()) {
    return weights;
  }

  std::map<std::string, double> species_fermi_sum;
  std::map<std::string, size_t> species_atom_counts;
  std::map<std::string, double> species_slope_sq_sum;
  std::map<std::string, size_t> species_slope_counts;
  double total_fermi_acc = 0.0;

  for (size_t i = 0; i < n_atoms; ++i) {
    const std::string &sym = ref_cell.atoms()[i].element().symbol;
    double mean_n = 0.0;
    const double slope_sq = computeAtomSlopeSquared(samples, i, mean_n);
    species_fermi_sum[sym] += mean_n;
    species_atom_counts[sym] += 1;
    total_fermi_acc += mean_n;
    if (slope_sq >= 0.0) {
      species_slope_sq_sum[sym] += slope_sq;
      species_slope_counts[sym] += 1;
    }
  }

  weights.total_fermi_dos = static_cast<real_t>(total_fermi_acc / static_cast<double>(n_atoms));
  for (const auto &[sym, count] : species_atom_counts) {
    const double n_alpha = species_fermi_sum[sym] / static_cast<double>(count);
    weights.fermi_dos_by_species[sym] = static_cast<real_t>(n_alpha);
    const size_t slope_count = species_slope_counts[sym];
    const double mean_grad_sq =
        (slope_count > 0) ? (species_slope_sq_sum[sym] / static_cast<double>(slope_count)) : 0.0;
    weights.eta_by_species[sym] =
        static_cast<real_t>((n_alpha > 1e-6) ? (mean_grad_sq / n_alpha) : 0.0);
  }
  return weights;
}

SpeciesCouplingWeights FiniteDisplacementCouplingEstimator::estimate(
    const core::Trajectory &traj, const mlip::MLIPInterface &model, const EliashbergParams &params,
    const std::atomic<bool> *cancel_flag) const {
  SpeciesCouplingWeights weights;
  if (traj.getFrameCount() == 0) {
    return weights;
  }

  const auto &ref_cell = traj.getFrame(0);
  const size_t n_atoms = ref_cell.atoms().size();
  if (n_atoms == 0) {
    return weights;
  }

  const auto ref_out = model.evaluate(ref_cell);
  if (ref_out.ldos_bins == 0 || ref_out.ldos.empty()) {
    return weights;
  }

  const auto energies = buildOutputEnergyGrid(ref_out);
  const auto ref_fermi =
      extractFrameFermiDos(ref_cell, ref_out, energies, params.electronic_temp_k);
  const double delta = std::max(1e-4, static_cast<double>(params.displacement_step));

  std::map<std::string, std::vector<size_t>> atoms_by_species;
  std::map<std::string, double> species_fermi_sum;
  double total_fermi_sum = 0.0;
  for (size_t i = 0; i < n_atoms; ++i) {
    const std::string &sym = ref_cell.atoms()[i].element().symbol;
    atoms_by_species[sym].push_back(i);
    species_fermi_sum[sym] += ref_fermi[i];
    total_fermi_sum += ref_fermi[i];
  }

  weights.total_fermi_dos = static_cast<real_t>(total_fermi_sum / static_cast<double>(n_atoms));
  const size_t max_per_species = std::max<size_t>(1, params.samples_per_species);

  for (const auto &[sym, atom_indices] : atoms_by_species) {
    const double n_alpha = species_fermi_sum[sym] / static_cast<double>(atom_indices.size());
    weights.fermi_dos_by_species[sym] = static_cast<real_t>(n_alpha);

    const size_t sample_count = std::min(atom_indices.size(), max_per_species);
    double grad_sq_acc = 0.0;
    size_t evaluated_count = 0;
    for (size_t sample_idx = 0; sample_idx < sample_count; ++sample_idx) {
      if (cancel_flag != nullptr && cancel_flag->load()) {
        return weights;
      }
      grad_sq_acc += evaluateAtomFiniteGradientSq(ref_cell, atom_indices[sample_idx], model,
                                                  energies, delta, params.electronic_temp_k);
      ++evaluated_count;
    }

    const double mean_grad_sq =
        (evaluated_count > 0) ? (grad_sq_acc / static_cast<double>(evaluated_count)) : 0.0;
    weights.eta_by_species[sym] =
        static_cast<real_t>((n_alpha > 1e-6) ? (mean_grad_sq / n_alpha) : 0.0);
  }
  return weights;
}

bool EliashbergCalculator::isConfigured() const {
  return model_ != nullptr || !mlip::OrbV3GNNModel::findDefaultModelPath().empty();
}

std::unique_ptr<IElectronPhononCouplingEstimator>
EliashbergCalculator::createCouplingEstimator(analysis::CouplingEstimationMode mode) {
  if (mode == analysis::CouplingEstimationMode::FiniteDisplacement) {
    return std::make_unique<FiniteDisplacementCouplingEstimator>();
  }
  return std::make_unique<TrajectoryFluctuationCouplingEstimator>();
}

analysis::SuperconductingDiagnostics
EliashbergCalculator::computeFromSpectrum(std::span<const real_t> omega_mev,
                                          std::span<const real_t> a2f, real_t mu_star,
                                          real_t omega_min_mev) noexcept {
  analysis::SuperconductingDiagnostics diag{};
  if (omega_mev.size() < 2 || omega_mev.size() != a2f.size()) {
    return diag;
  }

  double lambda = 0.0;
  double omega_log_mev = 0.0;
  double omega_2_mev = 0.0;
  computeEliashbergMoments(omega_mev, a2f, static_cast<double>(omega_min_mev), lambda,
                           omega_log_mev, omega_2_mev);
  populateCriticalTemperatures(lambda, omega_log_mev, omega_2_mev, static_cast<double>(mu_star),
                               diag);
  return diag;
}

EliashbergResult EliashbergCalculator::calculate(const core::Trajectory &traj,
                                                 const analysis::Histogram &vdos_vib,
                                                 const EliashbergParams &params,
                                                 const IElectronPhononCouplingEstimator *estimator,
                                                 const std::atomic<bool> *cancel_flag) {
  EliashbergResult result;
  if (params.model == nullptr || traj.getFrameCount() == 0 || vdos_vib.bins.size() < 2) {
    return result;
  }

  std::unique_ptr<IElectronPhononCouplingEstimator> owned_estimator;
  const IElectronPhononCouplingEstimator *active_estimator = estimator;
  if (active_estimator == nullptr) {
    owned_estimator = createCouplingEstimator(params.coupling_mode);
    active_estimator = owned_estimator.get();
  }

  result.coupling_weights = active_estimator->estimate(traj, *params.model, params, cancel_flag);
  if (cancel_flag != nullptr && cancel_flag->load()) {
    return result;
  }

  const double unit_to_mev = (vdos_vib.x_unit == "meV") ? 1.0 : K_THZ_TO_MEV;
  const size_t n_pts = vdos_vib.bins.size();
  std::vector<double> omega_mev_d(n_pts);
  std::vector<real_t> omega_mev_r(n_pts);
  for (size_t i = 0; i < n_pts; ++i) {
    omega_mev_d[i] = static_cast<double>(vdos_vib.bins[i]) * unit_to_mev;
    omega_mev_r[i] = static_cast<real_t>(omega_mev_d[i]);
  }

  std::vector<double> total_a2f(n_pts, 0.0);
  std::map<std::string, std::vector<real_t>> species_a2f;
  synthesizeSpeciesSpectra(vdos_vib, omega_mev_d, unit_to_mev, result.coupling_weights, total_a2f,
                           species_a2f);
  if (species_a2f.empty()) {
    synthesizeFallbackTotalSpectrum(vdos_vib, omega_mev_d, unit_to_mev, result.coupling_weights,
                                    total_a2f);
  }

  std::vector<real_t> total_a2f_r(n_pts, static_cast<real_t>(0.0));
  for (size_t i = 0; i < n_pts; ++i) {
    total_a2f_r[i] = static_cast<real_t>(total_a2f[i]);
  }

  buildSpectralHistogramMetadata(result.spectral_function, std::move(omega_mev_r),
                                 static_cast<int>(traj.getFrameCount()));
  result.spectral_function.partials["Total"] = total_a2f_r;
  for (auto &[sym, vals] : species_a2f) {
    result.spectral_function.partials[sym] = std::move(vals);
  }

  result.diagnostics = computeFromSpectrum(result.spectral_function.bins, total_a2f_r,
                                           params.mu_star, params.omega_min_mev);
  return result;
}

void EliashbergCalculator::calculateTrajectory(analysis::DistributionFunctions &dists,
                                               const core::Trajectory &traj,
                                               const analysis::AnalysisSettings &settings) const {
  if (model_ == nullptr || traj.getFrameCount() == 0) {
    return;
  }

  const auto &hists = dists.getAllHistograms();
  auto vdos_it = hists.find("VDOS");
  if (vdos_it == hists.end()) {
    vdos_it = hists.find("vDoS");
  }
  if (vdos_it == hists.end()) {
    return;
  }

  EliashbergParams params;
  params.model = model_;
  params.mu_star = settings.eliashberg_mu_star;
  params.coupling_mode = settings.eliashberg_coupling_mode;
  params.displacement_step = settings.eliashberg_disp_step;
  params.electronic_temp_k = settings.ldos_electronic_temp_k;

  auto res = calculate(traj, vdos_it->second, params, nullptr, settings.cancel_flag);
  if (!res.spectral_function.bins.empty()) {
    dists.addHistogram("Eliashberg", std::move(res.spectral_function));
  }
  if (res.diagnostics.has_eliashberg) {
    dists.setSuperconductingDiagnostics(res.diagnostics);
  }
}

} // namespace correlation::calculators
