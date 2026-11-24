/**
 * @file TDOSCalculator.cpp
 * @brief Implementation of Total Density of States (TDOS) and site-resolved LDOS calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/scattering/TDOSCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "calculators/scattering/LDOSPostProcessor.hpp"
#include "mlip/OrbV3GNNModel.hpp"

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace correlation::calculators {

bool TDOSCalculator::isConfigured() const {
  return model_ != nullptr || !correlation::mlip::OrbV3GNNModel::findDefaultModelPath().empty();
}

namespace {

const bool REGISTERED = CalculatorFactory::registerTypeSafe<TDOSCalculator>("TDOSCalculator");

[[nodiscard]] std::vector<real_t> buildEnergyGrid(real_t e_min, real_t e_max, size_t num_bins) {
  if (num_bins == 0) {
    return {};
  }
  std::vector<real_t> bins(num_bins);
  const real_t delta_e = (e_max - e_min) / static_cast<real_t>(num_bins);
  for (size_t i = 0; i < num_bins; ++i) {
    bins[i] = e_min + (static_cast<real_t>(i) + static_cast<real_t>(0.5)) * delta_e;
  }
  return bins;
}

void resolveEnergyBounds(const correlation::mlip::MLIPOutput &output, const TDOSParams &params,
                         real_t &grid_min, real_t &grid_max) noexcept {
  const bool use_model_grid = (output.ldos_e_min < output.ldos_e_max) &&
                              (params.e_min == static_cast<real_t>(-30.0)) &&
                              (params.e_max == static_cast<real_t>(15.0));
  grid_min = use_model_grid ? output.ldos_e_min : params.e_min;
  grid_max = use_model_grid ? output.ldos_e_max : params.e_max;
}

void populateHistogramMetadata(correlation::analysis::Histogram &hist, std::vector<real_t> bins) {
  hist.bins = std::move(bins);
  hist.title = "Total Density of States (TDOS)";
  hist.x_label = "Energy";
  hist.y_label = "DOS";
  hist.x_unit = "eV";
  hist.y_unit = "states/eV/atom";
  hist.description = "Total Density of States (TDOS)";
  hist.file_suffix = "_tdos";
}

void populateIdosHistogram(const correlation::analysis::Histogram &tdos_hist,
                           correlation::analysis::Histogram &idos_hist) {
  idos_hist.bins = tdos_hist.bins;
  idos_hist.title = "Integrated Density of States (IDOS)";
  idos_hist.x_label = "Energy";
  idos_hist.y_label = "IDOS";
  idos_hist.x_unit = "eV";
  idos_hist.y_unit = "states/atom";
  idos_hist.description = "Cumulative Integrated Density of States (IDOS)";
  idos_hist.file_suffix = "_idos";
  idos_hist.compute_count = tdos_hist.compute_count;

  for (const auto &[key, vals] : tdos_hist.partials) {
    idos_hist.partials[key] = LDOSPostProcessor::computeIntegratedDos(tdos_hist.bins, vals);
  }
}

void accumulateFrameLdos(const correlation::core::Cell &cell,
                         const correlation::mlip::MLIPOutput &output, size_t num_bins,
                         std::vector<double> &acc_total,
                         std::map<std::string, std::vector<double>> &acc_species) {
  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();
  if (num_atoms == 0) {
    return;
  }

  const double norm_factor = 1.0 / static_cast<double>(num_atoms);
  for (size_t i = 0; i < num_atoms; ++i) {
    const auto atom_ldos = output.getLdos(i);
    if (atom_ldos.size() < num_bins) {
      continue;
    }
    const std::string &elem_symbol = atoms[i].element().symbol;
    auto &species_buffer = acc_species[elem_symbol];
    if (species_buffer.size() < num_bins) {
      species_buffer.assign(num_bins, 0.0);
    }

    for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
      const double val = static_cast<double>(atom_ldos[bin_idx]) * norm_factor;
      acc_total[bin_idx] += val;
      species_buffer[bin_idx] += val;
    }
  }
}

void finalizeRawHistogram(correlation::analysis::Histogram &hist, double frame_scale,
                          const std::vector<double> &acc_total,
                          const std::map<std::string, std::vector<double>> &acc_species) {
  const size_t num_bins = hist.bins.size();
  auto &total_vec = hist.partials["total"];
  total_vec.resize(num_bins, static_cast<real_t>(0.0));
  for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
    total_vec[bin_idx] = static_cast<real_t>(acc_total[bin_idx] * frame_scale);
  }

  for (const auto &[species, vals] : acc_species) {
    auto &spec_vec = hist.partials[species];
    spec_vec.resize(num_bins, static_cast<real_t>(0.0));
    for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
      spec_vec[bin_idx] = static_cast<real_t>(vals[bin_idx] * frame_scale);
    }
  }
}

void applyPostProcessing(correlation::analysis::Histogram &hist, const TDOSParams &params) {
  if (params.broadening_kernel != correlation::analysis::LDOSBroadeningKernel::None &&
      params.broadening_sigma > static_cast<real_t>(0.0)) {
    for (auto &[key, vals] : hist.partials) {
      vals = LDOSPostProcessor::broadenSpectrum(hist.bins, vals, params.broadening_kernel,
                                                params.broadening_sigma);
    }
  }

  if (params.rebin_num_bins > 0 && params.rebin_num_bins != hist.bins.size() &&
      hist.bins.size() >= 2) {
    const auto dst_bins =
        buildEnergyGrid(hist.bins.front(), hist.bins.back(), params.rebin_num_bins);
    for (auto &[key, vals] : hist.partials) {
      vals = LDOSPostProcessor::rebinSpectrum(hist.bins, vals, dst_bins);
    }
    hist.bins = dst_bins;
  }
}

[[nodiscard]] std::vector<real_t>
computePerAtomFermiDos(const correlation::core::Cell &cell,
                       const correlation::mlip::MLIPOutput &output,
                       std::span<const real_t> native_energies, const TDOSParams &params) {
  const size_t num_atoms = cell.atoms().size();
  std::vector<real_t> per_atom(num_atoms, static_cast<real_t>(0.0));
  for (size_t i = 0; i < num_atoms; ++i) {
    const auto atom_ldos = output.getLdos(i);
    if (atom_ldos.size() < native_energies.size()) {
      continue;
    }
    const auto spectrum_span = atom_ldos.subspan(0, native_energies.size());
    per_atom[i] = LDOSPostProcessor::evaluateFermiDos(native_energies, spectrum_span,
                                                      params.electronic_temp_k);
  }
  return per_atom;
}

void populateElectronicDiagnostics(const correlation::analysis::Histogram &tdos_hist,
                                   const TDOSParams &params, std::vector<real_t> per_atom_fermi,
                                   correlation::analysis::ElectronicDiagnostics &diag) {
  const auto it_total = tdos_hist.partials.find("total");
  if (tdos_hist.bins.empty() || it_total == tdos_hist.partials.end()) {
    return;
  }

  diag.has_electronic_dos = true;
  diag.fermi_dos_total = LDOSPostProcessor::evaluateFermiDos(tdos_hist.bins, it_total->second,
                                                             params.electronic_temp_k);
  const auto gap =
      LDOSPostProcessor::detectBandGap(tdos_hist.bins, it_total->second, params.gap_threshold);
  diag.band_gap = gap.band_gap;
  diag.vbm_energy = gap.vbm_energy;
  diag.cbm_energy = gap.cbm_energy;
  diag.is_metallic = gap.is_metallic;

  for (const auto &[key, vals] : tdos_hist.partials) {
    if (key == "total") {
      continue;
    }
    diag.fermi_dos_species[key] =
        LDOSPostProcessor::evaluateFermiDos(tdos_hist.bins, vals, params.electronic_temp_k);
  }
  diag.per_atom_fermi_dos = std::move(per_atom_fermi);
}

[[nodiscard]] TDOSParams
buildParamsFromSettings(const correlation::mlip::MLIPInterface *model,
                        const correlation::analysis::AnalysisSettings &settings) noexcept {
  TDOSParams params;
  params.model = model;
  params.broadening_kernel = settings.ldos_broadening_kernel;
  params.broadening_sigma = settings.ldos_broadening_sigma;
  params.electronic_temp_k = settings.ldos_electronic_temp_k;
  params.gap_threshold = settings.ldos_gap_threshold;
  params.compute_idos = settings.ldos_compute_idos;
  return params;
}

void storeResultInDistributions(correlation::analysis::DistributionFunctions &dists,
                                LDOSResult result) {
  if (!result.tdos.bins.empty()) {
    dists.addHistogram("TDOS", std::move(result.tdos));
  }
  if (!result.idos.bins.empty()) {
    dists.addHistogram("IDOS", std::move(result.idos));
  }
  if (result.diagnostics.has_electronic_dos) {
    dists.setElectronicDiagnostics(std::move(result.diagnostics));
  }
}

} // namespace

void TDOSCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                    const correlation::analysis::AnalysisSettings &settings) const {
  if (model_ == nullptr) {
    return;
  }
  const auto params = buildParamsFromSettings(model_, settings);
  storeResultInDistributions(dists, calculateDetailed(dists.cell(), params));
}

void TDOSCalculator::calculateTrajectory(
    correlation::analysis::DistributionFunctions &dists, const correlation::core::Trajectory &traj,
    const correlation::analysis::AnalysisSettings &settings) const {
  if (model_ == nullptr) {
    return;
  }
  const auto params = buildParamsFromSettings(model_, settings);
  storeResultInDistributions(dists,
                             calculateTrajectoryDetailed(traj, params, settings.cancel_flag));
}

LDOSResult TDOSCalculator::calculateDetailed(const correlation::core::Cell &cell,
                                             const TDOSParams &params) {
  LDOSResult result;
  if (params.model == nullptr || cell.atoms().empty()) {
    return result;
  }

  const auto output = params.model->evaluate(cell);
  const size_t num_bins = output.ldos_bins;
  if (num_bins == 0 || output.ldos.empty()) {
    return result;
  }

  real_t grid_min = params.e_min;
  real_t grid_max = params.e_max;
  resolveEnergyBounds(output, params, grid_min, grid_max);

  auto native_bins = buildEnergyGrid(grid_min, grid_max, num_bins);
  auto per_atom_fermi = computePerAtomFermiDos(cell, output, native_bins, params);

  populateHistogramMetadata(result.tdos, std::move(native_bins));
  result.tdos.compute_count = 1;

  std::vector<double> acc_total(num_bins, 0.0);
  std::map<std::string, std::vector<double>> acc_species;
  accumulateFrameLdos(cell, output, num_bins, acc_total, acc_species);
  finalizeRawHistogram(result.tdos, 1.0, acc_total, acc_species);
  applyPostProcessing(result.tdos, params);

  populateElectronicDiagnostics(result.tdos, params, std::move(per_atom_fermi), result.diagnostics);
  if (params.compute_idos) {
    populateIdosHistogram(result.tdos, result.idos);
  }
  if (params.retain_per_atom_ldos) {
    result.per_atom_ldos = output.ldos;
    result.ldos_bins = num_bins;
  }
  return result;
}

correlation::analysis::Histogram TDOSCalculator::calculate(const correlation::core::Cell &cell,
                                                           const TDOSParams &params) {
  return calculateDetailed(cell, params).tdos;
}

LDOSResult TDOSCalculator::calculateTrajectoryDetailed(const correlation::core::Trajectory &traj,
                                                       const TDOSParams &params,
                                                       const std::atomic<bool> *cancel_flag) {
  LDOSResult result;
  const size_t total_frames = traj.getFrameCount();
  if (params.model == nullptr || total_frames == 0) {
    return result;
  }

  size_t num_bins = 0;
  real_t grid_min = params.e_min;
  real_t grid_max = params.e_max;
  std::vector<real_t> native_bins;
  std::vector<double> acc_total;
  std::map<std::string, std::vector<double>> acc_species;
  std::vector<double> acc_per_atom_fermi;
  size_t processed_frames = 0;

  for (size_t frame_idx = 0; frame_idx < total_frames; ++frame_idx) {
    if (cancel_flag != nullptr && cancel_flag->load()) {
      break;
    }

    const auto &cell = traj.getFrame(frame_idx);
    if (cell.atoms().empty()) {
      continue;
    }

    const auto output = params.model->evaluate(cell);
    if (output.ldos_bins == 0 || output.ldos.empty()) {
      continue;
    }

    if (num_bins == 0) {
      num_bins = output.ldos_bins;
      acc_total.assign(num_bins, 0.0);
      resolveEnergyBounds(output, params, grid_min, grid_max);
      native_bins = buildEnergyGrid(grid_min, grid_max, num_bins);
      acc_per_atom_fermi.assign(cell.atoms().size(), 0.0);
    }

    accumulateFrameLdos(cell, output, num_bins, acc_total, acc_species);
    const auto frame_fermi = computePerAtomFermiDos(cell, output, native_bins, params);
    if (frame_fermi.size() == acc_per_atom_fermi.size()) {
      for (size_t atom_idx = 0; atom_idx < frame_fermi.size(); ++atom_idx) {
        acc_per_atom_fermi[atom_idx] += static_cast<double>(frame_fermi[atom_idx]);
      }
    }
    ++processed_frames;
  }

  if (processed_frames == 0 || num_bins == 0) {
    return result;
  }

  populateHistogramMetadata(result.tdos, std::move(native_bins));
  result.tdos.compute_count = static_cast<int>(processed_frames);

  const double frame_scale = 1.0 / static_cast<double>(processed_frames);
  finalizeRawHistogram(result.tdos, frame_scale, acc_total, acc_species);
  applyPostProcessing(result.tdos, params);

  std::vector<real_t> avg_per_atom_fermi(acc_per_atom_fermi.size(), static_cast<real_t>(0.0));
  for (size_t i = 0; i < acc_per_atom_fermi.size(); ++i) {
    avg_per_atom_fermi[i] = static_cast<real_t>(acc_per_atom_fermi[i] * frame_scale);
  }
  populateElectronicDiagnostics(result.tdos, params, std::move(avg_per_atom_fermi),
                                result.diagnostics);
  if (params.compute_idos) {
    populateIdosHistogram(result.tdos, result.idos);
  }
  return result;
}

correlation::analysis::Histogram
TDOSCalculator::calculateTrajectory(const correlation::core::Trajectory &traj,
                                    const TDOSParams &params,
                                    const std::atomic<bool> *cancel_flag) {
  return calculateTrajectoryDetailed(traj, params, cancel_flag).tdos;
}

} // namespace correlation::calculators
