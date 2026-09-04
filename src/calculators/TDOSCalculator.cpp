/**
 * @file TDOSCalculator.cpp
 * @brief Implementation of Total Density of States (TDOS) calculator from MLIP predictions.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/TDOSCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <map>
#include <span>
#include <string>
#include <vector>

namespace correlation::calculators {

namespace {

const bool registered = CalculatorFactory::registerTypeSafe<TDOSCalculator>("TDOSCalculator");

/**
 * @brief Generates uniform energy bin centers between e_min and e_max.
 */
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

/**
 * @brief Populates standard metadata fields on a TDOS Histogram.
 */
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

/**
 * @brief Accumulates per-atom LDOS into total and species-resolved double-precision accumulators.
 */
void accumulateFrameLdos(const correlation::core::Cell &cell, const std::vector<std::vector<real_t>> &ldos,
                         size_t num_bins, std::vector<double> &acc_total,
                         std::map<std::string, std::vector<double>> &acc_species) {
  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();
  if (num_atoms == 0 || ldos.size() < num_atoms) {
    return;
  }

  const double norm_factor = 1.0 / static_cast<double>(num_atoms);

  for (size_t i = 0; i < num_atoms; ++i) {
    const auto &atom_ldos = ldos[i];
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

/**
 * @brief Finalizes accumulated double-precision sums into the target Histogram.
 */
void finalizeHistogram(correlation::analysis::Histogram &hist, double frame_scale,
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

} // namespace

void TDOSCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                    const correlation::analysis::AnalysisSettings & /*settings*/) const {
  const auto &cell = dists.cell();
  const TDOSParams params;
  auto hist = calculate(cell, params);
  if (!hist.bins.empty()) {
    dists.addHistogram("tdos", std::move(hist));
  }
}

void TDOSCalculator::calculateTrajectory(correlation::analysis::DistributionFunctions &dists,
                                         const correlation::core::Trajectory &traj,
                                         const correlation::analysis::AnalysisSettings &settings) const {
  const TDOSParams params;
  auto hist = calculateTrajectory(traj, params, settings.cancel_flag);
  if (!hist.bins.empty()) {
    dists.addHistogram("tdos", std::move(hist));
  }
}

correlation::analysis::Histogram TDOSCalculator::calculate(const correlation::core::Cell &cell,
                                                           const TDOSParams &params) {
  correlation::analysis::Histogram hist;
  if (params.model == nullptr || cell.atoms().empty()) {
    return hist;
  }

  const auto output = params.model->evaluate(cell);
  const size_t num_bins = output.ldos_bins;
  if (num_bins == 0 || output.ldos.empty()) {
    return hist;
  }

  populateHistogramMetadata(hist, buildEnergyGrid(params.e_min, params.e_max, num_bins));
  hist.compute_count = 1;

  std::vector<double> acc_total(num_bins, 0.0);
  std::map<std::string, std::vector<double>> acc_species;

  accumulateFrameLdos(cell, output.ldos, num_bins, acc_total, acc_species);
  finalizeHistogram(hist, 1.0, acc_total, acc_species);

  return hist;
}

correlation::analysis::Histogram TDOSCalculator::calculateTrajectory(const correlation::core::Trajectory &traj,
                                                                     const TDOSParams &params,
                                                                     const std::atomic<bool> *cancel_flag) {
  correlation::analysis::Histogram hist;
  const size_t total_frames = traj.getFrameCount();
  if (params.model == nullptr || total_frames == 0) {
    return hist;
  }

  size_t num_bins = 0;
  std::vector<double> acc_total;
  std::map<std::string, std::vector<double>> acc_species;
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
    }

    accumulateFrameLdos(cell, output.ldos, num_bins, acc_total, acc_species);
    ++processed_frames;
  }

  if (processed_frames == 0 || num_bins == 0) {
    return hist;
  }

  populateHistogramMetadata(hist, buildEnergyGrid(params.e_min, params.e_max, num_bins));
  hist.compute_count = static_cast<int>(processed_frames);

  const double frame_scale = 1.0 / static_cast<double>(processed_frames);
  finalizeHistogram(hist, frame_scale, acc_total, acc_species);

  return hist;
}

} // namespace correlation::calculators
