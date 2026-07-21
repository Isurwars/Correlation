/**
 * @file VDOSCalculator.cpp
 * @brief Implementation of vibrational density of states calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/VDOSCalculator.hpp"
#include "analysis/DynamicsAnalyzer.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"

#include <stdexcept>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
// NOLINTNEXTLINE(cert-err58-cpp)
const bool registered = CalculatorFactory::registerTypeSafe<VDOSCalculator>("VDOSCalculator");
} // namespace

void VDOSCalculator::calculateTrajectory(correlation::analysis::DistributionFunctions &dists,
                                         const correlation::core::Trajectory & /*traj*/,
                                         const correlation::analysis::AnalysisSettings & /*settings*/) const {
  const auto &all = dists.getAllHistograms();
  if (!all.contains("VACF")) {
    return; // VACF must be computed first
  }
  dists.addHistogram("VDOS", calculate(dists.getHistogram("VACF")));
}

correlation::analysis::Histogram VDOSCalculator::calculate(const correlation::analysis::Histogram &vacf_hist,
                                                           const VDOSParams &params) {
  const auto &vacf_data = vacf_hist.partials.at("Total");

  if (vacf_data.size() < 2) {
    throw std::logic_error("VACF data is too short for VDOS calculation.");
  }

  real_t const time_step = vacf_hist.bins[1] - vacf_hist.bins[0];

  auto [frequencies, intensities_real, intensities_imag] =
      correlation::analysis::DynamicsAnalyzer::calculateVDOS(vacf_data, time_step);

  size_t const num_points = frequencies.size();

  std::vector<real_t> combined_frequencies;
  std::vector<real_t> combined_frequencies_cmInv;
  std::vector<real_t> combined_frequencies_meV;
  std::vector<real_t> combined_intensities;
  combined_frequencies.reserve(2 * num_points);
  combined_frequencies_cmInv.reserve(2 * num_points);
  combined_frequencies_meV.reserve(2 * num_points);
  combined_intensities.reserve(2 * num_points);

  // Imaginary frequencies (negative axis, range: -params.max_imag_freq to 0)
  for (size_t i = num_points - 1; i > 0; --i) {
    if (frequencies[i] <= params.max_imag_freq) {
      combined_frequencies.push_back(-frequencies[i]);
      combined_frequencies_cmInv.push_back(-frequencies[i] * static_cast<real_t>(correlation::math::thz_to_cminv));
      combined_frequencies_meV.push_back(-frequencies[i] * static_cast<real_t>(correlation::math::thz_to_mev));
      combined_intensities.push_back(intensities_imag[i]);
    }
  }

  // Real frequencies (positive axis, range: 0 to +params.max_real_freq)
  for (size_t i = 0; i < num_points; ++i) {
    if (frequencies[i] <= params.max_real_freq) {
      combined_frequencies.push_back(frequencies[i]);
      combined_frequencies_cmInv.push_back(frequencies[i] * static_cast<real_t>(correlation::math::thz_to_cminv));
      combined_frequencies_meV.push_back(frequencies[i] * static_cast<real_t>(correlation::math::thz_to_mev));
      combined_intensities.push_back(intensities_real[i]);
    }
  }

  correlation::analysis::Histogram vdos_hist;
  vdos_hist.x_label = "ν";
  vdos_hist.title = "Vibrational Density of States";
  vdos_hist.description = "Vibrational Density of States";
  vdos_hist.file_suffix = "_VDOS";
  vdos_hist.y_label = "f(ω)";
  vdos_hist.x_unit = "THz";
  vdos_hist.y_unit = "arbitrary units";
  vdos_hist.bins = combined_frequencies;
  vdos_hist.partials["Total"] = combined_intensities;
  vdos_hist.partials["Frequency_cm_1"] = combined_frequencies_cmInv;
  vdos_hist.partials["Frequency_meV"] = combined_frequencies_meV;

  return vdos_hist;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
correlation::analysis::Histogram VDOSCalculator::calculate(const correlation::analysis::Histogram &vacf_hist,
                                                           real_t max_imag_freq, real_t max_real_freq) {
  return calculate(vacf_hist, VDOSParams{.max_imag_freq = max_imag_freq, .max_real_freq = max_real_freq});
}

} // namespace correlation::calculators
