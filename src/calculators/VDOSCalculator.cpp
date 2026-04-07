/**
 * @file VDOSCalculator.cpp
 * @brief Implementation of vibrational density of states calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "calculators/VDOSCalculator.hpp"
#include "DynamicsAnalyzer.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "math/Constants.hpp"

#include <stdexcept>

namespace {
bool registered = CalculatorFactory::instance().registerCalculator(
    std::make_unique<VDOSCalculator>());
} // namespace

void VDOSCalculator::calculateTrajectory(
    DistributionFunctions &df, const Trajectory &traj,
    const AnalysisSettings &settings) const {
  const auto &all = df.getAllHistograms();
  if (all.find("VACF") == all.end()) {
    return; // VACF must be computed first
  }
  df.addHistogram("VDOS", calculate(df.getHistogram("VACF")));
}

Histogram VDOSCalculator::calculate(const Histogram &vacf_hist) {
  const auto &vacf_data = vacf_hist.partials.at("Total");

  if (vacf_data.size() < 2) {
    throw std::logic_error("VACF data is too short for VDOS calculation.");
  }

  double dt = vacf_hist.bins[1] - vacf_hist.bins[0];

  auto [frequencies, intensities_real, intensities_imag] =
      DynamicsAnalyzer::calculateVDOS(vacf_data, dt);

  size_t num_points = frequencies.size();
  size_t total_points = 2 * num_points - 1;

  std::vector<double> combined_frequencies;
  std::vector<double> combined_frequencies_cmInv;
  std::vector<double> combined_frequencies_meV;
  std::vector<double> combined_intensities;
  combined_frequencies.reserve(total_points);
  combined_frequencies_cmInv.reserve(total_points);
  combined_frequencies_meV.reserve(total_points);
  combined_intensities.reserve(total_points);

  for (size_t i = num_points - 1; i > 0; --i) {
    combined_frequencies.push_back(-frequencies[i]);
    combined_frequencies_cmInv.push_back(
        -frequencies[i] * correlation::math::thz_to_cminv);
    combined_frequencies_meV.push_back(
        -frequencies[i] * correlation::math::thz_to_mev);
    combined_intensities.push_back(intensities_imag[i]);
  }

  for (size_t i = 0; i < num_points; ++i) {
    combined_frequencies.push_back(frequencies[i]);
    combined_frequencies_cmInv.push_back(
        frequencies[i] * correlation::math::thz_to_cminv);
    combined_frequencies_meV.push_back(
        frequencies[i] * correlation::math::thz_to_mev);
    combined_intensities.push_back(intensities_real[i]);
  }

  Histogram vdos_hist;
  vdos_hist.bin_label = "Frequency_THz";
  vdos_hist.bins = combined_frequencies;
  vdos_hist.partials["Total"] = combined_intensities;
  vdos_hist.partials["Frequency_cm_1"] = combined_frequencies_cmInv;
  vdos_hist.partials["Frequency_meV"] = combined_frequencies_meV;

  return vdos_hist;
}
