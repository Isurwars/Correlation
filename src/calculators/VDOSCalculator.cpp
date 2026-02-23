// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/VDOSCalculator.hpp"
#include "DynamicsAnalyzer.hpp"
#include "PhysicalData.hpp"
#include <stdexcept>

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
    combined_frequencies_cmInv.push_back(-frequencies[i] *
                                         constants::THz_to_cmInv);
    combined_frequencies_meV.push_back(-frequencies[i] * constants::THz_to_meV);
    combined_intensities.push_back(intensities_imag[i]);
  }

  for (size_t i = 0; i < num_points; ++i) {
    combined_frequencies.push_back(frequencies[i]);
    combined_frequencies_cmInv.push_back(frequencies[i] *
                                         constants::THz_to_cmInv);
    combined_frequencies_meV.push_back(frequencies[i] * constants::THz_to_meV);
    combined_intensities.push_back(intensities_real[i]);
  }

  Histogram vdos_hist;
  vdos_hist.bin_label = "Frequency (THz)";
  vdos_hist.bins = combined_frequencies;
  vdos_hist.partials["Total"] = combined_intensities;
  vdos_hist.partials["Frequency (cm-1)"] = combined_frequencies_cmInv;
  vdos_hist.partials["Frequency (meV)"] = combined_frequencies_meV;

  return vdos_hist;
}
