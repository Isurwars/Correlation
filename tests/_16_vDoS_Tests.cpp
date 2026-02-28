// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>
#include <vector>

#include "../include/DynamicsAnalyzer.hpp"
#include "../include/PhysicalData.hpp"

TEST(_16_vDoS_Tests, VDOSIsNonZeroAtZeroFrequencyForConstantVACF) {
  // Create a simple constant VACF (DC signal)
  // Fourier transform of a constant is a delta at f=0.
  // Since we define VDOS as the Fourier Transform of the VACF (without
  // frequency weighting), a constant VACF should result in a significant peak
  // at f=0 (DC component).

  std::vector<double> vacf(100, 1.0); // Constant VACF
  double dt = 1.0;

  auto [frequencies, intensities, _] =
      DynamicsAnalyzer::calculateVDOS(vacf, dt);

  ASSERT_FALSE(frequencies.empty());
  ASSERT_EQ(frequencies[0], 0.0);

  // The integral of a constant 1.0 over 100 points with dt=1.0 is roughly
  // 100.0. With windowing and normalization, it might vary, but it should
  // definitely be non-zero.
  EXPECT_GT(intensities[0], 1.0)
      << "VDOS at 0 THz should be non-zero for a constant VACF (DC component)";
}

TEST(_16_vDoS_Tests, PerfectSolidShowsSinglePeak) {
  // A perfect solid has a VACF represented by non-decaying cosine waves.
  // We model a single vibrational mode as a pure cosine wave.
  // The resulting VDOS (real part) should have a single maximum at that
  // specific frequency.

  double dt = 1.0; // 1 fs
  size_t num_frames = 2000;
  std::vector<double> vacf(num_frames);

  // Frequency of the pure cosine wave in THz
  double target_nu_thz = 25.0;

  for (size_t i = 0; i < num_frames; ++i) {
    double t_fs = i * dt;
    // 0.001 converts fs to ps for THz frequency
    vacf[i] = std::cos(2.0 * constants::pi * target_nu_thz * t_fs * 0.001);
  }

  auto [frequencies, intensities_real, intensities_imag] =
      DynamicsAnalyzer::calculateVDOS(vacf, dt);

  // Find the maximum intensity in the real part
  auto max_it =
      std::max_element(intensities_real.begin(), intensities_real.end());
  size_t max_idx = std::distance(intensities_real.begin(), max_it);

  // Delta nu is the frequency resolution
  double d_nu = frequencies[1] - frequencies[0];

  // Check if the maximum frequency matches the target frequency within an
  // acceptable resolution margin
  EXPECT_NEAR(frequencies[max_idx], target_nu_thz, 2.0 * d_nu)
      << "The calculated VDOS max frequency " << frequencies[max_idx]
      << " THz does not correspond to the input VACF frequency "
      << target_nu_thz << " THz.";

  EXPECT_GT(*max_it, 0.0) << "The peak real intensity should be positive.";
}

TEST(_16_vDoS_Tests, IdealGasShowsImaginaryPeak) {
  // An ideal gas (or a purely diffusive liquid) has a VACF that decays
  // exponentially without oscillations. A pure exponential decay e^(-at)
  // translates to a VDOS where the real part peaks at 0 THz and the imaginary
  // part peaks at frequency nu = a / 2pi.

  double dt = 1.0; // 1 fs
  size_t num_frames = 2000;
  std::vector<double> vacf(num_frames);

  // We choose an exponential decay coefficient 'a' such that the imaginary peak
  // is at 10.0 THz.
  double target_nu_peak = 10.0;
  double a = 2.0 * constants::pi * target_nu_peak; // a in THz (ps^-1)

  for (size_t i = 0; i < num_frames; ++i) {
    double t_ps = i * dt * 0.001; // fs to ps
    vacf[i] = std::exp(-a * t_ps);
  }

  auto [frequencies, intensities_real, intensities_imag] =
      DynamicsAnalyzer::calculateVDOS(vacf, dt);

  // 1. The Real part should peak at exactly 0 Hz for pure exponential decay
  auto max_real_it =
      std::max_element(intensities_real.begin(), intensities_real.end());
  size_t max_real_idx = std::distance(intensities_real.begin(), max_real_it);

  EXPECT_EQ(frequencies[max_real_idx], 0.0)
      << "The real VDOS for an exponential decay should have its maximum at 0 "
         "THz.";

  // 2. The Imaginary part should peak at target_nu_peak
  auto max_imag_it =
      std::max_element(intensities_imag.begin(), intensities_imag.end());
  size_t max_imag_idx = std::distance(intensities_imag.begin(), max_imag_it);

  double d_nu = frequencies[1] - frequencies[0];

  EXPECT_NEAR(frequencies[max_imag_idx], target_nu_peak, 2.0 * d_nu)
      << "The imaginary VDOS peak " << frequencies[max_imag_idx]
      << " THz does not correspond to the decay constant frequency "
      << target_nu_peak << " THz.";
}
