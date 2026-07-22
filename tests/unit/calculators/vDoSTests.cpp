// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DynamicsAnalyzer.hpp"
#include "math/Constants.hpp"
#include "math/Precision.hpp"

#include <algorithm>
#include <gtest/gtest.h>
#include <vector>

namespace correlation::analysis {

TEST(vDoSTests, VDOSIsNonZeroAtZeroFrequencyForConstantVACF) {
  // Create a simple constant VACF (DC signal)
  // Fourier transform of a constant is a delta at f=0.
  // Since we define VDOS as the Fourier Transform of the VACF (without
  // frequency weighting), a constant VACF should result in a significant peak
  // at f=0 (DC component).

  std::vector<real_t> vacf(100, 1.0); // Constant VACF
  real_t time_step = 1.0;

  auto [frequencies, intensities, unused] = DynamicsAnalyzer::calculateVDOS(vacf, time_step);

  ASSERT_FALSE(frequencies.empty());
  ASSERT_EQ(frequencies[0], 0.0);

  // The integral of a constant 1.0 over 100 points with dt=1.0 is roughly
  // 100.0. With windowing and normalization, it might vary, but it should
  // definitely be non-zero.
  EXPECT_GT(intensities[0], 1.0) << "VDOS at 0 THz should be non-zero for a constant VACF (DC component)";
}

TEST(vDoSTests, PerfectSolidShowsSinglePeak) {
  // A perfect solid has a VACF represented by non-decaying cosine waves.
  // We model a single vibrational mode as a pure cosine wave.
  // The resulting VDOS (real part) should have a single maximum at that
  // specific frequency.

  real_t time_step = 1.0; // 1 fs
  size_t num_frames = 2000;
  std::vector<real_t> vacf(num_frames);

  // Frequency of the pure cosine wave in THz
  real_t target_nu_thz = 25.0;

  for (size_t i = 0; i < num_frames; ++i) {
    real_t t_fs = static_cast<real_t>(i) * time_step;
    // 0.001 converts fs to ps for THz frequency
    vacf[i] = static_cast<real_t>(std::cos(2.0 * correlation::math::pi * target_nu_thz * t_fs * 0.001));
  }

  auto [frequencies, intensities_real, intensities_imag] = DynamicsAnalyzer::calculateVDOS(vacf, time_step);

  // Find the maximum intensity in the real part
  auto max_it = std::ranges::max_element(intensities_real);
  size_t max_idx = std::distance(intensities_real.begin(), max_it);

  // Delta nu is the frequency resolution
  real_t d_nu = frequencies[1] - frequencies[0];

  // Check if the maximum frequency matches the target frequency within an
  // acceptable resolution margin
  EXPECT_NEAR(frequencies[max_idx], target_nu_thz, 2.0 * d_nu)
      << "The calculated VDOS max frequency " << frequencies[max_idx]
      << " THz does not correspond to the input VACF frequency " << target_nu_thz << " THz.";

  EXPECT_GT(*max_it, 0.0) << "The peak real intensity should be positive.";
}

TEST(vDoSTests, IdealGasShowsImaginaryPeak) {
  // An ideal gas (or a purely diffusive liquid) has a VACF that decays
  // exponentially without oscillations. A pure exponential decay e^(-at)
  // translates to a VDOS where the real part peaks at 0 THz and the imaginary
  // part peaks at frequency nu = a / 2pi.

  real_t time_step = 1.0; // 1 fs
  size_t num_frames = 2000;
  std::vector<real_t> vacf(num_frames);

  // We choose an exponential decay coefficient 'a' such that the imaginary peak
  // is at 10.0 THz.
  real_t target_nu_peak = 10.0;
  auto vec_a = static_cast<real_t>(2.0 * correlation::math::pi * target_nu_peak); // a in THz (ps^-1)

  for (size_t i = 0; i < num_frames; ++i) {
    auto t_ps = static_cast<real_t>(static_cast<real_t>(i) * time_step * 0.001); // fs to ps
    vacf[i] = std::exp(-vec_a * t_ps);
  }

  auto [frequencies, intensities_real, intensities_imag] = DynamicsAnalyzer::calculateVDOS(vacf, time_step);

  // 1. The Real part should peak at exactly 0 Hz for pure exponential decay
  auto max_real_it = std::ranges::max_element(intensities_real);
  size_t max_real_idx = std::distance(intensities_real.begin(), max_real_it);

  EXPECT_NEAR(frequencies[max_real_idx], 0.0, correlation::is_single_precision ? 0.5 : 1e-6)
      << "The real VDOS for an exponential decay should have its maximum at 0 THz.";

  // 2. The Imaginary part should peak at target_nu_peak
  auto max_imag_it = std::ranges::max_element(intensities_imag);
  size_t max_imag_idx = std::distance(intensities_imag.begin(), max_imag_it);

  real_t d_nu = frequencies[1] - frequencies[0];

  EXPECT_NEAR(frequencies[max_imag_idx], target_nu_peak, 2.0 * d_nu)
      << "The imaginary VDOS peak " << frequencies[max_imag_idx]
      << " THz does not correspond to the decay constant frequency " << target_nu_peak << " THz.";
}
} // namespace correlation::analysis
