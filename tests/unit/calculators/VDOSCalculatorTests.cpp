// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DynamicsAnalyzer.hpp"
#include "calculators/temporal/VDOSCalculator.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/Constants.hpp"
#include "math/Precision.hpp"

#include <algorithm>
#include <gtest/gtest.h>
#include <vector>

using correlation::real_t;

namespace correlation::analysis {
namespace {

TEST(VDOSCalculatorTests, VDOSIsNonZeroAtZeroFrequencyForConstantVACF) {
  // Create a simple constant VACF (DC signal)
  // Fourier transform of a constant is a delta at f=0.
  // Since we define VDOS as the Fourier Transform of the VACF (without
  // frequency weighting), a constant VACF should result in a significant peak
  // at f=0 (DC component).

  const std::vector<real_t> vacf(100, 1.0); // Constant VACF
  const real_t time_step = 1.0;

  const auto [frequencies, intensities, unused] = DynamicsAnalyzer::calculateVDOS(vacf, time_step);

  ASSERT_FALSE(frequencies.empty());
  ASSERT_EQ(frequencies[0], 0.0);

  // The integral of a constant 1.0 over 100 points with dt=1.0 is roughly
  // 100.0. With windowing and normalization, it might vary, but it should
  // definitely be non-zero.
  EXPECT_GT(intensities[0], 1.0)
      << "VDOS at 0 THz should be non-zero for a constant VACF (DC component)";
}

TEST(VDOSCalculatorTests, PerfectSolidShowsSinglePeak) {
  // A perfect solid has a VACF represented by non-decaying cosine waves.
  // We model a single vibrational mode as a pure cosine wave.
  // The resulting VDOS (real part) should have a single maximum at that
  // specific frequency.

  const real_t time_step = 1.0; // 1 fs
  const size_t num_frames = 2000;
  std::vector<real_t> vacf(num_frames);

  // Frequency of the pure cosine wave in THz
  const real_t target_nu_thz = 25.0;

  for (size_t i = 0; i < num_frames; ++i) {
    const real_t t_fs = static_cast<real_t>(i) * time_step;
    // 0.001 converts fs to ps for THz frequency
    vacf[i] =
        static_cast<real_t>(std::cos(2.0 * correlation::math::PI * target_nu_thz * t_fs * 0.001));
  }

  const auto [frequencies, intensities_real, intensities_imag] =
      DynamicsAnalyzer::calculateVDOS(vacf, time_step);

  // Find the maximum intensity in the real part
  const auto max_it = std::ranges::max_element(intensities_real);
  const size_t max_idx = std::distance(intensities_real.begin(), max_it);

  // Delta nu is the frequency resolution
  const real_t d_nu = frequencies[1] - frequencies[0];

  // Check if the maximum frequency matches the target frequency within an
  // acceptable resolution margin
  EXPECT_NEAR(frequencies[max_idx], target_nu_thz, 2.0 * d_nu)
      << "The calculated VDOS max frequency " << frequencies[max_idx]
      << " THz does not correspond to the input VACF frequency " << target_nu_thz << " THz.";

  EXPECT_GT(*max_it, 0.0) << "The peak real intensity should be positive.";
}

TEST(VDOSCalculatorTests, IdealGasShowsImaginaryPeak) {
  // An ideal gas (or a purely diffusive liquid) has a VACF that decays
  // exponentially without oscillations. A pure exponential decay e^(-at)
  // translates to a VDOS where the real part peaks at 0 THz and the imaginary
  // part peaks at frequency nu = a / 2pi.

  const real_t time_step = 1.0; // 1 fs
  const size_t num_frames = 2000;
  std::vector<real_t> vacf(num_frames);

  // We choose an exponential decay coefficient 'a' such that the imaginary peak
  // is at 10.0 THz.
  const real_t target_nu_peak = 10.0;
  const auto vec_a =
      static_cast<real_t>(2.0 * correlation::math::PI * target_nu_peak); // a in THz (ps^-1)

  for (size_t i = 0; i < num_frames; ++i) {
    const auto t_ps = static_cast<real_t>(static_cast<real_t>(i) * time_step * 0.001); // fs to ps
    vacf[i] = std::exp(-vec_a * t_ps);
  }

  const auto [frequencies, intensities_real, intensities_imag] =
      DynamicsAnalyzer::calculateVDOS(vacf, time_step);

  // 1. The Real part should peak at exactly 0 Hz for pure exponential decay
  const auto max_real_it = std::ranges::max_element(intensities_real);
  const size_t max_real_idx = std::distance(intensities_real.begin(), max_real_it);

  EXPECT_NEAR(frequencies[max_real_idx], 0.0, correlation::is_single_precision ? 0.5 : 1e-6)
      << "The real VDOS for an exponential decay should have its maximum at 0 THz.";

  // 2. The Imaginary part should peak at target_nu_peak
  const auto max_imag_it = std::ranges::max_element(intensities_imag);
  const size_t max_imag_idx = std::distance(intensities_imag.begin(), max_imag_it);

  const real_t d_nu = frequencies[1] - frequencies[0];

  EXPECT_NEAR(frequencies[max_imag_idx], target_nu_peak, 2.0 * d_nu)
      << "The imaginary VDOS peak " << frequencies[max_imag_idx]
      << " THz does not correspond to the decay constant frequency " << target_nu_peak << " THz.";
}

TEST(VDOSCalculatorTests, ComputesRecommendedMaxFrequencyFromLightestAtom) {
  using correlation::calculators::VDOSCalculator;
  using correlation::core::Cell;

  // 1. Hydrogen cell (m ~ 1.008 amu) -> ν_max ~ 50.0 / sqrt(1.008) ~ 49.8 THz
  Cell h_cell;
  h_cell.addAtom("H", {0.0, 0.0, 0.0});
  const real_t h_max_freq = VDOSCalculator::computeRecommendedMaxFrequency(h_cell, 1.0);
  EXPECT_NEAR(h_max_freq, 49.8, 0.5);

  // 2. Carbon cell (m ~ 12.011 amu) -> ν_max ~ 50.0 / sqrt(12.011) ~ 14.4 THz
  Cell c_cell;
  c_cell.addAtom("C", {0.0, 0.0, 0.0});
  const real_t c_max_freq = VDOSCalculator::computeRecommendedMaxFrequency(c_cell, 1.0);
  EXPECT_NEAR(c_max_freq, 14.4, 0.5);

  // 3. Silicon cell (m ~ 28.085 amu) -> scaled ~ 9.4 THz -> clamped to 10.0 THz minimum
  Cell si_cell;
  si_cell.addAtom("Si", {0.0, 0.0, 0.0});
  const real_t si_max_freq = VDOSCalculator::computeRecommendedMaxFrequency(si_cell, 1.0);
  EXPECT_DOUBLE_EQ(si_max_freq, 10.0);

  // 4. Mixed cell with Si and H -> lightest atom (H) determines frequency ~ 49.8 THz
  Cell mixed_cell;
  mixed_cell.addAtom("Si", {0.0, 0.0, 0.0});
  mixed_cell.addAtom("H", {1.0, 0.0, 0.0});
  const real_t mixed_max_freq = VDOSCalculator::computeRecommendedMaxFrequency(mixed_cell, 1.0);
  EXPECT_NEAR(mixed_max_freq, 49.8, 0.5);

  // 5. Large time step (50 fs) -> Nyquist bound = 1000 / (2 * 50) = 10.0 THz
  const real_t nyquist_bounded_freq = VDOSCalculator::computeRecommendedMaxFrequency(h_cell, 50.0);
  EXPECT_DOUBLE_EQ(nyquist_bounded_freq, 10.0);

  // 6. Empty cell defaults to 10.0 THz
  Cell empty_cell;
  const real_t empty_max_freq = VDOSCalculator::computeRecommendedMaxFrequency(empty_cell, 1.0);
  EXPECT_DOUBLE_EQ(empty_max_freq, 10.0);
}

TEST(VDOSCalculatorTests, CalculateTrajectoryCalculatesHydrogenFrequenciesUpTo50THz) {
  using correlation::analysis::DistributionFunctions;
  using correlation::calculators::VDOSCalculator;
  using correlation::core::Cell;
  using correlation::core::Trajectory;

  Cell h_cell;
  h_cell.addAtom("H", {0.0, 0.0, 0.0});
  Trajectory traj({h_cell, h_cell}, 1.0);

  DistributionFunctions dists(h_cell);
  // Synthetic VACF with 50 frames
  Histogram vacf_hist;
  vacf_hist.bins.resize(50);
  for (size_t i = 0; i < 50; ++i) {
    vacf_hist.bins[i] = static_cast<real_t>(i); // dt = 1 fs
  }
  std::vector<real_t> vacf_data(50, 1.0);
  vacf_hist.partials["Total"] = vacf_data;
  dists.addHistogram("VACF", std::move(vacf_hist));

  VDOSCalculator calc;
  AnalysisSettings settings;
  calc.calculateTrajectory(dists, traj, settings);

  ASSERT_TRUE(dists.hasHistogram("VDOS"));
  const auto &vdos_hist = dists.getHistogram("VDOS");
  ASSERT_FALSE(vdos_hist.bins.empty());

  const real_t max_computed_freq = vdos_hist.bins.back();
  EXPECT_GT(max_computed_freq, 45.0)
      << "Hydrogen VDOS frequency should reach ~50 THz instead of being truncated to 10 THz";
  EXPECT_LE(max_computed_freq, 50.0);
}

} // namespace
} // namespace correlation::analysis
