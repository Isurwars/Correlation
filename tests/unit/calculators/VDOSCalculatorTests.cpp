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
  const std::vector<real_t> vacf(100, 1.0); // Constant VACF (pure DC / diffusive)
  const real_t time_step = 1.0;

  const auto [frequencies, intensities] = DynamicsAnalyzer::calculateVDOS(vacf, time_step);

  ASSERT_FALSE(frequencies.empty());
  ASSERT_EQ(frequencies[0], 0.0);
  EXPECT_GT(intensities[0], 1.0)
      << "VDOS at 0 THz should be non-zero for a constant VACF (DC component)";

  // Smooth monotonic decay near ν = 0 (no Filon spike / cliff at freq_idx = 1)
  EXPECT_GE(intensities[0], intensities[1]);
  for (const real_t val : intensities) {
    EXPECT_GE(val, static_cast<real_t>(0.0));
  }
}

TEST(VDOSCalculatorTests, PerfectSolidShowsSinglePeakAndNoZeroJump) {
  const real_t time_step = 1.0; // 1 fs
  const size_t num_frames = 2000;
  std::vector<real_t> vacf(num_frames);

  const real_t target_nu_thz = 25.0;
  for (size_t i = 0; i < num_frames; ++i) {
    const real_t t_fs = static_cast<real_t>(i) * time_step;
    vacf[i] =
        static_cast<real_t>(std::cos(2.0 * correlation::math::PI * target_nu_thz * t_fs * 0.001));
  }

  const auto [frequencies, intensities] = DynamicsAnalyzer::calculateVDOS(vacf, time_step);

  ASSERT_FALSE(frequencies.empty());
  EXPECT_GE(frequencies.front(), static_cast<real_t>(0.0));

  const auto max_it = std::ranges::max_element(intensities);
  const size_t max_idx = std::distance(intensities.begin(), max_it);
  const real_t d_nu = frequencies[1] - frequencies[0];

  EXPECT_NEAR(frequencies[max_idx], target_nu_thz, 2.0 * d_nu);
  EXPECT_GT(*max_it, 0.0);

  // Near-zero frequencies for a 25 THz harmonic solid must stay negligible (no Filon jump)
  EXPECT_LT(intensities[0], (*max_it) * static_cast<real_t>(1e-3));
  EXPECT_LT(intensities[1], (*max_it) * static_cast<real_t>(1e-3));
  EXPECT_LT(intensities[2], (*max_it) * static_cast<real_t>(1e-3));
}

TEST(VDOSCalculatorTests, DiffusiveLiquidPeaksSmoothlyAtZeroFrequency) {
  // Exponentially decaying VACF e^(-at) yields a Lorentzian VDOS peaking smoothly at ν = 0.
  const real_t time_step = 1.0;
  const size_t num_frames = 2000;
  std::vector<real_t> vacf(num_frames);
  const real_t decay_rate = static_cast<real_t>(2.0 * correlation::math::PI * 10.0);

  for (size_t i = 0; i < num_frames; ++i) {
    const auto t_ps = static_cast<real_t>(static_cast<real_t>(i) * time_step * 0.001);
    vacf[i] = std::exp(-decay_rate * t_ps);
  }

  const auto [frequencies, intensities] = DynamicsAnalyzer::calculateVDOS(vacf, time_step);
  ASSERT_GT(frequencies.size(), 10U);
  const auto max_it = std::ranges::max_element(intensities);
  const size_t max_idx = std::distance(intensities.begin(), max_it);

  EXPECT_EQ(max_idx, 0U) << "Pure exponential VACF must peak smoothly at 0 THz.";
  EXPECT_GE(intensities[0], intensities[1]);
  EXPECT_GE(intensities[1], intensities[2]);
}

TEST(VDOSCalculatorTests, AllWindowTypesAndAreaNormalization) {
  using correlation::calculators::VDOSCalculator;
  using correlation::calculators::VDOSParams;

  const real_t time_step = 1.0;
  const size_t num_frames = 512;
  Histogram vacf_hist;
  vacf_hist.bins.resize(num_frames);
  std::vector<real_t> vacf_data(num_frames);
  for (size_t i = 0; i < num_frames; ++i) {
    vacf_hist.bins[i] = static_cast<real_t>(i) * time_step;
    vacf_data[i] = static_cast<real_t>(std::cos(2.0 * correlation::math::PI * 15.0 * i * 0.001));
  }
  vacf_hist.partials["Total"] = vacf_data;

  for (const auto win : {VDOSWindowType::Hann, VDOSWindowType::Blackman, VDOSWindowType::Gaussian,
                         VDOSWindowType::None}) {
    VDOSParams params;
    params.max_real_freq = static_cast<real_t>(500.0); // Full Nyquist to verify sum rule
    params.window_type = win;
    params.zero_pad_factor = 4;
    params.normalize = true;

    const auto hist = VDOSCalculator::calculate(vacf_hist, params);
    ASSERT_FALSE(hist.bins.empty());
    EXPECT_DOUBLE_EQ(hist.bins.front(), 0.0); // Strictly real angular frequencies (ω >= 0)
    EXPECT_EQ(hist.x_label, "ω");
    EXPECT_EQ(hist.x_unit, "rad/ps");
    EXPECT_EQ(hist.y_label, "g_vib(ω)");
    EXPECT_EQ(hist.y_unit, "ps/rad");
    EXPECT_DOUBLE_EQ(hist.partials.at("Total").front(), 0.0); // Vibrational-only vanishes at ω = 0

    const auto &dos = hist.partials.at("Total");
    double area = 0.0;
    for (size_t i = 0; i + 1 < hist.bins.size(); ++i) {
      EXPECT_GE(dos[i], static_cast<real_t>(0.0));
      area += 0.5 * (static_cast<double>(dos[i]) + static_cast<double>(dos[i + 1])) *
              (static_cast<double>(hist.bins[i + 1]) - static_cast<double>(hist.bins[i]));
    }
    EXPECT_NEAR(area, 1.0, 1e-3);
  }
}

TEST(VDOSCalculatorTests, MassWeightedPartialVDOSSumsToTotal) {
  using correlation::calculators::VDOSCalculator;
  using correlation::calculators::VDOSParams;
  using correlation::core::Cell;

  Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("O", {static_cast<real_t>(1.6), 0.0, 0.0});

  const size_t num_frames = 256;
  Histogram vacf_hist;
  vacf_hist.bins.resize(num_frames);
  std::vector<real_t> vacf_si(num_frames);
  std::vector<real_t> vacf_o(num_frames);
  std::vector<real_t> vacf_tot(num_frames);
  for (size_t i = 0; i < num_frames; ++i) {
    vacf_hist.bins[i] = static_cast<real_t>(i);
    vacf_si[i] = static_cast<real_t>(std::cos(2.0 * correlation::math::PI * 8.0 * i * 0.001));
    vacf_o[i] = static_cast<real_t>(std::cos(2.0 * correlation::math::PI * 20.0 * i * 0.001));
    vacf_tot[i] = static_cast<real_t>(0.5) * (vacf_si[i] + vacf_o[i]);
  }
  vacf_hist.partials["Total"] = vacf_tot;
  vacf_hist.partials["Si"] = vacf_si;
  vacf_hist.partials["O"] = vacf_o;

  VDOSParams params;
  params.max_real_freq = static_cast<real_t>(50.0);
  const auto hist = VDOSCalculator::calculate(vacf_hist, params, &cell);

  ASSERT_TRUE(hist.partials.contains("Si"));
  ASSERT_TRUE(hist.partials.contains("O"));
  const auto &tot = hist.partials.at("Total");
  const auto &p_si = hist.partials.at("Si");
  const auto &p_o = hist.partials.at("O");
  for (size_t i = 0; i < tot.size(); ++i) {
    EXPECT_NEAR(p_si[i] + p_o[i], tot[i], 1e-5);
  }
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

  const real_t max_computed_omega = vdos_hist.bins.back();
  const real_t max_computed_thz = vdos_hist.partials.at("Frequency_THz").back();
  EXPECT_GT(max_computed_thz, 45.0)
      << "Hydrogen VDOS frequency should reach ~50 THz instead of being truncated to 10 THz";
  EXPECT_LE(max_computed_thz, 50.0);
  EXPECT_NEAR(max_computed_omega, max_computed_thz * correlation::math::TWO_PI, 1e-2);
  EXPECT_TRUE(dists.hasHistogram("VDOS_Translational"));
}

TEST(VDOSCalculatorTests, TwoPhaseDecompositionSeparatesLiquidVibrationalAndTranslationalModes) {
  using correlation::calculators::VDOSCalculator;
  using correlation::calculators::VDOSParams;
  using correlation::core::Cell;

  // Liquid-like VACF: diffusive exponential decay + damped cage oscillation (like liquid Bi)
  Cell cell({15.0, 15.0, 15.0, 90.0, 90.0, 90.0});
  for (int i = 0; i < 64; ++i) {
    cell.addAtom("Bi", {static_cast<real_t>(i % 4), static_cast<real_t>((i / 4) % 4),
                        static_cast<real_t>(i / 16)});
  }

  const size_t num_frames = 512;
  Histogram vacf_hist;
  vacf_hist.bins.resize(num_frames);
  std::vector<real_t> vacf_data(num_frames);
  for (size_t i = 0; i < num_frames; ++i) {
    const double t_fs = static_cast<double>(i);
    vacf_hist.bins[i] = static_cast<real_t>(t_fs);
    // Diffusive tail + 3 THz damped vibrational mode
    const double val =
        0.3 * std::exp(-t_fs / 80.0) +
        0.7 * std::exp(-t_fs / 120.0) * std::cos(2.0 * correlation::math::PI * 3.0 * t_fs * 0.001);
    vacf_data[i] = static_cast<real_t>(val);
  }
  vacf_hist.partials["Total"] = vacf_data;
  vacf_hist.partials["Bi"] = vacf_data;

  VDOSParams params;
  params.max_real_freq = static_cast<real_t>(500.0);
  const auto decomposed = VDOSCalculator::calculateDecomposed(vacf_hist, params, &cell);

  EXPECT_GT(decomposed.fluidicity, static_cast<real_t>(0.0));
  EXPECT_LE(decomposed.fluidicity, static_cast<real_t>(1.0));

  const auto &g_vib = decomposed.vibrational.partials.at("Total");
  const auto &g_trans = decomposed.translational.partials.at("Total");
  ASSERT_EQ(g_vib.size(), g_trans.size());

  // Vibrational modes strictly tend to 0 at ω = 0; translational modes peak at ω = 0
  EXPECT_DOUBLE_EQ(g_vib[0], 0.0);
  EXPECT_GT(g_trans[0], static_cast<real_t>(0.0));
  EXPECT_GE(g_trans[0], g_trans[1]);

  // Both spectra are independently area-normalized to 1 over [0, ω_Nyquist]
  double vib_area = 0.0;
  double trans_area = 0.0;
  for (size_t i = 0; i + 1 < decomposed.vibrational.bins.size(); ++i) {
    const double d_omega = static_cast<double>(decomposed.vibrational.bins[i + 1]) -
                           static_cast<double>(decomposed.vibrational.bins[i]);
    vib_area += 0.5 * (static_cast<double>(g_vib[i]) + static_cast<double>(g_vib[i + 1])) * d_omega;
    trans_area +=
        0.5 * (static_cast<double>(g_trans[i]) + static_cast<double>(g_trans[i + 1])) * d_omega;
  }
  EXPECT_NEAR(vib_area, 1.0, 1e-3);
  EXPECT_NEAR(trans_area, 1.0, 1e-3);
}

} // namespace
} // namespace correlation::analysis
