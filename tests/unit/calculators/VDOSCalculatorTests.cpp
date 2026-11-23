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
#include <cmath>
#include <gtest/gtest.h>
#include <iterator>
#include <vector>

using correlation::real_t;

namespace correlation::analysis {
namespace {

struct HarmonicVacfSpec {
  size_t num_frames{256};
  double freq_thz{10.0};
  double time_step{1.0};
  double scale{1.0};
};

std::vector<real_t> makeHarmonicVacf(const HarmonicVacfSpec &spec) {
  std::vector<real_t> vacf(spec.num_frames);
  for (size_t i = 0; i < spec.num_frames; ++i) {
    const double t_fs = static_cast<double>(i) * spec.time_step;
    vacf[i] = static_cast<real_t>(
        spec.scale * std::cos(2.0 * correlation::math::PI * spec.freq_thz * t_fs * 0.001));
  }
  return vacf;
}

Histogram makeHarmonicVacfHistogram(const HarmonicVacfSpec &spec) {
  Histogram vacf_hist;
  vacf_hist.bins.resize(spec.num_frames);
  for (size_t i = 0; i < spec.num_frames; ++i) {
    vacf_hist.bins[i] = static_cast<real_t>(static_cast<double>(i) * spec.time_step);
  }
  vacf_hist.partials["Total"] = makeHarmonicVacf(spec);
  return vacf_hist;
}

double computeTrapezoidalArea(const std::vector<real_t> &x_vals,
                              const std::vector<real_t> &y_vals) {
  double area = 0.0;
  for (size_t i = 0; i + 1 < x_vals.size(); ++i) {
    const double delta_x = static_cast<double>(x_vals[i + 1]) - static_cast<double>(x_vals[i]);
    area += 0.5 * (static_cast<double>(y_vals[i]) + static_cast<double>(y_vals[i + 1])) * delta_x;
  }
  return area;
}

void verifyNearZeroNegligible(const std::vector<real_t> &intensities, real_t peak_val) {
  const real_t threshold = peak_val * static_cast<real_t>(1e-3);
  EXPECT_LT(intensities[0], threshold);
  EXPECT_LT(intensities[1], threshold);
  EXPECT_LT(intensities[2], threshold);
}

void verifyVdosHistogramMetadata(const Histogram &hist) {
  ASSERT_FALSE(hist.bins.empty());
  EXPECT_DOUBLE_EQ(hist.bins.front(), 0.0); // Strictly real frequencies (THz >= 0)
  EXPECT_EQ(hist.x_label, "ω");
  EXPECT_EQ(hist.x_unit, "THz");
  EXPECT_EQ(hist.y_label, "g_vib(ω)");
  EXPECT_EQ(hist.y_unit, "THz⁻¹");
}

void verifyWindowNormalization(const Histogram &vacf_hist, VDOSWindowType win) {
  using correlation::calculators::VDOSCalculator;
  using correlation::calculators::VDOSParams;

  VDOSParams params;
  params.max_real_freq = static_cast<real_t>(500.0); // Full Nyquist to verify sum rule
  params.window_type = win;
  params.zero_pad_factor = 4;
  params.normalize = true;

  const auto hist = VDOSCalculator::calculate(vacf_hist, params);
  verifyVdosHistogramMetadata(hist);
  const auto &dos = hist.partials.at("Total");
  EXPECT_DOUBLE_EQ(dos.front(), 0.0); // Vibrational-only vanishes at ω = 0
  EXPECT_GE(*std::ranges::min_element(dos), static_cast<real_t>(0.0));
  EXPECT_NEAR(computeTrapezoidalArea(hist.bins, dos), 1.0, 1e-3);
}

void verifyTwoPhaseSpectra(const correlation::calculators::VDOSDecomposition &decomposed) {
  const auto &g_vib = decomposed.vibrational.partials.at("Total");
  const auto &g_trans = decomposed.translational.partials.at("Total");
  ASSERT_EQ(g_vib.size(), g_trans.size());

  // Vibrational modes strictly tend to 0 at ω = 0; translational modes peak at ω = 0
  EXPECT_DOUBLE_EQ(g_vib[0], 0.0);
  EXPECT_GT(g_trans[0], static_cast<real_t>(0.0));
  EXPECT_GE(g_trans[0], g_trans[1]);

  // Both spectra are independently area-normalized to 1 over [0, ν_Nyquist]
  EXPECT_NEAR(computeTrapezoidalArea(decomposed.vibrational.bins, g_vib), 1.0, 1e-3);
  EXPECT_NEAR(computeTrapezoidalArea(decomposed.translational.bins, g_trans), 1.0, 1e-3);
}

void verifyAuxiliaryBinUnits(const Histogram &norm_1x, const Histogram &norm_3x, size_t idx) {
  EXPECT_NEAR(norm_1x.partials.at("Total")[idx], norm_3x.partials.at("Total")[idx], 1e-4);
  const real_t thz = norm_1x.bins[idx];
  EXPECT_NEAR(norm_1x.partials.at("Frequency_cm_1")[idx],
              thz * static_cast<real_t>(correlation::math::THZ_TO_CMINV), 1e-2);
  EXPECT_NEAR(norm_1x.partials.at("Frequency_meV")[idx],
              thz * static_cast<real_t>(correlation::math::THZ_TO_MEV), 1e-3);
}

void verifyScaleInvarianceAndAuxiliaryUnits(const Histogram &norm_1x, const Histogram &norm_3x) {
  ASSERT_EQ(norm_1x.bins.size(), norm_3x.bins.size());
  EXPECT_FALSE(norm_1x.partials.contains("Frequency_THz"));
  for (size_t i = 0; i < norm_1x.bins.size(); ++i) {
    verifyAuxiliaryBinUnits(norm_1x, norm_3x, i);
  }
}

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
  const real_t target_nu_thz = 25.0;
  const auto vacf = makeHarmonicVacf(
      {.num_frames = num_frames, .freq_thz = target_nu_thz, .time_step = time_step});

  const auto [frequencies, intensities] = DynamicsAnalyzer::calculateVDOS(vacf, time_step);

  ASSERT_FALSE(frequencies.empty());
  EXPECT_GE(frequencies.front(), static_cast<real_t>(0.0));

  const auto max_it = std::ranges::max_element(intensities);
  const auto max_idx = static_cast<size_t>(std::distance(intensities.begin(), max_it));
  const real_t d_nu = frequencies[1] - frequencies[0];

  EXPECT_NEAR(frequencies[max_idx], target_nu_thz, 2.0 * d_nu);
  EXPECT_GT(*max_it, 0.0);
  verifyNearZeroNegligible(intensities, *max_it);
}

TEST(VDOSCalculatorTests, DiffusiveLiquidPeaksSmoothlyAtZeroFrequency) {
  // Exponentially decaying VACF e^(-at) yields a Lorentzian VDOS peaking smoothly at ν = 0.
  const real_t time_step = 1.0;
  const size_t num_frames = 2000;
  std::vector<real_t> vacf(num_frames);
  const auto decay_rate = static_cast<real_t>(2.0 * correlation::math::PI * 10.0);

  for (size_t i = 0; i < num_frames; ++i) {
    const auto t_ps = static_cast<real_t>(static_cast<double>(i) * time_step * 0.001);
    vacf[i] = std::exp(-decay_rate * t_ps);
  }

  const auto [frequencies, intensities] = DynamicsAnalyzer::calculateVDOS(vacf, time_step);
  ASSERT_GT(frequencies.size(), 10U);
  const auto max_it = std::ranges::max_element(intensities);
  const auto max_idx = static_cast<size_t>(std::distance(intensities.begin(), max_it));

  EXPECT_EQ(max_idx, 0U) << "Pure exponential VACF must peak smoothly at 0 THz.";
  EXPECT_GE(intensities[0], intensities[1]);
  EXPECT_GE(intensities[1], intensities[2]);
}

TEST(VDOSCalculatorTests, AllWindowTypesAndAreaNormalization) {
  const Histogram vacf_hist =
      makeHarmonicVacfHistogram({.num_frames = 512, .freq_thz = 15.0, .time_step = 1.0});

  for (const auto win : {VDOSWindowType::Hann, VDOSWindowType::Blackman, VDOSWindowType::Gaussian,
                         VDOSWindowType::None}) {
    verifyWindowNormalization(vacf_hist, win);
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
  const auto vacf_si = makeHarmonicVacf({.num_frames = num_frames, .freq_thz = 8.0});
  const auto vacf_o = makeHarmonicVacf({.num_frames = num_frames, .freq_thz = 20.0});
  Histogram vacf_hist;
  vacf_hist.bins.resize(num_frames);
  std::vector<real_t> vacf_tot(num_frames);
  for (size_t i = 0; i < num_frames; ++i) {
    vacf_hist.bins[i] = static_cast<real_t>(i);
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
  const Cell empty_cell;
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
  const Trajectory traj({h_cell, h_cell}, 1.0);

  DistributionFunctions dists(h_cell);
  // Synthetic VACF with 50 frames
  Histogram vacf_hist;
  vacf_hist.bins.resize(50);
  for (size_t i = 0; i < 50; ++i) {
    vacf_hist.bins[i] = static_cast<real_t>(i); // dt = 1 fs
  }
  const std::vector<real_t> vacf_data(50, 1.0);
  vacf_hist.partials["Total"] = vacf_data;
  dists.addHistogram("VACF", std::move(vacf_hist));

  const VDOSCalculator calc;
  const AnalysisSettings settings;
  calc.calculateTrajectory(dists, traj, settings);

  ASSERT_TRUE(dists.hasHistogram("VDOS"));
  const auto &vdos_hist = dists.getHistogram("VDOS");
  ASSERT_FALSE(vdos_hist.bins.empty());

  const real_t max_computed_thz = vdos_hist.bins.back();
  EXPECT_GT(max_computed_thz, 45.0)
      << "Hydrogen VDOS frequency should reach ~50 THz instead of being truncated to 10 THz";
  EXPECT_LE(max_computed_thz, 50.0);
  EXPECT_FALSE(vdos_hist.partials.contains("Frequency_THz"));
  EXPECT_TRUE(dists.hasHistogram("VDOS_Translational"));
}

TEST(VDOSCalculatorTests, TwoPhaseDecompositionSeparatesLiquidVibrationalAndTranslationalModes) {
  using correlation::calculators::VDOSCalculator;
  using correlation::calculators::VDOSParams;
  using correlation::core::Cell;

  // Liquid-like VACF: diffusive exponential decay + damped cage oscillation (like liquid Bi)
  Cell cell({15.0, 15.0, 15.0, 90.0, 90.0, 90.0});
  for (int i = 0; i < 64; ++i) {
    const int pos_x = i % 4;
    const int pos_y = (i / 4) % 4;
    const int pos_z = i / 16;
    cell.addAtom(
        "Bi", {static_cast<real_t>(pos_x), static_cast<real_t>(pos_y), static_cast<real_t>(pos_z)});
  }

  const size_t num_frames = 512;
  Histogram vacf_hist;
  vacf_hist.bins.resize(num_frames);
  std::vector<real_t> vacf_data(num_frames);
  for (size_t i = 0; i < num_frames; ++i) {
    const auto t_fs = static_cast<double>(i);
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
  verifyTwoPhaseSpectra(decomposed);
}

TEST(VDOSCalculatorTests, ThrowsWhenVACFTooShortAndSkipsMissingVACFInTrajectory) {
  using correlation::calculators::VDOSCalculator;
  using correlation::core::Cell;
  using correlation::core::Trajectory;

  Histogram short_vacf;
  short_vacf.bins = {0.0};
  short_vacf.partials["Total"] = {1.0};
  EXPECT_THROW((void)VDOSCalculator::calculate(short_vacf), std::logic_error);
  EXPECT_THROW((void)VDOSCalculator::calculateDecomposed(short_vacf), std::logic_error);

  Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {0.0, 0.0, 0.0});
  const Trajectory traj({cell, cell}, 1.0);
  DistributionFunctions dists(cell);

  const VDOSCalculator calc;
  const AnalysisSettings settings;
  calc.calculateTrajectory(dists, traj, settings);
  EXPECT_FALSE(dists.hasHistogram("VDOS"));
  EXPECT_FALSE(dists.hasHistogram("VDOS_Translational"));

  dists.addHistogram("VACF", std::move(short_vacf));
  calc.calculateTrajectory(dists, traj, settings);
  EXPECT_FALSE(dists.hasHistogram("VDOS"));
}

TEST(VDOSCalculatorTests, PureHarmonicSolidHasZeroFluidicityAndVanishingZeroFrequency) {
  using correlation::calculators::VDOSCalculator;
  using correlation::calculators::VDOSParams;
  using correlation::core::Cell;

  Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  const real_t target_thz = 20.0;
  Histogram vacf_hist = makeHarmonicVacfHistogram({.num_frames = 1024, .freq_thz = target_thz});
  vacf_hist.partials["Si"] = vacf_hist.partials.at("Total");

  VDOSParams params;
  params.max_real_freq = static_cast<real_t>(50.0);
  params.window_type = VDOSWindowType::Blackman;
  const auto decomposed = VDOSCalculator::calculateDecomposed(vacf_hist, params, &cell);

  EXPECT_LT(decomposed.fluidicity, static_cast<real_t>(1e-2));
  const auto &g_vib = decomposed.vibrational.partials.at("Total");
  EXPECT_DOUBLE_EQ(g_vib[0], 0.0);

  const auto max_it = std::ranges::max_element(g_vib);
  const auto peak_idx = static_cast<size_t>(std::distance(g_vib.begin(), max_it));
  const real_t d_freq = decomposed.vibrational.bins[1] - decomposed.vibrational.bins[0];
  EXPECT_NEAR(decomposed.vibrational.bins[peak_idx], target_thz, 2.0 * d_freq);
}

TEST(VDOSCalculatorTests, PureExponentialWithoutCellTriggersFallbackFluidicity) {
  using correlation::calculators::VDOSCalculator;
  using correlation::calculators::VDOSParams;

  const size_t num_frames = 1024;
  Histogram vacf_hist;
  vacf_hist.bins.resize(num_frames);
  std::vector<real_t> vacf_data(num_frames);
  for (size_t i = 0; i < num_frames; ++i) {
    vacf_hist.bins[i] = static_cast<real_t>(i);
    vacf_data[i] = static_cast<real_t>(std::exp(-static_cast<double>(i) / 40.0));
  }
  vacf_hist.partials["Total"] = vacf_data;

  VDOSParams params;
  params.max_real_freq = static_cast<real_t>(100.0);
  const auto decomposed = VDOSCalculator::calculateDecomposed(vacf_hist, params, nullptr);

  EXPECT_GT(decomposed.fluidicity, static_cast<real_t>(0.5));
  EXPECT_LE(decomposed.fluidicity, static_cast<real_t>(1.0));
  EXPECT_DOUBLE_EQ(decomposed.vibrational.partials.at("Total")[0], 0.0);
  const auto &g_trans = decomposed.translational.partials.at("Total");
  EXPECT_GT(g_trans[0], static_cast<real_t>(0.0));
  EXPECT_GE(g_trans[0], g_trans[1]);
  EXPECT_GE(g_trans[1], g_trans[2]);
}

TEST(VDOSCalculatorTests, ZeroPadFactorIncreasesResolutionWithoutShiftingPeak) {
  using correlation::calculators::VDOSCalculator;
  using correlation::calculators::VDOSParams;

  const real_t target_thz = 12.5;
  const Histogram vacf_hist =
      makeHarmonicVacfHistogram({.num_frames = 256, .freq_thz = target_thz});

  VDOSParams coarse_params;
  coarse_params.max_real_freq = static_cast<real_t>(50.0);
  coarse_params.zero_pad_factor = 1;
  const auto coarse = VDOSCalculator::calculate(vacf_hist, coarse_params);

  VDOSParams fine_params = coarse_params;
  fine_params.zero_pad_factor = 16;
  const auto fine = VDOSCalculator::calculate(vacf_hist, fine_params);

  EXPECT_GT(fine.bins.size(), coarse.bins.size());
  const real_t d_freq_coarse = coarse.bins[1] - coarse.bins[0];
  const real_t d_freq_fine = fine.bins[1] - fine.bins[0];
  EXPECT_LT(d_freq_fine, d_freq_coarse);

  const auto &fine_dos = fine.partials.at("Total");
  const auto fine_peak =
      static_cast<size_t>(std::distance(fine_dos.begin(), std::ranges::max_element(fine_dos)));
  EXPECT_NEAR(fine.bins[fine_peak], target_thz, 2.0 * d_freq_fine);
}

TEST(VDOSCalculatorTests, NormalizationScaleInvarianceAndAuxiliaryUnits) {
  using correlation::calculators::VDOSCalculator;
  using correlation::calculators::VDOSParams;

  const Histogram vacf_1x =
      makeHarmonicVacfHistogram({.num_frames = 256, .freq_thz = 10.0, .scale = 1.0});
  const Histogram vacf_3x =
      makeHarmonicVacfHistogram({.num_frames = 256, .freq_thz = 10.0, .scale = 3.0});

  VDOSParams norm_params;
  norm_params.max_real_freq = static_cast<real_t>(30.0);
  norm_params.normalize = true;
  const auto norm_1x = VDOSCalculator::calculate(vacf_1x, norm_params);
  const auto norm_3x = VDOSCalculator::calculate(vacf_3x, norm_params);

  verifyScaleInvarianceAndAuxiliaryUnits(norm_1x, norm_3x);

  VDOSParams raw_params = norm_params;
  raw_params.normalize = false;
  const auto raw_hist = VDOSCalculator::calculate(vacf_1x, raw_params);
  EXPECT_EQ(raw_hist.y_unit, "arbitrary units");
}

TEST(VDOSCalculatorTests, MassWeightingReflectsAtomicMassRatio) {
  using correlation::calculators::VDOSCalculator;
  using correlation::calculators::VDOSParams;
  using correlation::core::Cell;

  Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("H", {static_cast<real_t>(1.5), 0.0, 0.0});

  const size_t num_frames = 512;
  const auto vacf_si = makeHarmonicVacf({.num_frames = num_frames, .freq_thz = 10.0});
  const auto vacf_h = makeHarmonicVacf({.num_frames = num_frames, .freq_thz = 35.0});
  Histogram vacf_hist;
  vacf_hist.bins.resize(num_frames);
  std::vector<real_t> vacf_tot(num_frames);
  for (size_t i = 0; i < num_frames; ++i) {
    vacf_hist.bins[i] = static_cast<real_t>(i);
    vacf_tot[i] = static_cast<real_t>(0.5) * (vacf_si[i] + vacf_h[i]);
  }
  vacf_hist.partials["Total"] = vacf_tot;
  vacf_hist.partials["Si"] = vacf_si;
  vacf_hist.partials["H"] = vacf_h;

  VDOSParams params;
  params.max_real_freq = static_cast<real_t>(500.0);
  const auto hist = VDOSCalculator::calculate(vacf_hist, params, &cell);

  const double area_si = computeTrapezoidalArea(hist.bins, hist.partials.at("Si"));
  const double area_h = computeTrapezoidalArea(hist.bins, hist.partials.at("H"));
  // m_Si / (m_Si + m_H) ≈ 28.085 / 29.093 ≈ 0.9654
  EXPECT_NEAR(area_si, 28.085 / (28.085 + 1.008), 0.02);
  EXPECT_NEAR(area_h, 1.008 / (28.085 + 1.008), 0.02);
}

TEST(VDOSCalculatorTests, CalculateTrajectoryPropagatesSettings) {
  using correlation::calculators::VDOSCalculator;
  using correlation::core::Cell;
  using correlation::core::Trajectory;

  Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  const Trajectory traj({cell, cell}, 1.0);

  const auto make_dists_with_vacf = [&] {
    DistributionFunctions dists(cell);
    Histogram vacf_hist = makeHarmonicVacfHistogram({.num_frames = 256, .freq_thz = 5.0});
    vacf_hist.partials["Si"] = vacf_hist.partials.at("Total");
    dists.addHistogram("VACF", std::move(vacf_hist));
    return dists;
  };

  const VDOSCalculator calc;
  AnalysisSettings settings_low;
  settings_low.vdos_window = VDOSWindowType::Hann;
  settings_low.vdos_zero_pad = 1;
  auto dists_low = make_dists_with_vacf();
  calc.calculateTrajectory(dists_low, traj, settings_low);

  AnalysisSettings settings_high;
  settings_high.vdos_window = VDOSWindowType::Blackman;
  settings_high.vdos_zero_pad = 16;
  auto dists_high = make_dists_with_vacf();
  calc.calculateTrajectory(dists_high, traj, settings_high);

  ASSERT_TRUE(dists_low.hasHistogram("VDOS"));
  ASSERT_TRUE(dists_high.hasHistogram("VDOS"));
  EXPECT_GT(dists_high.getHistogram("VDOS").bins.size(),
            dists_low.getHistogram("VDOS").bins.size());
}

} // namespace
} // namespace correlation::analysis
