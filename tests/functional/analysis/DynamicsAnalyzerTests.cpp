// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "analysis/DynamicsAnalyzer.hpp"
#include "calculators/temporal/VACFCalculator.hpp"
#include "calculators/temporal/VDOSCalculator.hpp"
#include "core/Trajectory.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"
#include "readers/FileReader.hpp"

#include <algorithm>
#include <filesystem>
#include <gtest/gtest.h>
#include <string_view>
#include <vector>

using correlation::real_t;

namespace correlation::analysis {

// Assuming the test is run from the build directory or project root
// We need to locate the l-Bi.arc file
constexpr std::string_view EXAMPLE_FILE = "examples/l-Bi/l-Bi.arc";

namespace {

[[nodiscard]] std::filesystem::path resolveExampleFilePath(std::string_view rel_path) {
  for (const auto &prefix : {"", "../", "../../"}) {
    const auto cand_path = std::filesystem::path(prefix) / rel_path;
    if (std::filesystem::exists(cand_path)) {
      return std::filesystem::absolute(cand_path);
    }
  }
  return std::filesystem::absolute(rel_path);
}

void checkNonZeroVelocities(const std::vector<correlation::core::Cell> &frames) {
  correlation::real_t max_v_sq = 0.0;
  for (size_t idx = 0; idx < frames[10].atomCount(); ++idx) {
    const auto velocity = frames[10].atomVelocity(idx);
    max_v_sq = std::max(max_v_sq, correlation::math::dot(velocity, velocity));
  }
  EXPECT_GT(max_v_sq, 0.0) << "Particles should be moving";
}

} // namespace

TEST(DynamicsAnalyzerTests, CalculatesVACFFromExampletraj) {
  const auto file_path = resolveExampleFilePath(EXAMPLE_FILE);
  ASSERT_TRUE(std::filesystem::exists(file_path))
      << "Could not find example file: " << EXAMPLE_FILE;

  correlation::core::Trajectory traj =
      correlation::readers::readTrajectory(file_path.string(), correlation::readers::FileType::Arc);
  ASSERT_GT(traj.getFrameCount(), 0) << "correlation::core::Trajectory should not be empty";

  traj.calculateVelocities();

  const auto &frames = traj.getFrames();
  ASSERT_FALSE(frames.empty());
  EXPECT_NE(frames[1].atomVelocity(0).x(), 0.0);
  ASSERT_EQ(frames[0].atomCount(), traj.getFrames()[0].atomCount());

  checkNonZeroVelocities(frames);

  int const max_lag = 50;
  const std::vector<real_t> vacf =
      DynamicsAnalyzer::calculateVACF(traj, correlation::analysis::MaxFrames{max_lag});

  ASSERT_EQ(vacf.size(), max_lag + 1);
  EXPECT_GT(vacf[0], 0.0);

  const std::vector<real_t> norm_vacf =
      DynamicsAnalyzer::calculateNormalizedVACF(traj, correlation::analysis::MaxFrames{max_lag});

  ASSERT_EQ(norm_vacf.size(), max_lag + 1);
  EXPECT_NEAR(norm_vacf[0], 1.0, 1e-5) << "Normalized VACF should start at 1.0";
}

TEST(DynamicsAnalyzerTests, CalculatesVDOSCorrectly) {
  // 1. Create synthetic VACF data: a simple cosine wave
  // v(t) = cos(2 * pi * f0 * t)
  // VDOS should show a peak at f0

  real_t const time_step = 1.0;   // 1 fs
  real_t const frequency = 10.0;  // 10 THz frequency
  size_t const num_frames = 1000; // 1 ps total time

  std::vector<real_t> vacf(num_frames);
  for (size_t i = 0; i < num_frames; ++i) {
    real_t const time = static_cast<real_t>(i) * time_step;
    vacf[i] = static_cast<real_t>(std::cos(2.0 * correlation::math::PI * frequency * time * 0.001));
  }

  // 2. Calculate VDOS
  auto [frequencies, intensities] = DynamicsAnalyzer::calculateVDOS(vacf, time_step);

  ASSERT_FALSE(frequencies.empty());
  ASSERT_EQ(frequencies.size(), intensities.size());

  // 3. Find peak in real spectrum
  const auto max_it = std::ranges::max_element(intensities);
  size_t const peak_idx = std::distance(intensities.begin(), max_it);
  real_t const peak_freq = frequencies[peak_idx];

  // 4. Verify peak location
  EXPECT_NEAR(peak_freq, frequency, 0.5) << "VDOS Peak should be near the source frequency";
}

TEST(DynamicsAnalyzerTests, CalculatesMSDCorrectly) {
  correlation::core::Trajectory trajectory;

  // Frame 0
  correlation::core::Cell cell0({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell0.addAtom("Ar", {2.0, 0.0, 0.0});
  trajectory.addFrame(cell0);

  // Frame 1
  correlation::core::Cell cell1({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell1.addAtom("Ar", {2.5, 0.0, 0.0});
  trajectory.addFrame(cell1);

  // Frame 2
  correlation::core::Cell cell2({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell2.addAtom("Ar", {3.0, 0.0, 0.0});
  trajectory.addFrame(cell2);

  // Frame 3
  correlation::core::Cell cell3({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell3.addAtom("Ar", {3.5, 0.0, 0.0});
  trajectory.addFrame(cell3);

  auto msd = DynamicsAnalyzer::calculateMSD(trajectory, MaxFrames{2});
  ASSERT_EQ(msd.size(), 3);
  EXPECT_NEAR(msd[0], 0.0, 1e-6);
  EXPECT_NEAR(msd[1], 0.25, 1e-6);
  EXPECT_NEAR(msd[2], 1.0, 1e-6);
}

TEST(DynamicsAnalyzerTests, ComputesDiffusionCoefficientMSD) {
  const std::vector<real_t> time = {0.0, 1.0, 2.0, 3.0, 4.0};
  // Fit is done on the second half of the data.
  // time.size() / 2 = 2.
  // Second half indices are: 2, 3, 4.
  // time values: 2.0, 3.0, 4.0.
  // Let's choose msd = 6 * D * time.
  // If D = 0.5, slope should be 6 * 0.5 = 3.0.
  // So msd = 3.0 * time.
  // msd values: 6.0, 9.0, 12.0.
  const std::vector<real_t> msd = {0.0, 3.0, 6.0, 9.0, 12.0};

  const real_t d_coef = DynamicsAnalyzer::computeDiffusionCoefficientMSD(time, msd);
  EXPECT_NEAR(d_coef, 0.5, 1e-6);

  // Negative slope test should return 0.0
  const std::vector<real_t> msd_neg = {0.0, -3.0, -6.0, -9.0, -12.0};
  const real_t d_coef_neg = DynamicsAnalyzer::computeDiffusionCoefficientMSD(time, msd_neg);
  EXPECT_DOUBLE_EQ(d_coef_neg, 0.0);

  // Under minimum required points (needs at least 2 points in second half)
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeDiffusionCoefficientMSD({0.0}, {0.0}), 0.0);
}

TEST(DynamicsAnalyzerTests, ComputesDiffusionCoefficientVACF) {
  // D = 1/3 * integral of VACF from t=0 to t_max.
  // Let's use time = {0.0, 1.0, 2.0};
  // vacf = {1.0, 1.0, 1.0};
  // Trapezoidal integral:
  // Step 1: 0.5 * (1.0 + 1.0) * 1.0 = 1.0
  // Step 2: 0.5 * (1.0 + 1.0) * 1.0 = 1.0
  // Total integral = 2.0.
  // D = 2.0 / 3.0 ≈ 0.666667.
  const std::vector<real_t> time = {0.0, 1.0, 2.0};
  const std::vector<real_t> vacf = {1.0, 1.0, 1.0};

  const real_t d_coef = DynamicsAnalyzer::computeDiffusionCoefficientVACF(time, vacf);
  EXPECT_NEAR(d_coef, 2.0 / 3.0, 1e-6);

  // Time step <= 0.0 should return 0.0
  const std::vector<real_t> invalid_time = {0.0, 0.0, 1.0};
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeDiffusionCoefficientVACF(invalid_time, vacf), 0.0);

  // Mismatched sizes or empty vectors should return 0.0
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeDiffusionCoefficientVACF({}, {}), 0.0);
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeDiffusionCoefficientVACF({0.0}, {1.0}), 0.0);
}

TEST(DynamicsAnalyzerTests, ComputesRelaxationTime) {
  // relaxation_time = integral of normalized_vacf.
  // time = {0.0, 1.0, 2.0};
  // norm_vacf = {1.0, 0.5, 0.0};
  // Trapezoidal integral:
  // Step 1: 0.5 * (1.0 + 0.5) * 1.0 = 0.75
  // Step 2: 0.5 * (0.5 + 0.0) * 1.0 = 0.25
  // Total = 1.0.
  const std::vector<real_t> time = {0.0, 1.0, 2.0};
  const std::vector<real_t> norm_vacf = {1.0, 0.5, 0.0};

  const real_t tau = DynamicsAnalyzer::computeRelaxationTime(time, norm_vacf);
  EXPECT_NEAR(tau, 1.0, 1e-6);

  // Mismatched sizes or empty vectors should return 0.0
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeRelaxationTime({}, {}), 0.0);
}

TEST(DynamicsAnalyzerTests, HandlesEmptyAndInvalidTrajectories) {
  const correlation::core::Trajectory empty_traj;
  EXPECT_TRUE(DynamicsAnalyzer::calculateVACF(empty_traj, MaxFrames{5}).empty());
  EXPECT_TRUE(DynamicsAnalyzer::calculateMSD(empty_traj, MaxFrames{5}).empty());

  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {1.0, 1.0, 1.0});
  correlation::core::Trajectory trajectory;
  trajectory.addFrame(cell);
  trajectory.addFrame(cell);

  // start_frame >= end_frame
  EXPECT_TRUE(DynamicsAnalyzer::calculateVACF(trajectory, MaxFrames{1}, StartFrame{1}, EndFrame{1})
                  .empty());
  EXPECT_TRUE(
      DynamicsAnalyzer::calculateMSD(trajectory, MaxFrames{1}, StartFrame{1}, EndFrame{1}).empty());

  // Trajectory with 0 atoms per frame
  const correlation::core::Cell zero_atom_cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  correlation::core::Trajectory zero_atom_traj;
  zero_atom_traj.addFrame(zero_atom_cell);
  zero_atom_traj.addFrame(zero_atom_cell);
  EXPECT_TRUE(DynamicsAnalyzer::calculateVACF(zero_atom_traj, MaxFrames{1}).empty());
  EXPECT_TRUE(DynamicsAnalyzer::calculateMSD(zero_atom_traj, MaxFrames{1}).empty());

  // Single-frame trajectory: invalid for temporal autocorrelation
  correlation::core::Trajectory single_frame_traj;
  single_frame_traj.addFrame(cell);
  EXPECT_TRUE(DynamicsAnalyzer::calculateVACF(single_frame_traj, MaxFrames{1}).empty());
  EXPECT_TRUE(DynamicsAnalyzer::calculateMSD(single_frame_traj, MaxFrames{1}).empty());
}

TEST(DynamicsAnalyzerTests, HandlesVDOSBoundaryAndEdgeCaseInputs) {
  // Empty VACF input
  auto [freqs_empty, vdos_empty] = DynamicsAnalyzer::calculateVDOS({}, 0.5);
  EXPECT_TRUE(freqs_empty.empty());
  EXPECT_TRUE(vdos_empty.empty());

  // Non-positive time step
  auto [freqs_neg, vdos_neg] =
      DynamicsAnalyzer::calculateVDOS({1.0, static_cast<real_t>(0.8), 0.5}, -1.0);
  EXPECT_TRUE(freqs_neg.empty());

  auto [freqs_zero, vdos_zero] =
      DynamicsAnalyzer::calculateVDOS({1.0, static_cast<real_t>(0.8), 0.5}, 0.0);
  EXPECT_TRUE(freqs_zero.empty());

  // Single-frame VACF: must not crash or underflow
  auto [freqs_single, vdos_single] = DynamicsAnalyzer::calculateVDOS({1.0}, 0.5);
  EXPECT_FALSE(freqs_single.empty());
  EXPECT_EQ(freqs_single.size(), vdos_single.size());
}

TEST(DynamicsAnalyzerTests, LiquidBismuthVDOSDecomposesVibrationalAndTranslationalModes) {
  const auto file_path = resolveExampleFilePath(EXAMPLE_FILE);
  ASSERT_TRUE(std::filesystem::exists(file_path))
      << "Could not find example file: " << EXAMPLE_FILE;

  correlation::core::Trajectory traj =
      correlation::readers::readTrajectory(file_path.string(), correlation::readers::FileType::Arc);
  ASSERT_GT(traj.getFrameCount(), 10U);
  traj.calculateVelocities();

  DistributionFunctions dists(traj.firstFrame());
  AnalysisSettings settings;
  settings.vdos_window = VDOSWindowType::Hann;
  settings.vdos_zero_pad = 4;

  const correlation::calculators::VACFCalculator vacf_calc;
  vacf_calc.calculateTrajectory(dists, traj, settings);
  ASSERT_TRUE(dists.hasHistogram("VACF"));

  const correlation::calculators::VDOSCalculator vdos_calc;
  vdos_calc.calculateTrajectory(dists, traj, settings);
  ASSERT_TRUE(dists.hasHistogram("VDOS"));
  ASSERT_TRUE(dists.hasHistogram("VDOS_Translational"));

  const auto &vib = dists.getHistogram("VDOS");
  const auto &trans = dists.getHistogram("VDOS_Translational");
  EXPECT_EQ(vib.x_label, "ω");
  EXPECT_EQ(vib.x_unit, "meV");
  EXPECT_EQ(trans.x_unit, "meV");
  EXPECT_EQ(vib.y_label, "g_vib(ω)");
  EXPECT_EQ(trans.y_label, "g_trans(ω)");

  const auto &g_vib = vib.partials.at("Total");
  const auto &g_trans = trans.partials.at("Total");
  ASSERT_FALSE(g_vib.empty());
  ASSERT_EQ(g_vib.size(), g_trans.size());

  // Liquid Bi: vibrational modes strictly vanish at ω = 0; translational modes peak at ω = 0
  EXPECT_DOUBLE_EQ(g_vib[0], 0.0);
  EXPECT_GT(*std::ranges::max_element(g_vib), static_cast<real_t>(0.0));
  EXPECT_GT(g_trans[0], static_cast<real_t>(0.0));
  EXPECT_GE(g_trans[0], g_trans[1]);
}

} // namespace correlation::analysis
