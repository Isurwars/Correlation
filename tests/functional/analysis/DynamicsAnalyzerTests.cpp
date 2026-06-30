// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DynamicsAnalyzer.hpp"
#include "core/Trajectory.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"
#include "readers/FileReader.hpp"

#include <algorithm>
#include <filesystem>
#include <gtest/gtest.h>
#include <string_view>
#include <vector>

namespace correlation::analysis {

// Assuming the test is run from the build directory or project root
// We need to locate the l-Bi.arc file
constexpr std::string_view EXAMPLE_FILE = "examples/l-Bi/l-Bi.arc";

TEST(DynamicsAnalyzerTests, CalculatesVACFFromExampletraj) {
  // 1. Locate the file
  std::filesystem::path file_path = std::filesystem::absolute(EXAMPLE_FILE);

  if (!std::filesystem::exists(file_path)) {
    // Try going up levels if not found
    file_path = std::filesystem::path("../") / EXAMPLE_FILE;
    if (!std::filesystem::exists(file_path)) {
      file_path = std::filesystem::path("../../") / EXAMPLE_FILE;
    }
  }

  ASSERT_TRUE(std::filesystem::exists(file_path)) << "Could not find example file: " << EXAMPLE_FILE;

  // 2. Read correlation::core::Trajectory
  correlation::core::Trajectory traj =
      correlation::readers::readTrajectory(file_path.string(), correlation::readers::FileType::Arc);
  ASSERT_GT(traj.getFrameCount(), 0) << "correlation::core::Trajectory should not be empty";

  // 3. Calculate Velocities
  traj.calculateVelocities();

  const auto &frames = traj.getFrames();
  ASSERT_FALSE(frames.empty());
  EXPECT_NE(frames[1].atoms()[0].velocity().x(), 0.0);
  ASSERT_EQ(frames[0].atomCount(), traj.getFrames()[0].atomCount());

  // Check if we have some non-zero velocities (it's liquid Bi, particles move)
  double max_v_sq = 0.0;
  for (const auto &atom : frames[10].atoms()) { // Check some intermediate frame
    auto velocity = atom.velocity();
    max_v_sq = std::max(max_v_sq, correlation::math::dot(velocity, velocity));
  }
  EXPECT_GT(max_v_sq, 0.0) << "Particles should be moving";

  // 4. Calculate VACF
  int const max_lag = 50; // Calculate for 50 frames lag
  std::vector<double> vacf = DynamicsAnalyzer::calculateVACF(traj, correlation::analysis::MaxFrames{max_lag});

  ASSERT_EQ(vacf.size(), max_lag + 1);

  // C(0) should be positive (autocorrelation at t=0 is <v^2>)
  EXPECT_GT(vacf[0], 0.0);

  // 5. Calculate Normalized VACF
  std::vector<double> norm_vacf =
      DynamicsAnalyzer::calculateNormalizedVACF(traj, correlation::analysis::MaxFrames{max_lag});

  ASSERT_EQ(norm_vacf.size(), max_lag + 1);
  EXPECT_NEAR(norm_vacf[0], 1.0, 1e-5) << "Normalized VACF should start at 1.0";
}

TEST(DynamicsAnalyzerTests, CalculatesVDOSCorrectly) {
  // 1. Create synthetic VACF data: a simple cosine wave
  // v(t) = cos(2 * pi * f0 * t)
  // VDOS should show a peak at f0

  double const time_step = 1.0;   // 1 fs
  double const frequency = 10.0;  // 10 THz frequency
  size_t const num_frames = 1000; // 1 ps total time

  std::vector<double> vacf(num_frames);
  for (size_t i = 0; i < num_frames; ++i) {
    double const time = static_cast<double>(i) * time_step;
    vacf[i] = std::cos(2.0 * correlation::math::pi * frequency * time * 0.001);
  }

  // 2. Calculate VDOS
  auto [frequencies, intensities_real, intensities_imag] = DynamicsAnalyzer::calculateVDOS(vacf, time_step);

  ASSERT_FALSE(frequencies.empty());
  ASSERT_EQ(frequencies.size(), intensities_real.size());
  ASSERT_EQ(frequencies.size(), intensities_imag.size());

  // 3. Find peak in real part
  auto max_it = std::ranges::max_element(intensities_real);
  size_t const peak_idx = std::distance(intensities_real.begin(), max_it);
  double const peak_freq = frequencies[peak_idx];

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
  std::vector<double> time = {0.0, 1.0, 2.0, 3.0, 4.0};
  // Fit is done on the second half of the data.
  // time.size() / 2 = 2.
  // Second half indices are: 2, 3, 4.
  // time values: 2.0, 3.0, 4.0.
  // Let's choose msd = 6 * D * time.
  // If D = 0.5, slope should be 6 * 0.5 = 3.0.
  // So msd = 3.0 * time.
  // msd values: 6.0, 9.0, 12.0.
  std::vector<double> msd = {0.0, 3.0, 6.0, 9.0, 12.0};

  double d_coef = DynamicsAnalyzer::computeDiffusionCoefficientMSD(time, msd);
  EXPECT_NEAR(d_coef, 0.5, 1e-6);

  // Negative slope test should return 0.0
  std::vector<double> msd_neg = {0.0, -3.0, -6.0, -9.0, -12.0};
  double d_coef_neg = DynamicsAnalyzer::computeDiffusionCoefficientMSD(time, msd_neg);
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
  std::vector<double> time = {0.0, 1.0, 2.0};
  std::vector<double> vacf = {1.0, 1.0, 1.0};

  double d_coef = DynamicsAnalyzer::computeDiffusionCoefficientVACF(time, vacf);
  EXPECT_NEAR(d_coef, 2.0 / 3.0, 1e-6);

  // Time step <= 0.0 should return 0.0
  std::vector<double> invalid_time = {0.0, 0.0, 1.0};
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
  std::vector<double> time = {0.0, 1.0, 2.0};
  std::vector<double> norm_vacf = {1.0, 0.5, 0.0};

  double tau = DynamicsAnalyzer::computeRelaxationTime(time, norm_vacf);
  EXPECT_NEAR(tau, 1.0, 1e-6);

  // Mismatched sizes or empty vectors should return 0.0
  EXPECT_DOUBLE_EQ(DynamicsAnalyzer::computeRelaxationTime({}, {}), 0.0);
}

TEST(DynamicsAnalyzerTests, HandlesEmptyAndInvalidTrajectories) {
  correlation::core::Trajectory empty_traj;
  EXPECT_TRUE(DynamicsAnalyzer::calculateVACF(empty_traj, MaxFrames{5}).empty());
  EXPECT_TRUE(DynamicsAnalyzer::calculateMSD(empty_traj, MaxFrames{5}).empty());

  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {1.0, 1.0, 1.0});
  correlation::core::Trajectory trajectory;
  trajectory.addFrame(cell);
  trajectory.addFrame(cell);

  // start_frame >= end_frame
  EXPECT_TRUE(DynamicsAnalyzer::calculateVACF(trajectory, MaxFrames{1}, StartFrame{1}, EndFrame{1}).empty());
  EXPECT_TRUE(DynamicsAnalyzer::calculateMSD(trajectory, MaxFrames{1}, StartFrame{1}, EndFrame{1}).empty());
}

} // namespace correlation::analysis
