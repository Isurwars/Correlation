// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <filesystem>
#include <gtest/gtest.h>
#include <vector>

#include "../include/DynamicsAnalyzer.hpp"
#include "../include/FileIO.hpp"
#include "../include/Trajectory.hpp"

// Assuming the test is run from the build directory or project root
// We need to locate the l-Bi.arc file
const std::string EXAMPLE_FILE = "examples/l-Bi/l-Bi.arc";

TEST(Test11_DynamicsAnalyzer, CalculatesVACFFromExampletraj) {
  // 1. Locate the file
  std::filesystem::path file_path = std::filesystem::absolute(EXAMPLE_FILE);

  if (!std::filesystem::exists(file_path)) {
    // Try going up levels if not found
    file_path = std::filesystem::path("../" + EXAMPLE_FILE);
    if (!std::filesystem::exists(file_path)) {
      file_path = std::filesystem::path("../../" + EXAMPLE_FILE);
    }
  }

  ASSERT_TRUE(std::filesystem::exists(file_path))
      << "Could not find example file: " << EXAMPLE_FILE;

  // 2. Read Trajectory
  Trajectory traj =
      FileIO::readTrajectory(file_path.string(), FileIO::FileType::Arc);
  ASSERT_GT(traj.getFrameCount(), 0) << "Trajectory should not be empty";

  // 3. Calculate Velocities
  traj.calculateVelocities();

  const auto &velocities = traj.getVelocities();
  ASSERT_EQ(velocities.size(), traj.getFrameCount());
  ASSERT_EQ(velocities[0].size(), traj.getFrames()[0].atomCount());

  // Check if we have some non-zero velocities (it's liquid Bi, particles move)
  double max_v_sq = 0.0;
  for (const auto &v : velocities[10]) { // Check some intermediate frame
    max_v_sq = std::max(max_v_sq, linalg::dot(v, v));
  }
  EXPECT_GT(max_v_sq, 0.0) << "Particles should be moving";

  // 4. Calculate VACF
  int max_lag = 50; // Calculate for 50 frames lag
  std::vector<double> vacf = DynamicsAnalyzer::calculateVACF(traj, max_lag);

  ASSERT_EQ(vacf.size(), max_lag + 1);

  // C(0) should be positive (autocorrelation at t=0 is <v^2>)
  EXPECT_GT(vacf[0], 0.0);

  // 5. Calculate Normalized VACF
  std::vector<double> norm_vacf =
      DynamicsAnalyzer::calculateNormalizedVACF(traj, max_lag);

  ASSERT_EQ(norm_vacf.size(), max_lag + 1);
  EXPECT_NEAR(norm_vacf[0], 1.0, 1e-5) << "Normalized VACF should start at 1.0";
}

TEST(Test11_DynamicsAnalyzer, CalculatesVDOSCorrectly) {
  // 1. Create synthetic VACF data: a simple cosine wave
  // v(t) = cos(2 * pi * f0 * t)
  // VDOS should show a peak at f0

  double dt = 1.0;          // 1 fs
  double f0 = 10.0;         // 10 THz frequency
  size_t num_frames = 1000; // 1 ps total time

  std::vector<double> vacf(num_frames);
  const double PI = 3.14159265358979323846;

  for (size_t i = 0; i < num_frames; ++i) {
    double t = i * dt;
    // f0 is in THz (10^12 Hz), t in fs (10^-15 s). product is 10^-3
    // cos(2 * pi * f0 * 10^12 * t * 10^-15) = cos(2 * pi * f0 * t * 0.001)
    vacf[i] = std::cos(2.0 * PI * f0 * t * 0.001);
  }

  // 2. Calculate VDOS
  auto [frequencies, intensities_real, intensities_imag] =
      DynamicsAnalyzer::calculateVDOS(vacf, dt);

  ASSERT_FALSE(frequencies.empty());
  ASSERT_EQ(frequencies.size(), intensities_real.size());
  ASSERT_EQ(frequencies.size(), intensities_imag.size());

  // 3. Find peak in real part
  auto max_it =
      std::max_element(intensities_real.begin(), intensities_real.end());
  size_t peak_idx = std::distance(intensities_real.begin(), max_it);
  double peak_freq = frequencies[peak_idx];

  // 4. Verify peak location
  EXPECT_NEAR(peak_freq, f0, 0.5)
      << "VDOS Peak should be near the source frequency";
}
