// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>
#include <vector>

#include "../include/DynamicsAnalyzer.hpp"

TEST(Test14_vDoS, VDOSIsNonZeroAtZeroFrequencyForConstantVACF) {
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
