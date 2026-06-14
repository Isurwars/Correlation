// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "math/Smoothing.hpp"

#include <gtest/gtest.h>
#include <numeric>
#include <vector>

namespace correlation::testing {

using namespace correlation::math;

class SmoothingTests : public ::testing::Test {};

TEST_F(SmoothingTests, GenerateKernelNormalizesAndCalculatesCorrectly) {
  size_t const size = 11;
  double const dx = 0.1;
  double const sigma = 0.3;

  // Gaussian
  auto k_gauss = generateKernel(size, dx, sigma, KernelType::Gaussian);
  EXPECT_EQ(k_gauss.size(), size);
  double const sum_gauss = std::accumulate(k_gauss.begin(), k_gauss.end(), 0.0);
  EXPECT_NEAR(sum_gauss, 1.0, 1e-9);

  // Triweight
  auto k_tri = generateKernel(size, dx, sigma, KernelType::Triweight);
  EXPECT_EQ(k_tri.size(), size);
  double const sum_tri = std::accumulate(k_tri.begin(), k_tri.end(), 0.0);
  EXPECT_NEAR(sum_tri, 1.0, 1e-9);

  // Bump
  auto k_bump = generateKernel(size, dx, sigma, KernelType::Bump);
  EXPECT_EQ(k_bump.size(), size);
  double const sum_bump = std::accumulate(k_bump.begin(), k_bump.end(), 0.0);
  EXPECT_NEAR(sum_bump, 1.0, 1e-9);

  // Invalid kernel type throws invalid_argument
  EXPECT_THROW(generateKernel(size, dx, sigma, static_cast<KernelType>(-1)), std::invalid_argument);
}

TEST_F(SmoothingTests, KernelSmoothingBoundaryInputsAndErrors) {
  // Empty data
  std::vector<double> const empty_y;
  auto empty_res = KernelSmoothing(0.1, empty_y, 0.5, KernelType::Gaussian);
  EXPECT_TRUE(empty_res.empty());

  // Invalid dx
  std::vector<double> const y = {1.0, 2.0, 3.0, 4.0, 5.0};
  EXPECT_THROW(KernelSmoothing(0.0, y, 0.5, KernelType::Gaussian), std::invalid_argument);
  EXPECT_THROW(KernelSmoothing(-0.1, y, 0.5, KernelType::Gaussian), std::invalid_argument);

  // Invalid sigma
  EXPECT_THROW(KernelSmoothing(0.1, y, 0.0, KernelType::Gaussian), std::invalid_argument);
  EXPECT_THROW(KernelSmoothing(0.1, y, -0.5, KernelType::Gaussian), std::invalid_argument);
}

TEST_F(SmoothingTests, KernelSmoothingComputesConvolutions) {
  // Simple step signal
  std::vector<double> const y = {0.0, 0.0, 0.0, 10.0, 10.0, 10.0};
  double const dx = 0.2;
  double const sigma = 0.4;

  auto smoothed = KernelSmoothing(dx, y, sigma, KernelType::Gaussian);
  ASSERT_EQ(smoothed.size(), y.size());

  // Smoothing step function should blur the boundary.
  // The first element should be smoothed using left clamping (mostly 0.0, but slightly pulled up if kernel overlaps)
  // The transition should be intermediate.
  EXPECT_GT(smoothed[3], smoothed[0]);
  EXPECT_LT(smoothed[3], 10.0);
}

TEST_F(SmoothingTests, CompatibilityOverloadDerivesBinWidth) {
  std::vector<double> const r = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
  std::vector<double> const y = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  double const sigma = 0.5;

  auto smoothed_r = KernelSmoothing(r, y, sigma, KernelType::Gaussian);
  auto smoothed_dx = KernelSmoothing(0.2, y, sigma, KernelType::Gaussian);

  ASSERT_EQ(smoothed_r.size(), smoothed_dx.size());
  for (size_t i = 0; i < smoothed_r.size(); ++i) {
    EXPECT_DOUBLE_EQ(smoothed_r[i], smoothed_dx[i]);
  }

  // Size mismatch returns empty
  std::vector<double> const r_bad = {0.0, 0.2};
  EXPECT_TRUE(KernelSmoothing(r_bad, y, sigma, KernelType::Gaussian).empty());

  // Size < 2 returns empty
  std::vector<double> const r_short = {0.0};
  std::vector<double> const y_short = {1.0};
  EXPECT_TRUE(KernelSmoothing(r_short, y_short, sigma, KernelType::Gaussian).empty());
}

} // namespace correlation::testing
