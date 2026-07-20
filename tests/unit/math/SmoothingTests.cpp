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
namespace {
class SmoothingTests : public ::testing::Test {};
} // namespace

TEST_F(SmoothingTests, GenerateKernelNormalizesAndCalculatesCorrectly) {
  size_t const size = 11;
  real_t const bin_width = 0.1;
  real_t const sigma = 0.3;

  // Gaussian
  auto k_gauss = generateKernel({size, bin_width, sigma, KernelType::Gaussian});
  EXPECT_EQ(k_gauss.size(), size);
  real_t const sum_gauss = std::accumulate(k_gauss.begin(), k_gauss.end(), 0.0);
  EXPECT_NEAR(sum_gauss, 1.0, 1e-9);

  // Triweight
  auto k_tri = generateKernel({size, bin_width, sigma, KernelType::Triweight});
  EXPECT_EQ(k_tri.size(), size);
  real_t const sum_tri = std::accumulate(k_tri.begin(), k_tri.end(), 0.0);
  EXPECT_NEAR(sum_tri, 1.0, 1e-9);

  // Bump
  auto k_bump = generateKernel({size, bin_width, sigma, KernelType::Bump});
  EXPECT_EQ(k_bump.size(), size);
  real_t const sum_bump = std::accumulate(k_bump.begin(), k_bump.end(), 0.0);
  EXPECT_NEAR(sum_bump, 1.0, 1e-9);

  // Epanechnikov
  auto k_epan = generateKernel({size, bin_width, sigma, KernelType::Epanechnikov});
  EXPECT_EQ(k_epan.size(), size);
  real_t const sum_epan = std::accumulate(k_epan.begin(), k_epan.end(), 0.0);
  EXPECT_NEAR(sum_epan, 1.0, 1e-9);

  // Cosine
  auto k_cos = generateKernel({size, bin_width, sigma, KernelType::Cosine});
  EXPECT_EQ(k_cos.size(), size);
  real_t const sum_cos = std::accumulate(k_cos.begin(), k_cos.end(), 0.0);
  EXPECT_NEAR(sum_cos, 1.0, 1e-9);

  // Biweight
  auto k_biw = generateKernel({size, bin_width, sigma, KernelType::Biweight});
  EXPECT_EQ(k_biw.size(), size);
  real_t const sum_biw = std::accumulate(k_biw.begin(), k_biw.end(), 0.0);
  EXPECT_NEAR(sum_biw, 1.0, 1e-9);

  // Invalid kernel type throws invalid_argument
  EXPECT_THROW((void)generateKernel({size, bin_width, sigma, static_cast<KernelType>(-1)}), std::invalid_argument);
}

TEST_F(SmoothingTests, KernelSmoothingBoundaryInputsAndErrors) {
  // Empty data
  std::vector<real_t> const empty_y;
  auto empty_res = KernelSmoothing(0.1, empty_y, 0.5, KernelType::Gaussian);
  EXPECT_TRUE(empty_res.empty());

  // Invalid bin_width
  std::vector<real_t> const y_values = {1.0, 2.0, 3.0, 4.0, 5.0};
  EXPECT_THROW(KernelSmoothing(0.0, y_values, 0.5, KernelType::Gaussian), std::invalid_argument);
  EXPECT_THROW(KernelSmoothing(-0.1, y_values, 0.5, KernelType::Gaussian), std::invalid_argument);

  // Invalid sigma
  EXPECT_THROW(KernelSmoothing(0.1, y_values, 0.0, KernelType::Gaussian), std::invalid_argument);
  EXPECT_THROW(KernelSmoothing(0.1, y_values, -0.5, KernelType::Gaussian), std::invalid_argument);
}

TEST_F(SmoothingTests, KernelSmoothingComputesConvolutions) {
  // Simple step signal
  std::vector<real_t> const y_values = {0.0, 0.0, 0.0, 10.0, 10.0, 10.0};
  real_t const bin_width = 0.2;
  real_t const sigma = 0.4;

  auto smoothed = KernelSmoothing(bin_width, y_values, sigma, KernelType::Gaussian);
  ASSERT_EQ(smoothed.size(), y_values.size());

  // Smoothing step function should blur the boundary.
  // The first element should be smoothed using left clamping (mostly 0.0, but slightly pulled up if kernel overlaps)
  // The transition should be intermediate.
  EXPECT_GT(smoothed[3], smoothed[0]);
  EXPECT_LT(smoothed[3], 10.0);
}

TEST_F(SmoothingTests, CompatibilityOverloadDerivesBinWidth) {
  std::vector<real_t> const r_values = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
  std::vector<real_t> const y_values = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  real_t const sigma = 0.5;

  auto smoothed_r = KernelSmoothing(r_values, y_values, sigma, KernelType::Gaussian);
  auto smoothed_dx = KernelSmoothing(0.2, y_values, sigma, KernelType::Gaussian);

  ASSERT_EQ(smoothed_r.size(), smoothed_dx.size());
  for (size_t i = 0; i < smoothed_r.size(); ++i) {
    EXPECT_DOUBLE_EQ(smoothed_r[i], smoothed_dx[i]);
  }

  // Size mismatch returns empty
  std::vector<real_t> const r_bad = {0.0, 0.2};
  EXPECT_TRUE(KernelSmoothing(r_bad, y_values, sigma, KernelType::Gaussian).empty());

  // Size < 2 returns empty
  std::vector<real_t> const r_short = {0.0};
  std::vector<real_t> const y_short = {1.0};
  EXPECT_TRUE(KernelSmoothing(r_short, y_short, sigma, KernelType::Gaussian).empty());
}

} // namespace correlation::testing
