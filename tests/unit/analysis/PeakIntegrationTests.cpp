/**
 * @file PeakIntegrationTests.cpp
 * @brief Unit tests for PeakIntegrationService.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/PeakIntegrationService.hpp"
#include "math/Constants.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <vector>

using namespace correlation::analysis;

class PeakIntegrationTests : public ::testing::Test {
protected:
  static auto generateGaussianProfile(real_t mu = 2.5, real_t sigma = 0.15, real_t amp = 4.0,
                                      real_t x_min = 0.0, real_t x_max = 5.0, real_t dx = 0.01)
      -> std::pair<std::vector<real_t>, std::vector<real_t>> {
    std::vector<real_t> x_bins;
    std::vector<real_t> y_vals;
    const auto n_pts = static_cast<std::size_t>((x_max - x_min) / dx);
    x_bins.reserve(n_pts);
    y_vals.reserve(n_pts);

    for (std::size_t i = 0; i < n_pts; ++i) {
      const real_t x = x_min + static_cast<real_t>(i) * dx;
      x_bins.push_back(x);
      const real_t diff = (x - mu) / sigma;
      const real_t y = amp * std::exp(static_cast<real_t>(-0.5) * diff * diff);
      y_vals.push_back(y);
    }
    return {x_bins, y_vals};
  }
};

TEST_F(PeakIntegrationTests, IntegratesGaussianPeakAccurately) {
  constexpr real_t mu = 2.5;
  constexpr real_t sigma = 0.12;
  constexpr real_t amp = 5.0;
  const auto [x_bins, y_vals] = generateGaussianProfile(mu, sigma, amp);

  const auto result =
      PeakIntegrationService::integrateRange(x_bins, y_vals, mu - 4.0 * sigma, mu + 4.0 * sigma);

  const real_t analytical_area =
      amp * sigma * std::sqrt(static_cast<real_t>(2.0 * correlation::math::PI));
  EXPECT_NEAR(result.area, analytical_area, 0.02);
  EXPECT_NEAR(result.centroid, mu, 0.01);
  EXPECT_NEAR(result.r_peak, mu, 0.01);
  EXPECT_NEAR(result.peak_height, amp, 0.01);

  // FWHM of Gaussian: 2 * sqrt(2 * ln(2)) * sigma ≈ 2.35482 * sigma
  const real_t expected_fwhm = static_cast<real_t>(2.35482) * sigma;
  EXPECT_NEAR(result.fwhm, expected_fwhm, 0.01);
}

TEST_F(PeakIntegrationTests, ComputesCoordinationNumberForRadialShell) {
  // Model a single coordination shell of 8 neighbors at r = 2.48 Å in a BCC cell
  constexpr real_t r0 = 2.48;
  constexpr real_t sigma = 0.08;
  constexpr real_t rho = 0.085; // macroscopic number density (atoms / Å³)

  // In radial shell: CN = 4πρ ∫ r² g(r) dr.
  // If g(r) = A * exp(-(r-r0)² / 2σ²), then ∫ r² g(r) dr ≈ r0² A σ sqrt(2π).
  // For CN = 8.0: A = 8.0 / (4πρ r0² σ sqrt(2π))
  constexpr real_t four_pi = correlation::math::FOUR_PI;
  constexpr real_t sqrt_2pi = static_cast<real_t>(std::sqrt(2.0 * correlation::math::PI));
  const real_t amp = static_cast<real_t>(8.0) / (four_pi * rho * r0 * r0 * sigma * sqrt_2pi);

  const auto [x_bins, y_vals] = generateGaussianProfile(r0, sigma, amp, 0.0, 5.0, 0.005);

  const auto result = PeakIntegrationService::integrateRange(
      x_bins, y_vals, r0 - 3.5 * sigma, r0 + 3.5 * sigma, rho, /*is_radial=*/true);

  ASSERT_TRUE(result.coordination_number.has_value());
  EXPECT_NEAR(*result.coordination_number, 8.0, 0.1);
}

TEST_F(PeakIntegrationTests, AutoDetectPeakNearFindsTargetPeak) {
  // Two distinct peaks at r = 1.6 Å and r = 3.2 Å
  std::vector<real_t> x_bins;
  std::vector<real_t> y_vals;
  for (int i = 0; i <= 500; ++i) {
    const real_t x = static_cast<real_t>(i) * 0.01;
    x_bins.push_back(x);
    const real_t d1 = (x - 1.6) / 0.1;
    const real_t d2 = (x - 3.2) / 0.15;
    const real_t y = 3.0 * std::exp(-0.5 * d1 * d1) + 2.0 * std::exp(-0.5 * d2 * d2);
    y_vals.push_back(y);
  }

  // Click near first peak (e.g. at 1.55)
  const auto peak1 = PeakIntegrationService::autoDetectPeakNear(x_bins, y_vals, 1.55);
  ASSERT_TRUE(peak1.has_value());
  EXPECT_NEAR(peak1->r_peak, 1.6, 0.02);
  EXPECT_LT(peak1->r_min, 1.6);
  EXPECT_GT(peak1->r_max, 1.6);

  // Click near second peak (e.g. at 3.25)
  const auto peak2 = PeakIntegrationService::autoDetectPeakNear(x_bins, y_vals, 3.25);
  ASSERT_TRUE(peak2.has_value());
  EXPECT_NEAR(peak2->r_peak, 3.2, 0.02);
  EXPECT_LT(peak2->r_min, 3.2);
  EXPECT_GT(peak2->r_max, 3.2);
}

TEST_F(PeakIntegrationTests, HandlesEmptyAndSinglePointGracefully) {
  const std::vector<real_t> empty_bins;
  const std::vector<real_t> empty_vals;

  const auto res_empty = PeakIntegrationService::integrateRange(empty_bins, empty_vals, 0.0, 1.0);
  EXPECT_DOUBLE_EQ(res_empty.area, 0.0);
  EXPECT_FALSE(res_empty.coordination_number.has_value());

  const std::vector<real_t> single_bin = {1.0};
  const std::vector<real_t> single_val = {2.0};
  const auto res_single = PeakIntegrationService::integrateRange(single_bin, single_val, 0.5, 1.5);
  EXPECT_DOUBLE_EQ(res_single.area, 0.0);

  const auto auto_none = PeakIntegrationService::autoDetectPeakNear(empty_bins, empty_vals, 1.0);
  EXPECT_FALSE(auto_none.has_value());
}

TEST_F(PeakIntegrationTests, HandlesInvertedRangeGracefully) {
  const auto [x_bins, y_vals] = generateGaussianProfile(2.0, 0.1, 3.0);
  // Pass r_min = 2.5, r_max = 1.5 (inverted)
  const auto result = PeakIntegrationService::integrateRange(x_bins, y_vals, 2.5, 1.5);
  EXPECT_GT(result.area, 0.0);
  EXPECT_NEAR(result.r_peak, 2.0, 0.02);
}
