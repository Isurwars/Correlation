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
using correlation::real_t;

namespace {

struct GaussianProfileConfig {
  real_t mean{static_cast<real_t>(2.5)};
  real_t sigma{static_cast<real_t>(0.15)};
  real_t amp{static_cast<real_t>(4.0)};
  real_t x_min{static_cast<real_t>(0.0)};
  real_t x_max{static_cast<real_t>(5.0)};
  real_t step_size{static_cast<real_t>(0.01)};
};

class PeakIntegrationTests : public ::testing::Test {
protected:
  static auto generateGaussianProfile(const GaussianProfileConfig &config = {})
      -> std::pair<std::vector<real_t>, std::vector<real_t>> {
    std::vector<real_t> x_bins;
    std::vector<real_t> y_vals;
    const auto n_pts = static_cast<std::size_t>((config.x_max - config.x_min) / config.step_size);
    x_bins.reserve(n_pts);
    y_vals.reserve(n_pts);

    for (std::size_t idx = 0; idx < n_pts; ++idx) {
      const real_t x_pos = config.x_min + static_cast<real_t>(idx) * config.step_size;
      x_bins.push_back(x_pos);
      const real_t diff = (x_pos - config.mean) / config.sigma;
      const real_t y_val = config.amp * std::exp(static_cast<real_t>(-0.5) * diff * diff);
      y_vals.push_back(y_val);
    }
    return {x_bins, y_vals};
  }
};

void verifyDetectedPeak(const std::optional<PeakIntegralResult> &peak, real_t expected_center) {
  ASSERT_TRUE(peak.has_value());
  EXPECT_NEAR(peak->r_peak, expected_center, static_cast<real_t>(0.02));
  EXPECT_LT(peak->r_min, expected_center);
  EXPECT_GT(peak->r_max, expected_center);
}

TEST_F(PeakIntegrationTests, IntegratesGaussianPeakAccurately) {
  constexpr real_t MEAN = 2.5;
  constexpr real_t SIGMA = 0.12;
  constexpr real_t AMP = 5.0;
  const auto [x_bins, y_vals] = generateGaussianProfile({.mean = MEAN, .sigma = SIGMA, .amp = AMP});

  const auto result = PeakIntegrationService::integrateRange(
      x_bins, y_vals, MEAN - static_cast<real_t>(4.0) * SIGMA,
      MEAN + static_cast<real_t>(4.0) * SIGMA);

  const real_t analytical_area =
      AMP * SIGMA * std::sqrt(static_cast<real_t>(2.0 * correlation::math::PI));
  EXPECT_NEAR(result.area, analytical_area, 0.02);
  EXPECT_NEAR(result.centroid, MEAN, 0.01);
  EXPECT_NEAR(result.r_peak, MEAN, 0.01);
  EXPECT_NEAR(result.peak_height, AMP, 0.01);

  // FWHM of Gaussian: 2 * sqrt(2 * ln(2)) * sigma ≈ 2.35482 * sigma
  const real_t expected_fwhm = static_cast<real_t>(2.35482) * SIGMA;
  EXPECT_NEAR(result.fwhm, expected_fwhm, 0.01);
}

TEST_F(PeakIntegrationTests, ComputesCoordinationNumberForRadialShell) {
  // Model a single coordination shell of 8 neighbors at r = 2.48 Å in a BCC cell
  constexpr real_t SHELL_RADIUS = 2.48;
  constexpr real_t SIGMA = 0.08;
  constexpr real_t RHO = 0.085; // macroscopic number density (atoms / Å³)

  // In radial shell: CN = 4πρ ∫ r² g(r) dr.
  // If g(r) = A * exp(-(r-r0)² / 2σ²), then ∫ r² g(r) dr ≈ r0² A σ sqrt(2π).
  // For CN = 8.0: A = 8.0 / (4πρ r0² σ sqrt(2π))
  constexpr real_t FOUR_PI = correlation::math::FOUR_PI;
  const auto sqrt_2pi = static_cast<real_t>(std::sqrt(2.0 * correlation::math::PI));
  const auto amp =
      static_cast<real_t>(8.0) / (FOUR_PI * RHO * SHELL_RADIUS * SHELL_RADIUS * SIGMA * sqrt_2pi);

  const auto [x_bins, y_vals] = generateGaussianProfile({
      .mean = SHELL_RADIUS,
      .sigma = SIGMA,
      .amp = amp,
      .x_min = static_cast<real_t>(0.0),
      .x_max = static_cast<real_t>(5.0),
      .step_size = static_cast<real_t>(0.005),
  });

  const auto result = PeakIntegrationService::integrateRange(
      x_bins, y_vals, SHELL_RADIUS - static_cast<real_t>(3.5) * SIGMA,
      SHELL_RADIUS + static_cast<real_t>(3.5) * SIGMA, RHO, /*is_radial=*/true);

  ASSERT_TRUE(result.coordination_number.has_value());
  EXPECT_NEAR(*result.coordination_number, 8.0, 0.1);
}

TEST_F(PeakIntegrationTests, AutoDetectPeakNearFindsTargetPeak) {
  // Two distinct peaks at r = 1.6 Å and r = 3.2 Å
  std::vector<real_t> x_bins;
  std::vector<real_t> y_vals;
  for (int idx = 0; idx <= 500; ++idx) {
    const real_t x_pos = static_cast<real_t>(idx) * static_cast<real_t>(0.01);
    x_bins.push_back(x_pos);
    const real_t diff1 = (x_pos - static_cast<real_t>(1.6)) / static_cast<real_t>(0.1);
    const real_t diff2 = (x_pos - static_cast<real_t>(3.2)) / static_cast<real_t>(0.15);
    const real_t y_val =
        static_cast<real_t>(3.0) * std::exp(static_cast<real_t>(-0.5) * diff1 * diff1) +
        static_cast<real_t>(2.0) * std::exp(static_cast<real_t>(-0.5) * diff2 * diff2);
    y_vals.push_back(y_val);
  }

  // Click near first peak (e.g. at 1.55)
  const auto peak1 =
      PeakIntegrationService::autoDetectPeakNear(x_bins, y_vals, static_cast<real_t>(1.55));
  verifyDetectedPeak(peak1, static_cast<real_t>(1.6));

  // Click near second peak (e.g. at 3.25)
  const auto peak2 =
      PeakIntegrationService::autoDetectPeakNear(x_bins, y_vals, static_cast<real_t>(3.25));
  verifyDetectedPeak(peak2, static_cast<real_t>(3.2));
}

TEST_F(PeakIntegrationTests, HandlesEmptyAndSinglePointGracefully) {
  const std::vector<real_t> empty_bins;
  const std::vector<real_t> empty_vals;

  const auto res_empty = PeakIntegrationService::integrateRange(
      empty_bins, empty_vals, static_cast<real_t>(0.0), static_cast<real_t>(1.0));
  EXPECT_DOUBLE_EQ(res_empty.area, 0.0);
  EXPECT_FALSE(res_empty.coordination_number.has_value());

  const std::vector<real_t> single_bin = {static_cast<real_t>(1.0)};
  const std::vector<real_t> single_val = {static_cast<real_t>(2.0)};
  const auto res_single = PeakIntegrationService::integrateRange(
      single_bin, single_val, static_cast<real_t>(0.5), static_cast<real_t>(1.5));
  EXPECT_DOUBLE_EQ(res_single.area, 0.0);

  const auto auto_none =
      PeakIntegrationService::autoDetectPeakNear(empty_bins, empty_vals, static_cast<real_t>(1.0));
  EXPECT_FALSE(auto_none.has_value());
}

TEST_F(PeakIntegrationTests, HandlesInvertedRangeGracefully) {
  const auto [x_bins, y_vals] = generateGaussianProfile({
      .mean = static_cast<real_t>(2.0),
      .sigma = static_cast<real_t>(0.1),
      .amp = static_cast<real_t>(3.0),
  });
  // Pass r_min = 2.5, r_max = 1.5 (inverted)
  const auto result = PeakIntegrationService::integrateRange(
      x_bins, y_vals, static_cast<real_t>(2.5), static_cast<real_t>(1.5));
  EXPECT_GT(result.area, 0.0);
  EXPECT_NEAR(result.r_peak, 2.0, 0.02);
}

} // namespace
