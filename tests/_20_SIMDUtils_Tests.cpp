// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "math/SIMDUtils.hpp"
#include <cmath>
#include <gtest/gtest.h>
#include <random>
#include <vector>

// Test fixture for SIMDUtils
class _20_SIMDUtils_Tests : public ::testing::Test {
protected:
  // Random number generator setup for generating test data
  std::mt19937 gen{1337}; // Fixed seed for reproducibility
  std::uniform_real_distribution<double> dist{-10.0, 10.0};
  std::uniform_real_distribution<double> pos_dist{0.1, 15.0}; // For distances

  // Helper to generate a vector of random doubles
  std::vector<double> generateRandomData(size_t size,
                                         bool positive_only = false) {
    std::vector<double> data(size);
    for (size_t i = 0; i < size; ++i) {
      data[i] = positive_only ? pos_dist(gen) : dist(gen);
    }
    return data;
  }
};

// -----------------------------------------------------------------------------
// Test: sinc_integral (Fourier Transform core)
// -----------------------------------------------------------------------------
TEST_F(_20_SIMDUtils_Tests, SincIntegralMatchesScalar) {
  // Test both exact multiples of SIMD width and non-multiples to hit the tail
  // logic
  const std::vector<size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100, 1024, 1025};
  const double Q = 2.5;

  for (size_t size : sizes) {
    std::vector<double> integrand = generateRandomData(size);
    std::vector<double> rbins =
        generateRandomData(size, true); // r must be positive
    std::vector<double> scratch(size, 0.0);

    // Calculate expected result using the scalar fallback logic
    double expected_acc = 0.0;
    for (size_t j = 0; j < size; ++j) {
      expected_acc += integrand[j] * std::sin(Q * rbins[j]);
    }

    // Call the SIMD implementation
    double actual_acc = correlation::math::sinc_integral(
        Q, integrand.data(), rbins.data(), scratch.data(), size);

    // Assert with a small tolerance due to potential floating point reordering
    // in SIMD
    EXPECT_NEAR(actual_acc, expected_acc, 1e-9) << "Failed for size: " << size;
  }
}

// -----------------------------------------------------------------------------
// Test: compute_dsq_block (Distance Squared)
// -----------------------------------------------------------------------------
TEST_F(_20_SIMDUtils_Tests, ComputeDsqBlockMatchesScalar) {
  const std::vector<size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100, 1024, 1025};

  // Single atom A
  double ax = dist(gen);
  double ay = dist(gen);
  double az = dist(gen);

  for (size_t size : sizes) {
    // Block of atoms B
    std::vector<double> bx = generateRandomData(size);
    std::vector<double> by = generateRandomData(size);
    std::vector<double> bz = generateRandomData(size);

    correlation::math::PositionBlock block{bx.data(), by.data(), bz.data(),
                                           size};

    // Output array, initialized to -1 to detect unwritten values
    std::vector<double> actual_dsq(size, -1.0);

    // Calculate expected scalar result
    std::vector<double> expected_dsq(size);
    for (size_t k = 0; k < size; ++k) {
      expected_dsq[k] =
          correlation::math::dist_sq_scalar(ax, ay, az, bx[k], by[k], bz[k]);
    }

    // Call SIMD implementation
    correlation::math::compute_dsq_block(ax, ay, az, block, actual_dsq.data());

    for (size_t k = 0; k < size; ++k) {
      EXPECT_NEAR(actual_dsq[k], expected_dsq[k], 1e-9)
          << "Failed at index " << k << " for size: " << size;
    }
  }
}

// -----------------------------------------------------------------------------
// Test: normalize_rdf_bins
// -----------------------------------------------------------------------------
TEST_F(_20_SIMDUtils_Tests, NormalizeRDFBinsMatchesScalar) {
  const std::vector<size_t> sizes = {
      2,  5,  8,  9,
      16, 17, 33, 100}; // Need at least size 2 since bin 0 is skipped

  const double g_norm = 1.5;
  const double inv_Ni_dr = 0.1;
  const double inv_Nj_dr = 0.2;
  const double pi4_rho_j = 3.14;

  for (size_t size : sizes) {
    std::vector<double> H = generateRandomData(size, true);
    std::vector<double> rbins = generateRandomData(size, true);
    // Ensure rbins[0] is very small/zero as happens in practice
    rbins[0] = 0.0;

    std::vector<double> actual_g(size, 0.0);
    std::vector<double> actual_G(size, 0.0);
    std::vector<double> actual_J(size, 0.0);
    std::vector<double> actual_Jinv(size, 0.0);

    // Call SIMD implementation
    correlation::math::normalize_rdf_bins(
        H.data(), rbins.data(), g_norm, inv_Ni_dr, inv_Nj_dr, pi4_rho_j,
        actual_g.data(), actual_G.data(), actual_J.data(), actual_Jinv.data(),
        size);

    // Compare with scalar logic from index 1
    for (size_t k = 1; k < size; ++k) {
      const double r = rbins[k];
      if (r < 1e-9)
        continue;

      const double expected_g = H[k] * g_norm / (r * r);
      const double expected_G = pi4_rho_j * r * (expected_g - 1.0);
      const double expected_J = H[k] * inv_Ni_dr;
      const double expected_Jinv = H[k] * inv_Nj_dr;

      EXPECT_NEAR(actual_g[k], expected_g, 1e-9)
          << "g failed at index " << k << " for size: " << size;
      EXPECT_NEAR(actual_G[k], expected_G, 1e-9)
          << "G failed at index " << k << " for size: " << size;
      EXPECT_NEAR(actual_J[k], expected_J, 1e-9)
          << "J failed at index " << k << " for size: " << size;
      EXPECT_NEAR(actual_Jinv[k], expected_Jinv, 1e-9)
          << "Jinv failed at index " << k << " for size: " << size;
    }
  }
}

// -----------------------------------------------------------------------------
// Test: scale_bins
// -----------------------------------------------------------------------------
TEST_F(_20_SIMDUtils_Tests, ScaleBinsMatchesScalar) {
  const std::vector<size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100};
  const double scale_factor = 2.5;

  for (size_t size : sizes) {
    std::vector<double> data = generateRandomData(size);
    std::vector<double> expected = data;

    for (auto &val : expected) {
      val *= scale_factor;
    }

    correlation::math::scale_bins(data.data(), scale_factor, size);

    for (size_t k = 0; k < size; ++k) {
      EXPECT_NEAR(data[k], expected[k], 1e-9)
          << "Failed at index " << k << " for size: " << size;
    }
  }
}

// -----------------------------------------------------------------------------
// Test: dot_block
// -----------------------------------------------------------------------------
TEST_F(_20_SIMDUtils_Tests, DotBlockMatchesScalar) {
  const std::vector<size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100};

  double v1x = dist(gen);
  double v1y = dist(gen);
  double v1z = dist(gen);

  for (size_t size : sizes) {
    std::vector<double> v2x = generateRandomData(size);
    std::vector<double> v2y = generateRandomData(size);
    std::vector<double> v2z = generateRandomData(size);

    std::vector<double> actual_out(size, 0.0);

    std::vector<double> expected_out(size);
    for (size_t k = 0; k < size; ++k) {
      expected_out[k] = v1x * v2x[k] + v1y * v2y[k] + v1z * v2z[k];
    }

    correlation::math::dot_block(v1x, v1y, v1z, v2x.data(), v2y.data(),
                                 v2z.data(), actual_out.data(), size);

    for (size_t k = 0; k < size; ++k) {
      EXPECT_NEAR(actual_out[k], expected_out[k], 1e-9)
          << "Failed at index " << k << " for size: " << size;
    }
  }
}
