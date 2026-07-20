#include "math/Precision.hpp"
#include "math/SIMDUtils.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <random>
#include <vector>

namespace {
using correlation::real_t;

// Test fixture for SIMDUtils
class SIMDUtilsTests : public ::testing::Test {
public:
  // Random number generator setup for generating test data
  // NOLINTNEXTLINE(cert-msc51-cpp,cert-msc32-c,bugprone-random-generator-seed)
  std::mt19937 gen{1337}; // Fixed seed for reproducibility
  std::uniform_real_distribution<double> dist_d{-10.0, 10.0};
  std::uniform_real_distribution<double> pos_dist_d{0.1, 15.0}; // For distances

  std::uniform_real_distribution<real_t> dist{-10.0, 10.0};
  std::uniform_real_distribution<real_t> pos_dist{0.1, 15.0}; // For distances

protected:
  // Helper to generate a vector of random real_t
  std::vector<real_t> generateRandomData(size_t size, bool positive_only = false) {
    std::vector<real_t> data(size);
    for (size_t i = 0; i < size; ++i) {
      data[i] = positive_only ? pos_dist(gen) : dist(gen);
    }
    return data;
  }

  // Helper to generate a vector of random doubles
  std::vector<double> generateRandomDataDouble(size_t size, bool positive_only = false) {
    std::vector<double> data(size);
    for (size_t i = 0; i < size; ++i) {
      data[i] = positive_only ? pos_dist_d(gen) : dist_d(gen);
    }
    return data;
  }
};
} // namespace

// -----------------------------------------------------------------------------
// Test: sinc_integral (Fourier Transform core)
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, SincIntegralMatchesScalar) {
  // Test both exact multiples of SIMD width and non-multiples to hit the tail
  // logic
  const std::vector<size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100, 1024, 1025};
  const double Q_value = 2.5;

  for (size_t const size : sizes) {
    std::vector<double> integrand = generateRandomDataDouble(size);
    std::vector<double> rbins = generateRandomDataDouble(size, true); // r must be positive
    std::vector<double> scratch(size, 0.0);

    // Calculate expected result using the scalar fallback logic
    double expected_acc = 0.0;
    for (size_t j = 0; j < size; ++j) {
      expected_acc += integrand[j] * std::sin(Q_value * rbins[j]);
    }

    // Call the SIMD implementation
    double const actual_acc =
        correlation::math::sinc_integral(Q_value, integrand.data(), rbins.data(), scratch.data(), size);

    // Assert with a small tolerance due to potential floating point reordering
    // in SIMD
    EXPECT_NEAR(actual_acc, expected_acc, 1e-5) << "Failed for size: " << size;
  }
}

// -----------------------------------------------------------------------------
// Test: simd_dot (Dot Product core)
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, SimdDotMatchesScalar) {
  const std::vector<size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100, 1024, 1025};

  for (size_t const size : sizes) {
    std::vector<double> vec_a = generateRandomDataDouble(size);
    std::vector<double> vec_b = generateRandomDataDouble(size);

    // Calculate expected result using the scalar fallback logic
    double expected_acc = 0.0;
    for (size_t j = 0; j < size; ++j) {
      expected_acc += vec_a[j] * vec_b[j];
    }

    // Call the SIMD implementation
    double const actual_acc = correlation::math::simd_dot(vec_a.data(), vec_b.data(), size);

    // Assert with a small tolerance due to potential floating point reordering in SIMD
    EXPECT_NEAR(actual_acc, expected_acc, 1e-5) << "Failed for size: " << size;
  }
}

// -----------------------------------------------------------------------------
// Test: compute_dsq_block (Distance Squared)
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, ComputeDsqBlockMatchesScalar) {
  const std::vector<size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100, 1024, 1025};

  // Single atom A
  real_t const vec_ax = dist(gen);
  real_t const vec_ay = dist(gen);
  real_t const vec_az = dist(gen);

  for (size_t const size : sizes) {
    // Block of atoms B
    std::vector<real_t> vec_bx = generateRandomData(size);
    std::vector<real_t> vec_by = generateRandomData(size);
    std::vector<real_t> vec_bz = generateRandomData(size);

    correlation::math::PositionBlock const block{
        .x = vec_bx.data(), .y = vec_by.data(), .z = vec_bz.data(), .count = size};

    // Output array, initialized to -1 to detect unwritten values
    std::vector<real_t> actual_dsq(size, static_cast<real_t>(-1.0));

    // Calculate expected scalar result
    std::vector<real_t> expected_dsq(size);
    for (size_t k = 0; k < size; ++k) {
      expected_dsq[k] = correlation::math::dist_sq_scalar({.x = vec_ax, .y = vec_ay, .z = vec_az},
                                                          {.x = vec_bx[k], .y = vec_by[k], .z = vec_bz[k]});
    }

    // Call SIMD implementation
    correlation::math::compute_dsq_block(vec_ax, vec_ay, vec_az, block, actual_dsq.data());

    for (size_t k = 0; k < size; ++k) {
      EXPECT_NEAR(actual_dsq[k], expected_dsq[k], 1e-4) << "Failed at index " << k << " for size: " << size;
    }
  }
}

// -----------------------------------------------------------------------------
// Test: normalize_rdf_bins
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, NormalizeRDFBinsMatchesScalar) {
  const std::vector<size_t> sizes = {2, 5, 8, 9, 16, 17, 33, 100}; // Need at least size 2 since bin 0 is skipped

  const double g_norm = 1.5;
  const double inv_Ni_dr = 0.1;
  const double inv_Nj_dr = 0.2;
  const double pi4_rho_j = 3.14;

  for (size_t const size : sizes) {
    std::vector<double> vec_H = generateRandomDataDouble(size, true);
    std::vector<double> rbins = generateRandomDataDouble(size, true);
    // Ensure rbins[0] is very small/zero as happens in practice
    rbins[0] = 0.0;

    std::vector<double> actual_g(size, 0.0);
    std::vector<double> actual_G(size, 0.0);
    std::vector<double> actual_J(size, 0.0);
    std::vector<double> actual_Jinv(size, 0.0);

    // Call SIMD implementation
    correlation::math::normalize_rdf_bins(vec_H.data(), rbins.data(), g_norm, inv_Ni_dr, inv_Nj_dr, pi4_rho_j,
                                          actual_g.data(), actual_G.data(), actual_J.data(), actual_Jinv.data(), size);

    // Compare with scalar logic from index 1
    for (size_t k = 1; k < size; ++k) {
      const double r_val = rbins[k];
      if (r_val < 1e-9) {
        continue;
      }

      const double expected_g = vec_H[k] * g_norm / (r_val * r_val);
      const double expected_G = pi4_rho_j * r_val * (expected_g - 1.0);
      const double expected_J = vec_H[k] * inv_Ni_dr;
      const double expected_Jinv = vec_H[k] * inv_Nj_dr;

      EXPECT_NEAR(actual_g[k], expected_g, 1e-9) << "g failed at index " << k << " for size: " << size;
      EXPECT_NEAR(actual_G[k], expected_G, 1e-9) << "G failed at index " << k << " for size: " << size;
      EXPECT_NEAR(actual_J[k], expected_J, 1e-9) << "J failed at index " << k << " for size: " << size;
      EXPECT_NEAR(actual_Jinv[k], expected_Jinv, 1e-9) << "Jinv failed at index " << k << " for size: " << size;
    }
  }
}

// -----------------------------------------------------------------------------
// Test: scale_bins
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, ScaleBinsMatchesScalar) {
  const std::vector<size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100};
  const double scale_factor = 2.5;

  for (size_t const size : sizes) {
    std::vector<double> data = generateRandomDataDouble(size);
    std::vector<double> expected = data;

    for (auto &val : expected) {
      val *= scale_factor;
    }

    correlation::math::scale_bins(data.data(), scale_factor, size);

    for (size_t k = 0; k < size; ++k) {
      EXPECT_NEAR(data[k], expected[k], 1e-9) << "Failed at index " << k << " for size: " << size;
    }
  }
}

// -----------------------------------------------------------------------------
// Test: dot_block
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, DotBlockMatchesScalar) {
  const std::vector<size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100};

  double const v1x = dist_d(gen);
  double const v1y = dist_d(gen);
  double const v1z = dist_d(gen);

  for (size_t const size : sizes) {
    std::vector<double> v2x = generateRandomDataDouble(size);
    std::vector<double> v2y = generateRandomDataDouble(size);
    std::vector<double> v2z = generateRandomDataDouble(size);

    std::vector<double> actual_out(size, 0.0);

    std::vector<double> expected_out(size);
    for (size_t k = 0; k < size; ++k) {
      expected_out[k] = v1x * v2x[k] + v1y * v2y[k] + v1z * v2z[k];
    }

    correlation::math::dot_block(v1x, v1y, v1z, v2x.data(), v2y.data(), v2z.data(), actual_out.data(), size);

    for (size_t k = 0; k < size; ++k) {
      EXPECT_NEAR(actual_out[k], expected_out[k], 1e-9) << "Failed at index " << k << " for size: " << size;
    }
  }
}

// -----------------------------------------------------------------------------
// Test: Kahan Summation Precision
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, KahanSummationPrecision) {
  // Construct inputs where standard summation would truncate small contributions
  std::vector<double> vec_a = {1.0, 1e-16, 1e-16};
  std::vector<double> vec_b = {1.0, 1.0, 1.0};

  double std_sum = 1.0;
  std_sum += vec_a[1] * vec_b[1];
  std_sum += vec_a[2] * vec_b[2];

  double kahan_sum = correlation::math::simd_dot(vec_a.data(), vec_b.data(), 3);

  // Standard addition of 1e-16 to 1.0 rounds down to 1.0 immediately
  EXPECT_EQ(std_sum, 1.0);

  // Kahan compensated summation preserves the accumulated low-order bits,
  // yielding a result strictly greater than 1.0
  EXPECT_GT(kahan_sum, 1.0);
  EXPECT_DOUBLE_EQ(kahan_sum, 1.0000000000000002);
}
