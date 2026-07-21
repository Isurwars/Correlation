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

  std::uniform_real_distribution<float> dist_f{-10.0F, 10.0F};
  std::uniform_real_distribution<float> pos_dist_f{0.1F, 15.0F};

  std::uniform_real_distribution<real_t> dist{-10.0, 10.0};
  std::uniform_real_distribution<real_t> pos_dist{0.1, 15.0}; // For distances

protected:
  // Helper to generate a vector of random real_t
  std::vector<real_t> generateRandomData(std::size_t size, bool positive_only = false) {
    std::vector<real_t> data(size);
    for (std::size_t idx = 0; idx < size; ++idx) {
      data[idx] = positive_only ? pos_dist(gen) : dist(gen);
    }
    return data;
  }

  // Helper to generate a vector of random doubles
  std::vector<double> generateRandomDataDouble(std::size_t size, bool positive_only = false) {
    std::vector<double> data(size);
    for (std::size_t idx = 0; idx < size; ++idx) {
      data[idx] = positive_only ? pos_dist_d(gen) : dist_d(gen);
    }
    return data;
  }

  // Helper to generate a vector of random floats
  std::vector<float> generateRandomDataFloat(std::size_t size, bool positive_only = false) {
    std::vector<float> data(size);
    for (std::size_t idx = 0; idx < size; ++idx) {
      data[idx] = positive_only ? pos_dist_f(gen) : dist_f(gen);
    }
    return data;
  }
};
} // namespace

// -----------------------------------------------------------------------------
// Test: sinc_integral (Fourier Transform core)
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, SincIntegralMatchesScalar) {
  const std::vector<std::size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100, 1024, 1025};
  const double q_value = 2.5;

  for (const std::size_t size : sizes) {
    std::vector<double> integrand = generateRandomDataDouble(size);
    std::vector<double> rbins = generateRandomDataDouble(size, true); // r must be positive
    std::vector<double> scratch(size, 0.0);

    double expected_acc = 0.0;
    for (std::size_t idx = 0; idx < size; ++idx) {
      expected_acc += integrand[idx] * std::sin(q_value * rbins[idx]);
    }

    const double actual_acc =
        correlation::math::sinc_integral(q_value, integrand.data(), rbins.data(), scratch.data(), size);

    EXPECT_NEAR(actual_acc, expected_acc, 1e-5) << "Failed for size: " << size;
  }
}

TEST_F(SIMDUtilsTests, SincIntegralFloatMatchesScalar) {
  const std::vector<std::size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100};
  const float q_value = 2.5F;

  for (const std::size_t size : sizes) {
    std::vector<float> integrand = generateRandomDataFloat(size);
    std::vector<float> rbins = generateRandomDataFloat(size, true);
    std::vector<float> scratch(size, 0.0F);

    float expected_acc = 0.0F;
    for (std::size_t idx = 0; idx < size; ++idx) {
      expected_acc += integrand[idx] * std::sin(q_value * rbins[idx]);
    }

    const float actual_acc =
        correlation::math::sinc_integral(q_value, integrand.data(), rbins.data(), scratch.data(), size);

    EXPECT_NEAR(actual_acc, expected_acc, 1e-3F) << "Float sinc_integral failed for size: " << size;
  }
}

// -----------------------------------------------------------------------------
// Test: simd_dot (Dot Product core)
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, SimdDotMatchesScalar) {
  const std::vector<std::size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100, 1024, 1025};

  for (const std::size_t size : sizes) {
    std::vector<double> vec_a = generateRandomDataDouble(size);
    std::vector<double> vec_b = generateRandomDataDouble(size);

    double expected_acc = 0.0;
    for (std::size_t idx = 0; idx < size; ++idx) {
      expected_acc += vec_a[idx] * vec_b[idx];
    }

    const double actual_acc = correlation::math::simd_dot(vec_a.data(), vec_b.data(), size);

    EXPECT_NEAR(actual_acc, expected_acc, 1e-5) << "Failed for size: " << size;
  }
}

TEST_F(SIMDUtilsTests, SimdDotFloatMatchesScalar) {
  const std::vector<std::size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100};

  for (const std::size_t size : sizes) {
    std::vector<float> vec_a = generateRandomDataFloat(size);
    std::vector<float> vec_b = generateRandomDataFloat(size);

    float expected_acc = 0.0F;
    for (std::size_t idx = 0; idx < size; ++idx) {
      expected_acc += vec_a[idx] * vec_b[idx];
    }

    const float actual_acc = correlation::math::simd_dot(vec_a.data(), vec_b.data(), size);

    EXPECT_NEAR(actual_acc, expected_acc, 1e-3F) << "Float simd_dot failed for size: " << size;
  }
}

// -----------------------------------------------------------------------------
// Test: compute_dsq_block (Distance Squared)
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, ComputeDsqBlockMatchesScalar) {
  const std::vector<std::size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100, 1024, 1025};

  const real_t vec_ax = dist(gen);
  const real_t vec_ay = dist(gen);
  const real_t vec_az = dist(gen);

  for (const std::size_t size : sizes) {
    std::vector<real_t> vec_bx = generateRandomData(size);
    std::vector<real_t> vec_by = generateRandomData(size);
    std::vector<real_t> vec_bz = generateRandomData(size);

    const correlation::math::PositionBlock block{
        .x = vec_bx.data(), .y = vec_by.data(), .z = vec_bz.data(), .count = size};

    std::vector<real_t> actual_dsq(size, static_cast<real_t>(-1.0));

    std::vector<real_t> expected_dsq(size);
    for (std::size_t idx = 0; idx < size; ++idx) {
      expected_dsq[idx] = correlation::math::dist_sq_scalar({.x = vec_ax, .y = vec_ay, .z = vec_az},
                                                            {.x = vec_bx[idx], .y = vec_by[idx], .z = vec_bz[idx]});
    }

    correlation::math::compute_dsq_block(vec_ax, vec_ay, vec_az, block, actual_dsq.data());

    for (std::size_t idx = 0; idx < size; ++idx) {
      EXPECT_NEAR(actual_dsq[idx], expected_dsq[idx], 1e-4) << "Failed at index " << idx << " for size: " << size;
    }
  }
}

TEST_F(SIMDUtilsTests, ComputeDsqBlockFloatMatchesScalar) {
  const std::vector<std::size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100};

  const float vec_ax = dist_f(gen);
  const float vec_ay = dist_f(gen);
  const float vec_az = dist_f(gen);

  for (const std::size_t size : sizes) {
    std::vector<float> vec_bx = generateRandomDataFloat(size);
    std::vector<float> vec_by = generateRandomDataFloat(size);
    std::vector<float> vec_bz = generateRandomDataFloat(size);

    const correlation::math::PositionBlockT<float> block{
        .x = vec_bx.data(), .y = vec_by.data(), .z = vec_bz.data(), .count = size};

    std::vector<float> actual_dsq(size, -1.0F);
    std::vector<float> expected_dsq(size);
    for (std::size_t idx = 0; idx < size; ++idx) {
      expected_dsq[idx] = correlation::math::dist_sq_scalar<float>({.x = vec_ax, .y = vec_ay, .z = vec_az},
                                                                   {.x = vec_bx[idx], .y = vec_by[idx], .z = vec_bz[idx]});
    }

    correlation::math::compute_dsq_block(vec_ax, vec_ay, vec_az, block, actual_dsq.data());

    for (std::size_t idx = 0; idx < size; ++idx) {
      EXPECT_NEAR(actual_dsq[idx], expected_dsq[idx], 1e-3F) << "Float dsq failed at index " << idx << " for size: " << size;
    }
  }
}

// -----------------------------------------------------------------------------
// Test: normalize_rdf_bins
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, NormalizeRDFBinsMatchesScalar) {
  const std::vector<std::size_t> sizes = {2, 5, 8, 9, 16, 17, 33, 100};

  const double g_norm = 1.5;
  const double inv_Ni_dr = 0.1;
  const double inv_Nj_dr = 0.2;
  const double pi4_rho_j = 3.14;

  for (const std::size_t size : sizes) {
    std::vector<double> vec_h = generateRandomDataDouble(size, true);
    std::vector<double> rbins = generateRandomDataDouble(size, true);
    rbins[0] = 0.0;

    std::vector<double> actual_g(size, 0.0);
    std::vector<double> actual_g_cap(size, 0.0);
    std::vector<double> actual_j(size, 0.0);
    std::vector<double> actual_jinv(size, 0.0);

    correlation::math::normalize_rdf_bins(vec_h.data(), rbins.data(), g_norm, inv_Ni_dr, inv_Nj_dr, pi4_rho_j,
                                          actual_g.data(), actual_g_cap.data(), actual_j.data(), actual_jinv.data(), size);

    for (std::size_t idx = 1; idx < size; ++idx) {
      const double r_val = rbins[idx];
      if (r_val < 1e-9) {
        continue;
      }

      const double expected_g = vec_h[idx] * g_norm / (r_val * r_val);
      const double expected_g_cap = pi4_rho_j * r_val * (expected_g - 1.0);
      const double expected_j = vec_h[idx] * inv_Ni_dr;
      const double expected_jinv = vec_h[idx] * inv_Nj_dr;

      EXPECT_NEAR(actual_g[idx], expected_g, 1e-9) << "g failed at index " << idx << " for size: " << size;
      EXPECT_NEAR(actual_g_cap[idx], expected_g_cap, 1e-9) << "G failed at index " << idx << " for size: " << size;
      EXPECT_NEAR(actual_j[idx], expected_j, 1e-9) << "J failed at index " << idx << " for size: " << size;
      EXPECT_NEAR(actual_jinv[idx], expected_jinv, 1e-9) << "Jinv failed at index " << idx << " for size: " << size;
    }
  }
}

TEST_F(SIMDUtilsTests, NormalizeRDFBinsFloatMatchesScalar) {
  const std::vector<std::size_t> sizes = {2, 5, 8, 9, 16, 17, 33, 100};

  const float g_norm = 1.5F;
  const float inv_Ni_dr = 0.1F;
  const float inv_Nj_dr = 0.2F;
  const float pi4_rho_j = 3.14F;

  for (const std::size_t size : sizes) {
    std::vector<float> vec_h = generateRandomDataFloat(size, true);
    std::vector<float> rbins = generateRandomDataFloat(size, true);
    rbins[0] = 0.0F;

    std::vector<float> actual_g(size, 0.0F);
    std::vector<float> actual_g_cap(size, 0.0F);
    std::vector<float> actual_j(size, 0.0F);
    std::vector<float> actual_jinv(size, 0.0F);

    correlation::math::normalize_rdf_bins(vec_h.data(), rbins.data(), g_norm, inv_Ni_dr, inv_Nj_dr, pi4_rho_j,
                                          actual_g.data(), actual_g_cap.data(), actual_j.data(), actual_jinv.data(), size);

    for (std::size_t idx = 1; idx < size; ++idx) {
      const float r_val = rbins[idx];
      if (r_val < 1e-9F) {
        continue;
      }

      const float expected_g = vec_h[idx] * g_norm / (r_val * r_val);
      const float expected_g_cap = pi4_rho_j * r_val * (expected_g - 1.0F);
      const float expected_j = vec_h[idx] * inv_Ni_dr;
      const float expected_jinv = vec_h[idx] * inv_Nj_dr;

      EXPECT_NEAR(actual_g[idx], expected_g, 1e-4F) << "Float g failed at index " << idx << " for size: " << size;
      EXPECT_NEAR(actual_g_cap[idx], expected_g_cap, 1e-4F) << "Float G failed at index " << idx << " for size: " << size;
      EXPECT_NEAR(actual_j[idx], expected_j, 1e-4F) << "Float J failed at index " << idx << " for size: " << size;
      EXPECT_NEAR(actual_jinv[idx], expected_jinv, 1e-4F) << "Float Jinv failed at index " << idx << " for size: " << size;
    }
  }
}

// -----------------------------------------------------------------------------
// Test: scale_bins
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, ScaleBinsMatchesScalar) {
  const std::vector<std::size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100};
  const double scale_factor = 2.5;

  for (const std::size_t size : sizes) {
    std::vector<double> data = generateRandomDataDouble(size);
    std::vector<double> expected = data;

    for (auto &val : expected) {
      val *= scale_factor;
    }

    correlation::math::scale_bins(data.data(), scale_factor, size);

    for (std::size_t idx = 0; idx < size; ++idx) {
      EXPECT_NEAR(data[idx], expected[idx], 1e-9) << "Failed at index " << idx << " for size: " << size;
    }
  }
}

TEST_F(SIMDUtilsTests, ScaleBinsFloatMatchesScalar) {
  const std::vector<std::size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100};
  const float scale_factor = 2.5F;

  for (const std::size_t size : sizes) {
    std::vector<float> data = generateRandomDataFloat(size);
    std::vector<float> expected = data;

    for (auto &val : expected) {
      val *= scale_factor;
    }

    correlation::math::scale_bins(data.data(), scale_factor, size);

    for (std::size_t idx = 0; idx < size; ++idx) {
      EXPECT_NEAR(data[idx], expected[idx], 1e-4F) << "Float scale_bins failed at index " << idx << " for size: " << size;
    }
  }
}

// -----------------------------------------------------------------------------
// Test: debye_sum
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, DebyeSumDoubleAndFloat) {
  const std::vector<std::size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100};
  const double q_d = 1.2;
  const float q_f = 1.2F;

  for (const std::size_t size : sizes) {
    std::vector<double> dists_d = generateRandomDataDouble(size, true);
    std::vector<double> scratch_d(size, 0.0);
    const double actual_d = correlation::math::debye_sum(q_d, dists_d.data(), scratch_d.data(), size);

    double expected_d = 0.0;
    for (std::size_t idx = 0; idx < size; ++idx) {
      const double val_x = q_d * dists_d[idx];
      expected_d += (val_x < 1e-4) ? (1.0 - (val_x * val_x) / 6.0) : (std::sin(val_x) / val_x);
    }
    EXPECT_NEAR(actual_d, expected_d, 1e-5) << "Debye sum double failed for size: " << size;

    std::vector<float> dists_f = generateRandomDataFloat(size, true);
    std::vector<float> scratch_f(size, 0.0F);
    const float actual_f = correlation::math::debye_sum(q_f, dists_f.data(), scratch_f.data(), size);

    float expected_f = 0.0F;
    for (std::size_t idx = 0; idx < size; ++idx) {
      const float val_x = q_f * dists_f[idx];
      expected_f += (val_x < 1e-4F) ? (1.0F - (val_x * val_x) / 6.0F) : (std::sin(val_x) / val_x);
    }
    EXPECT_NEAR(actual_f, expected_f, 1e-3F) << "Debye sum float failed for size: " << size;
  }
}


// -----------------------------------------------------------------------------
// Test: dot_block
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, DotBlockMatchesScalar) {
  const std::vector<std::size_t> sizes = {1, 4, 7, 8, 15, 16, 33, 100};

  const double v1x = dist_d(gen);
  const double v1y = dist_d(gen);
  const double v1z = dist_d(gen);

  for (const std::size_t size : sizes) {
    std::vector<double> v2x = generateRandomDataDouble(size);
    std::vector<double> v2y = generateRandomDataDouble(size);
    std::vector<double> v2z = generateRandomDataDouble(size);

    std::vector<double> actual_out(size, 0.0);

    std::vector<double> expected_out(size);
    for (std::size_t idx = 0; idx < size; ++idx) {
      expected_out[idx] = v1x * v2x[idx] + v1y * v2y[idx] + v1z * v2z[idx];
    }

    correlation::math::dot_block(v1x, v1y, v1z, v2x.data(), v2y.data(), v2z.data(), actual_out.data(), size);

    for (std::size_t idx = 0; idx < size; ++idx) {
      EXPECT_NEAR(actual_out[idx], expected_out[idx], 1e-9) << "Failed at index " << idx << " for size: " << size;
    }
  }
}

// -----------------------------------------------------------------------------
// Test: Kahan Summation Precision
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, KahanSummationPrecision) {
  std::vector<double> vec_a = {1.0, 1e-16, 1e-16};
  std::vector<double> vec_b = {1.0, 1.0, 1.0};

  double std_sum = 1.0;
  std_sum += vec_a[1] * vec_b[1];
  std_sum += vec_a[2] * vec_b[2];

  const double kahan_sum = correlation::math::simd_dot(vec_a.data(), vec_b.data(), 3);

  EXPECT_EQ(std_sum, 1.0);
  EXPECT_GT(kahan_sum, 1.0);
  EXPECT_DOUBLE_EQ(kahan_sum, 1.0000000000000002);
}

