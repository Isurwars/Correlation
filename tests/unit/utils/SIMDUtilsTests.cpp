#include "math/Precision.hpp"
#include "math/SIMDUtils.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <random>
#include <vector>

namespace {
using correlation::real_t;

struct DummyPos {
  double x_val;
  double y_val;
  double z_val;
  [[nodiscard]] DummyPos position() const { return *this; }
  [[nodiscard]] double x() const { return x_val; }
  [[nodiscard]] double y() const { return y_val; }
  [[nodiscard]] double z() const { return z_val; }
};

// Test fixture for SIMDUtils
class SIMDUtilsTests : public ::testing::Test {
public:
  // Random number generator setup for generating test data
  std::seed_seq seed{1337};
  std::mt19937 gen{seed};
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

    const double actual_acc = correlation::math::sinc_integral(correlation::math::SincIntegralParams<double>{
        .q_magnitude = q_value,
        .integrand = integrand.data(),
        .radial_bins = rbins.data(),
        .sinqr_scratch = scratch.data(),
        .count = size,
    });

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

    const float actual_acc = correlation::math::sinc_integral(correlation::math::SincIntegralParams<float>{
        .q_magnitude = q_value,
        .integrand = integrand.data(),
        .radial_bins = rbins.data(),
        .sinqr_scratch = scratch.data(),
        .count = size,
    });

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
        .x = vec_bx.data(),
        .y = vec_by.data(),
        .z = vec_bz.data(),
        .count = size,
    };

    std::vector<real_t> actual_dsq(size, static_cast<real_t>(-1.0));

    std::vector<real_t> expected_dsq(size);
    for (std::size_t idx = 0; idx < size; ++idx) {
      expected_dsq[idx] = correlation::math::dist_sq_scalar(
          {
              .x = vec_ax,
              .y = vec_ay,
              .z = vec_az,
          },
          {
              .x = vec_bx[idx],
              .y = vec_by[idx],
              .z = vec_bz[idx],
          });
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
        .x = vec_bx.data(),
        .y = vec_by.data(),
        .z = vec_bz.data(),
        .count = size,
    };

    std::vector<float> actual_dsq(size, -1.0F);
    std::vector<float> expected_dsq(size);
    for (std::size_t idx = 0; idx < size; ++idx) {
      expected_dsq[idx] = correlation::math::dist_sq_scalar<float>(
          {
              .x = vec_ax,
              .y = vec_ay,
              .z = vec_az,
          },
          {
              .x = vec_bx[idx],
              .y = vec_by[idx],
              .z = vec_bz[idx],
          });
    }

    correlation::math::compute_dsq_block(vec_ax, vec_ay, vec_az, block, actual_dsq.data());

    for (std::size_t idx = 0; idx < size; ++idx) {
      EXPECT_NEAR(actual_dsq[idx], expected_dsq[idx], 1e-3F)
          << "Float dsq failed at index " << idx << " for size: " << size;
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

    correlation::math::normalize_rdf_bins(correlation::math::RDFNormalizationParams<double>{
        .hist_data = vec_h.data(),
        .radial_bins = rbins.data(),
        .g_norm = g_norm,
        .inv_Ni_dr = inv_Ni_dr,
        .inv_Nj_dr = inv_Nj_dr,
        .pi4_rho_j = pi4_rho_j,
        .g_out = actual_g.data(),
        .G_out = actual_g_cap.data(),
        .J_out = actual_j.data(),
        .Jinv_out = actual_jinv.data(),
        .count = size,
    });

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

    correlation::math::normalize_rdf_bins(correlation::math::RDFNormalizationParams<float>{
        .hist_data = vec_h.data(),
        .radial_bins = rbins.data(),
        .g_norm = g_norm,
        .inv_Ni_dr = inv_Ni_dr,
        .inv_Nj_dr = inv_Nj_dr,
        .pi4_rho_j = pi4_rho_j,
        .g_out = actual_g.data(),
        .G_out = actual_g_cap.data(),
        .J_out = actual_j.data(),
        .Jinv_out = actual_jinv.data(),
        .count = size,
    });

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
      EXPECT_NEAR(actual_g_cap[idx], expected_g_cap, 1e-4F)
          << "Float G failed at index " << idx << " for size: " << size;
      EXPECT_NEAR(actual_j[idx], expected_j, 1e-4F) << "Float J failed at index " << idx << " for size: " << size;
      EXPECT_NEAR(actual_jinv[idx], expected_jinv, 1e-4F)
          << "Float Jinv failed at index " << idx << " for size: " << size;
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
      EXPECT_NEAR(data[idx], expected[idx], 1e-4F)
          << "Float scale_bins failed at index " << idx << " for size: " << size;
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

// -----------------------------------------------------------------------------
// Test: normalize_rdf_bins Bin 0 Initialization Bug Fix
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, NormalizeRdfBinsZeroInitializesBinZeroDouble) {
  constexpr std::size_t size = 16;
  std::vector<double> hist(size, 5.0);
  std::vector<double> rbins(size);
  for (std::size_t i = 0; i < size; ++i) {
    rbins[i] = static_cast<double>(i) * 0.1;
  }

  // Pre-fill outputs with non-zero garbage to ensure bin 0 is explicitly zeroed
  std::vector<double> g_out(size, 999.0);
  std::vector<double> G_out(size, 999.0);
  std::vector<double> J_out(size, 999.0);
  std::vector<double> Jinv_out(size, 999.0);

  correlation::math::normalize_rdf_bins(correlation::math::RDFNormalizationParams<double>{
      .hist_data = hist.data(),
      .radial_bins = rbins.data(),
      .g_norm = 1.0,
      .inv_Ni_dr = 0.5,
      .inv_Nj_dr = 0.5,
      .pi4_rho_j = 12.56,
      .g_out = g_out.data(),
      .G_out = G_out.data(),
      .J_out = J_out.data(),
      .Jinv_out = Jinv_out.data(),
      .count = size,
  });

  EXPECT_DOUBLE_EQ(g_out[0], 0.0);
  EXPECT_DOUBLE_EQ(G_out[0], 0.0);
  EXPECT_DOUBLE_EQ(J_out[0], 0.0);
  EXPECT_DOUBLE_EQ(Jinv_out[0], 0.0);
}

TEST_F(SIMDUtilsTests, NormalizeRdfBinsZeroInitializesBinZeroFloat) {
  constexpr std::size_t size = 16;
  std::vector<float> hist(size, 5.0F);
  std::vector<float> rbins(size);
  for (std::size_t i = 0; i < size; ++i) {
    rbins[i] = static_cast<float>(i) * 0.1F;
  }

  // Pre-fill outputs with non-zero garbage
  std::vector<float> g_out(size, 999.0F);
  std::vector<float> G_out(size, 999.0F);
  std::vector<float> J_out(size, 999.0F);
  std::vector<float> Jinv_out(size, 999.0F);

  correlation::math::normalize_rdf_bins(correlation::math::RDFNormalizationParams<float>{
      .hist_data = hist.data(),
      .radial_bins = rbins.data(),
      .g_norm = 1.0F,
      .inv_Ni_dr = 0.5F,
      .inv_Nj_dr = 0.5F,
      .pi4_rho_j = 12.56F,
      .g_out = g_out.data(),
      .G_out = G_out.data(),
      .J_out = J_out.data(),
      .Jinv_out = Jinv_out.data(),
      .count = size,
  });

  EXPECT_FLOAT_EQ(g_out[0], 0.0F);
  EXPECT_FLOAT_EQ(G_out[0], 0.0F);
  EXPECT_FLOAT_EQ(J_out[0], 0.0F);
  EXPECT_FLOAT_EQ(Jinv_out[0], 0.0F);
}

// -----------------------------------------------------------------------------
// Test: Template Dispatchers
// -----------------------------------------------------------------------------
TEST_F(SIMDUtilsTests, TemplateDispatchersMatchFloat) {
  const std::vector<float> vec_a = {1.0F, 2.0F, 3.0F, 4.0F};
  const std::vector<float> vec_b = {0.5F, 1.5F, 2.5F, 3.5F};
  const float dot = correlation::math::simd_dot<float>(vec_a.data(), vec_b.data(), vec_a.size());
  EXPECT_NEAR(dot, 25.0F, 1e-5F);
}

TEST_F(SIMDUtilsTests, MillerPhaseSumParamsStructMatchesScalar) {
  const std::vector<double> c_1 = {1.0, 0.5, 0.2};
  const std::vector<double> s_1 = {0.0, 0.866, 0.9798};
  const std::vector<double> c_2 = {0.8, 0.6, 1.0};
  const std::vector<double> s_2 = {0.6, 0.8, 0.0};
  const std::vector<double> c_3 = {0.5, 0.5, 0.5};
  const std::vector<double> s_3 = {0.866, 0.866, 0.866};

  correlation::math::MillerPhaseSumResult<double> res;
  correlation::math::miller_phase_sum(
      correlation::math::MillerPhaseSumParams<double>{
          .cos1 = c_1.data(),
          .sin1 = s_1.data(),
          .cos2 = c_2.data(),
          .sin2 = s_2.data(),
          .cos3 = c_3.data(),
          .sin3 = s_3.data(),
          .count = c_1.size(),
      },
      res);

  double expected_cos = 0.0;
  double expected_sin = 0.0;
  for (std::size_t i = 0; i < c_1.size(); ++i) {
    const double c12 = c_1[i] * c_2[i] - s_1[i] * s_2[i];
    const double s12 = s_1[i] * c_2[i] + c_1[i] * s_2[i];
    expected_cos += c12 * c_3[i] - s12 * s_3[i];
    expected_sin += s12 * c_3[i] + c12 * s_3[i];
  }

  EXPECT_NEAR(res.cos_sum, expected_cos, 1e-5);
  EXPECT_NEAR(res.sin_sum, expected_sin, 1e-5);
}

TEST_F(SIMDUtilsTests, MillerPhaseSumResultStructMatchesScalar) {
  const std::vector<double> c_1 = {1.0, 0.5, 0.2};
  const std::vector<double> s_1 = {0.0, 0.866, 0.9798};
  const std::vector<double> c_2 = {0.8, 0.6, 1.0};
  const std::vector<double> s_2 = {0.6, 0.8, 0.0};
  const std::vector<double> c_3 = {0.5, 0.5, 0.5};
  const std::vector<double> s_3 = {0.866, 0.866, 0.866};

  const auto params = correlation::math::MillerPhaseSumParams<double>{
      .cos1 = c_1.data(),
      .sin1 = s_1.data(),
      .cos2 = c_2.data(),
      .sin2 = s_2.data(),
      .cos3 = c_3.data(),
      .sin3 = s_3.data(),
      .count = c_1.size(),
  };

  const correlation::math::MillerPhaseSumResult<double> res = correlation::math::miller_phase_sum(params);

  double expected_cos = 0.0;
  double expected_sin = 0.0;
  for (std::size_t i = 0; i < c_1.size(); ++i) {
    const double c12 = c_1[i] * c_2[i] - s_1[i] * s_2[i];
    const double s12 = s_1[i] * c_2[i] + c_1[i] * s_2[i];
    expected_cos += c12 * c_3[i] - s12 * s_3[i];
    expected_sin += s12 * c_3[i] + c12 * s_3[i];
  }

  EXPECT_NEAR(res.cos_sum, expected_cos, 1e-5);
  EXPECT_NEAR(res.sin_sum, expected_sin, 1e-5);
}

TEST_F(SIMDUtilsTests, FillPositionBlockParamsStructWorks) {
  const std::vector<DummyPos> dummy_atoms = {{
                                                 .x_val = 1.0,
                                                 .y_val = 2.0,
                                                 .z_val = 3.0,
                                             },
                                             {
                                                 .x_val = 4.0,
                                                 .y_val = 5.0,
                                                 .z_val = 6.0,
                                             },
                                             {
                                                 .x_val = 7.0,
                                                 .y_val = 8.0,
                                                 .z_val = 9.0,
                                             }};
  std::vector<double> x_s;
  std::vector<double> y_s;
  std::vector<double> z_s;

  const std::size_t count =
      correlation::math::fill_position_block(correlation::math::FillPositionBlockParams<std::vector<DummyPos>, double>{
          .atoms = &dummy_atoms,
          .begin_idx = 0,
          .end_idx = dummy_atoms.size(),
          .x_s = &x_s,
          .y_s = &y_s,
          .z_s = &z_s,
      });

  EXPECT_EQ(count, 3);
  EXPECT_EQ(x_s.size(), 3);
  EXPECT_DOUBLE_EQ(x_s[0], 1.0);
  EXPECT_DOUBLE_EQ(y_s[1], 5.0);
  EXPECT_DOUBLE_EQ(z_s[2], 9.0);
}

TEST_F(SIMDUtilsTests, ComplexExpSumParamsStructMatchesScalar) {
  const std::vector<float> x_s = {1.0F, 2.0F, 3.0F};
  const std::vector<float> y_s = {0.5F, 1.5F, 2.5F};
  const std::vector<float> z_s = {0.1F, 0.2F, 0.3F};
  const float q_x = 0.5F;
  const float q_y = 1.0F;
  const float q_z = 1.5F;

  correlation::math::ComplexExpSumResult<float> result;
  correlation::math::complex_exp_sum(
      correlation::math::ComplexExpSumParams<float>{
          .q_x = q_x,
          .q_y = q_y,
          .q_z = q_z,
          .x_s = x_s.data(),
          .y_s = y_s.data(),
          .z_s = z_s.data(),
          .count = x_s.size(),
      },
      result);

  float expected_cos = 0.0F;
  float expected_sin = 0.0F;
  for (std::size_t i = 0; i < x_s.size(); ++i) {
    const float phase = q_x * x_s[i] + q_y * y_s[i] + q_z * z_s[i];
    expected_cos += std::cos(phase);
    expected_sin += std::sin(phase);
  }

  EXPECT_NEAR(result.cos_sum, expected_cos, 1e-5F);
  EXPECT_NEAR(result.sin_sum, expected_sin, 1e-5F);
}
