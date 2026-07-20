// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only

#include "math/LinearAlgebra.hpp"
#include <cmath>
#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::math;

namespace {
class PrecisionTests : public ::testing::Test {
protected:
  static constexpr double d_x1 = 1.2345678901234567;
  static constexpr double d_y1 = -2.3456789012345678;
  static constexpr double d_z1 = 3.4567890123456789;
  static constexpr double d_x2 = 9.8765432109876543;
  static constexpr double d_y2 = 8.7654321098765432;
  static constexpr double d_z2 = -7.6543210987654321;

  const Vector3<double> vd1{d_x1, d_y1, d_z1};
  const Vector3<double> vd2{d_x2, d_y2, d_z2};

  const Vector3<float> vf1{static_cast<float>(d_x1), static_cast<float>(d_y1), static_cast<float>(d_z1)};
  const Vector3<float> vf2{static_cast<float>(d_x2), static_cast<float>(d_y2), static_cast<float>(d_z2)};
};
} // namespace

TEST_F(PrecisionTests, Vector3AdditionAndSubtraction) {
  // Vector Addition
  double const add_expected = d_x1 + d_x2;
  double const add_double = (vd1 + vd2).x();
  float const add_float = (vf1 + vf2).x();
  EXPECT_NEAR(add_double, add_expected, 1e-14);
  EXPECT_NEAR(add_float, static_cast<float>(add_expected), 1e-5F);

  // Vector Subtraction
  double const sub_expected = d_y1 - d_y2;
  double const sub_double = (vd1 - vd2).y();
  float const sub_float = (vf1 - vf2).y();
  EXPECT_NEAR(sub_double, sub_expected, 1e-14);
  EXPECT_NEAR(sub_float, static_cast<float>(sub_expected), 1e-5F);
}

TEST_F(PrecisionTests, Vector3DotProduct) {
  double const dot_expected = d_x1 * d_x2 + d_y1 * d_y2 + d_z1 * d_z2;
  double const dot_double = dot(vd1, vd2);
  float const dot_float = dot(vf1, vf2);
  EXPECT_NEAR(dot_double, dot_expected, 1e-12);
  EXPECT_NEAR(dot_float, static_cast<float>(dot_expected), 1e-4F);
}

TEST_F(PrecisionTests, Vector3CrossProduct) {
  double const cross_x_expected = d_y1 * d_z2 - d_z1 * d_y2;
  double const cross_x_double = cross(vd1, vd2).x();
  float const cross_x_float = cross(vf1, vf2).x();
  EXPECT_NEAR(cross_x_double, cross_x_expected, 1e-12);
  EXPECT_NEAR(cross_x_float, static_cast<float>(cross_x_expected), 1e-4F);
}

TEST_F(PrecisionTests, Vector3NormAndDistance) {
  double const norm_expected = std::sqrt(d_x1 * d_x1 + d_y1 * d_y1 + d_z1 * d_z1);
  double const norm_double = norm(vd1);
  float const norm_float = norm(vf1);
  EXPECT_NEAR(norm_double, norm_expected, 1e-14);
  EXPECT_NEAR(norm_float, static_cast<float>(norm_expected), 1e-5F);

  Vector3<double> const diff_d = vd1 - vd2;
  double const dist_expected = norm(diff_d);
  double const dist_double = distance(vd1, vd2);
  float const dist_float = distance(vf1, vf2);
  EXPECT_NEAR(dist_double, dist_expected, 1e-14);
  EXPECT_NEAR(dist_float, static_cast<float>(dist_expected), 1e-5F);
}

TEST_F(PrecisionTests, Matrix3Determinant) {
  Vector3<double> const m_c0_d(2.12345678901234, 0.54321098765432, -1.23456789012345);
  Vector3<double> const m_c1_d(-0.98765432109876, 4.32109876543210, 2.34567890123456);
  Vector3<double> const m_c2_d(3.45678901234567, -2.10987654321098, 5.67890123456789);
  Matrix3<double> const m_double(m_c0_d, m_c1_d, m_c2_d);

  Vector3<float> const vf_c0(static_cast<float>(m_c0_d.x()), static_cast<float>(m_c0_d.y()),
                             static_cast<float>(m_c0_d.z()));
  Vector3<float> const vf_c1(static_cast<float>(m_c1_d.x()), static_cast<float>(m_c1_d.y()),
                             static_cast<float>(m_c1_d.z()));
  Vector3<float> const vf_c2(static_cast<float>(m_c2_d.x()), static_cast<float>(m_c2_d.y()),
                             static_cast<float>(m_c2_d.z()));
  Matrix3<float> const m_float(vf_c0, vf_c1, vf_c2);

  double const det_expected = determinant(m_double);
  double const det_double = determinant(m_double);
  float const det_float = determinant(m_float);
  EXPECT_NEAR(det_double, det_expected, 1e-8);
  EXPECT_NEAR(det_float, static_cast<float>(det_expected), 1e-4F);
}

TEST_F(PrecisionTests, Matrix3Inversion) {
  Vector3<double> const m_c0_d(2.12345678901234, 0.54321098765432, -1.23456789012345);
  Vector3<double> const m_c1_d(-0.98765432109876, 4.32109876543210, 2.34567890123456);
  Vector3<double> const m_c2_d(3.45678901234567, -2.10987654321098, 5.67890123456789);
  Matrix3<double> const m_double(m_c0_d, m_c1_d, m_c2_d);

  Vector3<float> const vf_c0(static_cast<float>(m_c0_d.x()), static_cast<float>(m_c0_d.y()),
                             static_cast<float>(m_c0_d.z()));
  Vector3<float> const vf_c1(static_cast<float>(m_c1_d.x()), static_cast<float>(m_c1_d.y()),
                             static_cast<float>(m_c1_d.z()));
  Vector3<float> const vf_c2(static_cast<float>(m_c2_d.x()), static_cast<float>(m_c2_d.y()),
                             static_cast<float>(m_c2_d.z()));
  Matrix3<float> const m_float(vf_c0, vf_c1, vf_c2);

  auto const m_double_inv = invert(m_double);
  auto const m_float_inv = invert(m_float);

  auto const m_double_prod = m_double * m_double_inv;
  auto const m_float_prod = m_float * m_float_inv;

  double inv_err_double = 0.0;
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      double expected_cell = (i == j) ? 1.0 : 0.0;
      inv_err_double += std::abs(m_double_prod(i, j) - expected_cell);
    }
  }

  float inv_err_float = 0.0F;
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      float expected_cell = (i == j) ? 1.0F : 0.0F;
      inv_err_float += std::abs(m_float_prod(i, j) - expected_cell);
    }
  }

  EXPECT_LT(inv_err_double, 1e-14);
  EXPECT_LT(inv_err_float, 1e-5F);
}

TEST_F(PrecisionTests, AccumulationPrecisionLoss) {
  double sum_double = 1.0;
  float sum_float = 1.0F;
  constexpr double delta = 1e-7;
  constexpr int n_iter = 1000000;

  for (int i = 0; i < n_iter; ++i) {
    sum_double += delta;
    sum_float += static_cast<float>(delta);
  }

  double const sum_expected = 1.0 + n_iter * delta;
  EXPECT_NEAR(sum_double, sum_expected, 1e-10);
  EXPECT_TRUE(std::abs(sum_float - static_cast<float>(sum_expected)) > 1e-4F);
}

} // namespace correlation::testing
