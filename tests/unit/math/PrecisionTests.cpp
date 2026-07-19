// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only

#include "math/LinearAlgebra.hpp"
#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::math;

namespace {
class PrecisionTests : public ::testing::Test {};
} // namespace

TEST_F(PrecisionTests, Vector3SingleVsDoublePrecision) {
  // We define identical vectors in float and double
  Vector3<double> const vd_1(1.23456789012345, -2.34567890123456, 3.45678901234567);
  Vector3<float> const vf_1(1.23456789012345F, -2.34567890123456F, 3.45678901234567F);

  Vector3<double> const vd_2(9.87654321098765, 8.76543210987654, -7.65432109876543);
  Vector3<float> const vf_2(9.87654321098765F, 8.76543210987654F, -7.65432109876543F);

  // 1. Addition precision
  auto const vd_sum = vd_1 + vd_2;
  auto const vf_sum = vf_1 + vf_2;
  EXPECT_NEAR(vd_sum.x(), vf_sum.x(), 1e-6);
  EXPECT_NEAR(vd_sum.y(), vf_sum.y(), 1e-6);
  EXPECT_NEAR(vd_sum.z(), vf_sum.z(), 1e-6);

  // 2. Dot product precision
  double const vd_dot = dot(vd_1, vd_2);
  float const vf_dot = dot(vf_1, vf_2);
  EXPECT_NEAR(vd_dot, vf_dot, 1e-5);

  // 3. Cross product precision
  auto const vd_cross = cross(vd_1, vd_2);
  auto const vf_cross = cross(vf_1, vf_2);
  EXPECT_NEAR(vd_cross.x(), vf_cross.x(), 1e-5);
  EXPECT_NEAR(vd_cross.y(), vf_cross.y(), 1e-5);
  EXPECT_NEAR(vd_cross.z(), vf_cross.z(), 1e-5);

  // 4. Norm precision
  double const vd_norm = norm(vd_1);
  float const vf_norm = norm(vf_1);
  EXPECT_NEAR(vd_norm, vf_norm, 1e-6);
}

TEST_F(PrecisionTests, Matrix3SingleVsDoublePrecision) {
  Vector3<double> const vd_c0(2.123456789, 0.543210987, -1.23456789);
  Vector3<double> const vd_c1(-0.987654321, 4.321098765, 2.345678901);
  Vector3<double> const vd_c2(3.456789012, -2.109876543, 5.678901234);
  Matrix3<double> const m_double(vd_c0, vd_c1, vd_c2);

  Vector3<float> const vf_c0(2.123456789F, 0.543210987F, -1.23456789F);
  Vector3<float> const vf_c1(-0.987654321F, 4.321098765F, 2.345678901F);
  Vector3<float> const vf_c2(3.456789012F, -2.109876543F, 5.678901234F);
  Matrix3<float> const m_float(vf_c0, vf_c1, vf_c2);

  // 1. Determinant precision
  double const det_d = determinant(m_double);
  float const det_f = determinant(m_float);
  EXPECT_NEAR(det_d, det_f, 1e-5);

  // 2. Inversion precision
  auto const m_double_inv = invert(m_double);
  auto const m_float_inv = invert(m_float);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      EXPECT_NEAR(m_double_inv(i, j), m_float_inv(i, j), 1e-5);
    }
  }
}

} // namespace correlation::testing
