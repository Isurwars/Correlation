// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "math/LinearAlgebra.hpp"

#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::math;
namespace {
class Vector3Tests : public ::testing::Test {};
} // namespace

TEST_F(Vector3Tests, Vector3ConstructorsAndAccessors) {
  // Default constructor
  Vector3<double> v_1;
  EXPECT_DOUBLE_EQ(v_1.x(), 0.0);
  EXPECT_DOUBLE_EQ(v_1.y(), 0.0);
  EXPECT_DOUBLE_EQ(v_1.z(), 0.0);
  EXPECT_TRUE(v_1.empty());

  // Parameterized constructor
  Vector3<double> v_2(1.5, -2.5, 3.0);
  EXPECT_DOUBLE_EQ(v_2.x(), 1.5);
  EXPECT_DOUBLE_EQ(v_2.y(), -2.5);
  EXPECT_DOUBLE_EQ(v_2.z(), 3.0);
  EXPECT_FALSE(v_2.empty());

  // Array constructor
  std::array<double, 3> const arr = {10.0, 20.0, 30.0};
  Vector3<double> v_3(arr);
  EXPECT_DOUBLE_EQ(v_3[0], 10.0);
  EXPECT_DOUBLE_EQ(v_3[1], 20.0);
  EXPECT_DOUBLE_EQ(v_3[2], 30.0);
  EXPECT_EQ(v_3.array(), arr);

  // Mutable access
  v_3.x() = 1.0;
  v_3[1] = 2.0;
  v_3(2) = 3.0;
  EXPECT_DOUBLE_EQ(v_3.x(), 1.0);
  EXPECT_DOUBLE_EQ(v_3.y(), 2.0);
  EXPECT_DOUBLE_EQ(v_3.z(), 3.0);
}

TEST_F(Vector3Tests, Vector3Arithmetic) {
  Vector3<double> const vec_a(1.0, 2.0, 3.0);
  Vector3<double> const vec_b(4.0, 5.0, 6.0);

  // Addition
  auto vec_c = vec_a + vec_b;
  EXPECT_DOUBLE_EQ(vec_c.x(), 5.0);
  EXPECT_DOUBLE_EQ(vec_c.y(), 7.0);
  EXPECT_DOUBLE_EQ(vec_c.z(), 9.0);

  // Subtraction
  auto vec_d = vec_a - vec_b;
  EXPECT_DOUBLE_EQ(vec_d.x(), -3.0);
  EXPECT_DOUBLE_EQ(vec_d.y(), -3.0);
  EXPECT_DOUBLE_EQ(vec_d.z(), -3.0);

  // Scalar multiplication / division
  auto vec_e = vec_a * 2.0;
  EXPECT_DOUBLE_EQ(vec_e.x(), 2.0);
  EXPECT_DOUBLE_EQ(vec_e.y(), 4.0);
  EXPECT_DOUBLE_EQ(vec_e.z(), 6.0);

  auto vec_f = vec_b / 2.0;
  EXPECT_DOUBLE_EQ(vec_f.x(), 2.0);
  EXPECT_DOUBLE_EQ(vec_f.y(), 2.5);
  EXPECT_DOUBLE_EQ(vec_f.z(), 3.0);

  // In-place operators
  Vector3<double> vec_g = vec_a;
  vec_g += vec_b;
  EXPECT_DOUBLE_EQ(vec_g.x(), 5.0);

  vec_g -= vec_b;
  EXPECT_DOUBLE_EQ(vec_g.x(), 1.0);

  // Dot product operator
  double const dot_ab = vec_a * vec_b;
  EXPECT_DOUBLE_EQ(dot_ab, 32.0);
}

TEST_F(Vector3Tests, Vector3FreeFunctions) {
  Vector3<double> const vec_a(1.0, 2.0, 3.0);
  Vector3<double> const vec_b(4.0, -5.0, 6.0);

  // dot
  EXPECT_DOUBLE_EQ(dot(vec_a, vec_b), 12.0);

  // cross
  auto vec_cross = cross(vec_a, vec_b);
  EXPECT_DOUBLE_EQ(vec_cross.x(), 2.0 * 6.0 - 3.0 * (-5.0)); // 27
  EXPECT_DOUBLE_EQ(vec_cross.y(), 3.0 * 4.0 - 1.0 * 6.0);    // 6
  EXPECT_DOUBLE_EQ(vec_cross.z(), 1.0 * (-5.0) - 2.0 * 4.0); // -13

  // norm_sq and norm
  Vector3<double> const vec_v(3.0, 4.0, 0.0);
  EXPECT_DOUBLE_EQ(norm_sq(vec_v), 25.0);
  EXPECT_DOUBLE_EQ(norm(vec_v), 5.0);

  // normalize
  auto vec_vn = normalize(vec_v);
  EXPECT_DOUBLE_EQ(vec_vn.x(), 0.6);
  EXPECT_DOUBLE_EQ(vec_vn.y(), 0.8);
  EXPECT_DOUBLE_EQ(vec_vn.z(), 0.0);

  // normalize singular vector throws
  Vector3<double> const zero_v;
  EXPECT_THROW((void)normalize(zero_v), std::domain_error);

  Vector3<double> const tiny_v(1e-301, 1e-301, 1e-301);
  EXPECT_THROW((void)normalize(tiny_v), std::domain_error);
}

TEST_F(Vector3Tests, Vector3DivisionByZero) {
  Vector3<double> const vec_v(1.0, 2.0, 3.0);
  auto result = vec_v / 0.0;

  // Division by zero should produce infinity
  EXPECT_TRUE(std::isinf(result.x()));
  EXPECT_TRUE(std::isinf(result.y()));
  EXPECT_TRUE(std::isinf(result.z()));

  // Zero / zero should produce NaN
  Vector3<double> const zero_v;
  auto nan_result = zero_v / 0.0;
  EXPECT_TRUE(std::isnan(nan_result.x()));
}

TEST_F(Vector3Tests, Vector3UnaryNegation) {
  Vector3<double> const vec_v(1.0, -2.5, 3.0);
  auto vec_neg = -vec_v;
  EXPECT_DOUBLE_EQ(vec_neg.x(), -1.0);
  EXPECT_DOUBLE_EQ(vec_neg.y(), 2.5);
  EXPECT_DOUBLE_EQ(vec_neg.z(), -3.0);

  // Double negation should return original
  auto vec_double_neg = -(-vec_v);
  EXPECT_DOUBLE_EQ(vec_double_neg.x(), vec_v.x());
  EXPECT_DOUBLE_EQ(vec_double_neg.y(), vec_v.y());
  EXPECT_DOUBLE_EQ(vec_double_neg.z(), vec_v.z());

  // Negation of zero vector
  Vector3<double> const vec_zero;
  auto vec_neg_zero = -vec_zero;
  EXPECT_TRUE(vec_neg_zero.empty());
}

TEST_F(Vector3Tests, Vector3InPlaceScalarOps) {
  Vector3<double> vec_v1(2.0, 4.0, 6.0);
  vec_v1 *= 3.0;
  EXPECT_DOUBLE_EQ(vec_v1.x(), 6.0);
  EXPECT_DOUBLE_EQ(vec_v1.y(), 12.0);
  EXPECT_DOUBLE_EQ(vec_v1.z(), 18.0);

  vec_v1 /= 2.0;
  EXPECT_DOUBLE_EQ(vec_v1.x(), 3.0);
  EXPECT_DOUBLE_EQ(vec_v1.y(), 6.0);
  EXPECT_DOUBLE_EQ(vec_v1.z(), 9.0);

  // Multiply by zero
  vec_v1 *= 0.0;
  EXPECT_TRUE(vec_v1.empty());
}

TEST_F(Vector3Tests, DistanceFunction) {
  Vector3<double> const vec_v1(1.0, 0.0, 0.0);
  Vector3<double> const vec_v2(4.0, 0.0, 0.0);
  EXPECT_DOUBLE_EQ(distance(vec_v1, vec_v2), 3.0);

  // Distance is symmetric
  EXPECT_DOUBLE_EQ(distance(vec_v1, vec_v2), distance(vec_v2, vec_v1));

  // Distance to self is zero
  EXPECT_DOUBLE_EQ(distance(vec_v1, vec_v1), 0.0);

  // 3-4-5 triangle
  Vector3<double> const origin(0.0, 0.0, 0.0);
  Vector3<double> const point(3.0, 4.0, 0.0);
  EXPECT_DOUBLE_EQ(distance(origin, point), 5.0);
}

TEST_F(Vector3Tests, Vector3Equality) {
  Vector3<double> const vec_v1(1.0, 2.0, 3.0);
  Vector3<double> const vec_v2(1.0, 2.0, 3.0);
  Vector3<double> const vec_v3(1.0, 2.0, 3.1);

  EXPECT_TRUE(vec_v1 == vec_v2);
  EXPECT_FALSE(vec_v1 != vec_v2);

  EXPECT_FALSE(vec_v1 == vec_v3);
  EXPECT_TRUE(vec_v1 != vec_v3);

  // Zero vectors
  Vector3<double> const vec_z1;
  Vector3<double> const vec_z2;
  EXPECT_TRUE(vec_z1 == vec_z2);
}

TEST_F(Vector3Tests, Vector3BeginEnd) {
  Vector3<double> const vec_v(10.0, 20.0, 30.0);

  // Verify pointer-based iteration
  const double *iter = vec_v.begin();
  EXPECT_DOUBLE_EQ(*iter, 10.0);
  ++iter;
  EXPECT_DOUBLE_EQ(*iter, 20.0);
  ++iter;
  EXPECT_DOUBLE_EQ(*iter, 30.0);
  ++iter;
  EXPECT_EQ(iter, vec_v.end());

  // Range-for loop reads
  double sum = 0.0;
  for (double val : vec_v) {
    sum += val;
  }
  EXPECT_DOUBLE_EQ(sum, 60.0);

  // Mutable iteration
  Vector3<double> vec_v2(1.0, 2.0, 3.0);
  for (double &val : vec_v2) {
    val *= 10.0;
  }
  EXPECT_DOUBLE_EQ(vec_v2.x(), 10.0);
  EXPECT_DOUBLE_EQ(vec_v2.y(), 20.0);
  EXPECT_DOUBLE_EQ(vec_v2.z(), 30.0);

  // std::distance between begin and end
  EXPECT_EQ(vec_v.end() - vec_v.begin(), 3);
}

TEST_F(Vector3Tests, ScalarTimesVector) {
  Vector3<double> const vec_v(2.0, 3.0, 4.0);

  // scalar * vector should equal vector * scalar
  auto lhs = 5.0 * vec_v;
  auto rhs = vec_v * 5.0;
  EXPECT_DOUBLE_EQ(lhs.x(), rhs.x());
  EXPECT_DOUBLE_EQ(lhs.y(), rhs.y());
  EXPECT_DOUBLE_EQ(lhs.z(), rhs.z());

  // Works with int scalar
  auto int_result = 3 * vec_v;
  EXPECT_DOUBLE_EQ(int_result.x(), 6.0);
  EXPECT_DOUBLE_EQ(int_result.y(), 9.0);
  EXPECT_DOUBLE_EQ(int_result.z(), 12.0);

  // Works with float scalar
  auto float_result = 2.0F * Vector3<float>(1.0F, 2.0F, 3.0F);
  EXPECT_FLOAT_EQ(float_result.x(), 2.0F);
  EXPECT_FLOAT_EQ(float_result.y(), 4.0F);
  EXPECT_FLOAT_EQ(float_result.z(), 6.0F);
}

TEST_F(Vector3Tests, CrossProductProperties) {
  Vector3<double> const vec_1(1.0, 2.0, 3.0);
  Vector3<double> const vec_2(4.0, -5.0, 6.0);

  // Anti-commutativity: a x b = -(b x a)
  auto axb = cross(vec_1, vec_2);
  auto bxa = cross(vec_2, vec_1);
  EXPECT_DOUBLE_EQ(axb.x(), -bxa.x());
  EXPECT_DOUBLE_EQ(axb.y(), -bxa.y());
  EXPECT_DOUBLE_EQ(axb.z(), -bxa.z());

  // Orthogonality: (a x b) . a = 0 and (a x b) . b = 0
  EXPECT_NEAR(dot(axb, vec_1), 0.0, 1e-15);
  EXPECT_NEAR(dot(axb, vec_2), 0.0, 1e-15);

  // Self cross product is zero
  auto axa = cross(vec_1, vec_1);
  EXPECT_NEAR(axa.x(), 0.0, 1e-15);
  EXPECT_NEAR(axa.y(), 0.0, 1e-15);
  EXPECT_NEAR(axa.z(), 0.0, 1e-15);

  // Basis vectors: x × y = z
  Vector3<double> const vec_ex(1.0, 0.0, 0.0);
  Vector3<double> const vec_ey(0.0, 1.0, 0.0);
  Vector3<double> const vec_ez(0.0, 0.0, 1.0);
  auto vec_x_cross_y = cross(vec_ex, vec_ey);
  EXPECT_DOUBLE_EQ(vec_x_cross_y.x(), 0.0);
  EXPECT_DOUBLE_EQ(vec_x_cross_y.y(), 0.0);
  EXPECT_DOUBLE_EQ(vec_x_cross_y.z(), 1.0);

  // Lagrange's identity: |a x b|^2 = |a|^2*|b|^2 - (a.b)^2
  double const lhs = norm_sq(axb);
  double const rhs = norm_sq(vec_1) * norm_sq(vec_2) - dot(vec_1, vec_2) * dot(vec_1, vec_2);
  EXPECT_NEAR(lhs, rhs, 1e-10);
}

TEST_F(Vector3Tests, Vector3IntType) {
  Vector3<int> const vec_a(1, 2, 3);
  Vector3<int> const vec_b(4, 5, 6);

  auto sum = vec_a + vec_b;
  EXPECT_EQ(sum.x(), 5);
  EXPECT_EQ(sum.y(), 7);
  EXPECT_EQ(sum.z(), 9);

  auto diff = vec_a - vec_b;
  EXPECT_EQ(diff.x(), -3);
  EXPECT_EQ(diff.y(), -3);
  EXPECT_EQ(diff.z(), -3);

  auto neg = -vec_a;
  EXPECT_EQ(neg.x(), -1);
  EXPECT_EQ(neg.y(), -2);
  EXPECT_EQ(neg.z(), -3);

  auto scaled = vec_a * 3;
  EXPECT_EQ(scaled.x(), 3);
  EXPECT_EQ(scaled.y(), 6);
  EXPECT_EQ(scaled.z(), 9);

  EXPECT_EQ(dot(vec_a, vec_b), 32);

  auto vec_cross = cross(vec_a, vec_b);
  EXPECT_EQ(vec_cross.x(), 2 * 6 - 3 * 5); // -3
  EXPECT_EQ(vec_cross.y(), 3 * 4 - 1 * 6); // 6
  EXPECT_EQ(vec_cross.z(), 1 * 5 - 2 * 4); // -3

  Vector3<int> const vec_a2(1, 2, 3);
  EXPECT_TRUE(vec_a == vec_a2);
  EXPECT_FALSE(vec_a != vec_a2);
  EXPECT_FALSE(vec_a == vec_b);

  Vector3<int> const zero;
  EXPECT_TRUE(zero.empty());
  EXPECT_FALSE(vec_a.empty());
}

TEST_F(Vector3Tests, Vector3FloatType) {
  Vector3<float> const vec_a(1.0F, 2.0F, 3.0F);
  Vector3<float> const vec_b(4.0F, 5.0F, 6.0F);

  auto sum = vec_a + vec_b;
  EXPECT_FLOAT_EQ(sum.x(), 5.0F);
  EXPECT_FLOAT_EQ(sum.y(), 7.0F);
  EXPECT_FLOAT_EQ(sum.z(), 9.0F);

  EXPECT_FLOAT_EQ(dot(vec_a, vec_b), 32.0F);
  EXPECT_FLOAT_EQ(norm_sq(vec_a), 14.0F);

  auto normalized = normalize(vec_a);
  EXPECT_NEAR(norm(normalized), 1.0F, 1e-6F);
}

TEST_F(Vector3Tests, NormProperties) {
  Vector3<double> const vec_v(3.0, 4.0, 0.0);

  EXPECT_DOUBLE_EQ(norm_sq(vec_v), norm(vec_v) * norm(vec_v));

  Vector3<double> const zero;
  EXPECT_DOUBLE_EQ(norm(zero), 0.0);
  EXPECT_DOUBLE_EQ(norm_sq(zero), 0.0);

  Vector3<double> const neg(-3.0, -4.0, -5.0);
  EXPECT_GE(norm(neg), 0.0);

  Vector3<double> const vec_a(1.0, 2.0, 3.0);
  Vector3<double> const vec_b(4.0, -1.0, 2.0);
  EXPECT_LE(norm(vec_a + vec_b), norm(vec_a) + norm(vec_b) + 1e-15);

  double const scalar = -3.5;
  EXPECT_NEAR(norm(vec_v * scalar), std::abs(scalar) * norm(vec_v), 1e-15);

  auto normalized = normalize(vec_v);
  EXPECT_NEAR(norm(normalized), 1.0, 1e-15);
}

} // namespace correlation::testing
