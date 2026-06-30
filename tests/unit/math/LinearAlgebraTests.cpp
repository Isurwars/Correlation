// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "math/LinearAlgebra.hpp"

#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::math;
namespace {
class LinearAlgebraTests : public ::testing::Test {};
} // namespace
// --- Vector3 Tests ---

TEST_F(LinearAlgebraTests, Vector3ConstructorsAndAccessors) {
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

TEST_F(LinearAlgebraTests, Vector3Arithmetic) {
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

TEST_F(LinearAlgebraTests, Vector3FreeFunctions) {
  Vector3<double> const vec_a(1.0, 2.0, 3.0);
  Vector3<double> const vec_b(4.0, -5.0, 6.0);

  // dot
  EXPECT_DOUBLE_EQ(dot(vec_a, vec_b), 12.0);

  // cross
  auto vec_cross = cross(vec_a, vec_b);
  EXPECT_DOUBLE_EQ(vec_cross.x(), 2.0 * 6.0 - 3.0 * (-5.0)); // 27
  EXPECT_DOUBLE_EQ(vec_cross.y(), 3.0 * 4.0 - 1.0 * 6.0);    // 6
  EXPECT_DOUBLE_EQ(vec_cross.z(), 1.0 * (-5.0) - 2.0 * 4.0); // -1 3

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

// --- Matrix3 Tests ---

TEST_F(LinearAlgebraTests, Matrix3ConstructorsAndAccessors) {
  // Default constructor (zeros)
  Matrix3<double> mat_1;
  EXPECT_DOUBLE_EQ(mat_1(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(mat_1.trace(), 0.0);

  // Column constructor
  Vector3<double> const vec_c0(1.0, 2.0, 3.0);
  Vector3<double> const vec_c1(4.0, 5.0, 6.0);
  Vector3<double> const vec_c2(7.0, 8.0, 9.0);
  Matrix3<double> mat_2(vec_c0, vec_c1, vec_c2);

  EXPECT_DOUBLE_EQ(mat_2(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(mat_2(1, 0), 2.0); // Row 1, Col 0
  EXPECT_DOUBLE_EQ(mat_2(0, 1), 4.0); // Row 0, Col 1
  EXPECT_DOUBLE_EQ(mat_2(2, 2), 9.0);
  EXPECT_DOUBLE_EQ(mat_2.trace(), 15.0);

  // Mutable access
  mat_2[0][0] = -1.0;
  EXPECT_DOUBLE_EQ(mat_2(0, 0), -1.0);
}

TEST_F(LinearAlgebraTests, Matrix3Operations) {
  Vector3<double> const vec_c0(1.0, 0.0, 0.0);
  Vector3<double> const vec_c1(0.0, 2.0, 0.0);
  Vector3<double> const vec_c2(0.0, 0.0, 3.0);
  Matrix3<double> const mat_m(vec_c0, vec_c1, vec_c2);

  // Matrix * Vector
  Vector3<double> const vec_v(1.0, 2.0, 3.0);
  auto column_vec = mat_m * vec_v;
  EXPECT_DOUBLE_EQ(column_vec.x(), 1.0);
  EXPECT_DOUBLE_EQ(column_vec.y(), 4.0);
  EXPECT_DOUBLE_EQ(column_vec.z(), 9.0);

  // Matrix * Matrix
  auto mat_prod = mat_m * mat_m;
  EXPECT_DOUBLE_EQ(mat_prod(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(mat_prod(1, 1), 4.0);
  EXPECT_DOUBLE_EQ(mat_prod(2, 2), 9.0);

  // Scalar multiplication
  auto m_scaled = mat_m * 2.0;
  EXPECT_DOUBLE_EQ(m_scaled(1, 1), 4.0);
}

TEST_F(LinearAlgebraTests, DeterminantAndInversion) {
  // Simple orthogonal matrix
  Vector3<double> const vec_c0(2.0, 0.0, 0.0);
  Vector3<double> const vec_c1(0.0, 4.0, 0.0);
  Vector3<double> const vec_c2(0.0, 0.0, 5.0);
  Matrix3<double> const mat_m(vec_c0, vec_c1, vec_c2);

  EXPECT_DOUBLE_EQ(determinant(mat_m), 40.0);

  auto minv = invert(mat_m);
  EXPECT_DOUBLE_EQ(minv(0, 0), 0.5);
  EXPECT_DOUBLE_EQ(minv(1, 1), 0.25);
  EXPECT_DOUBLE_EQ(minv(2, 2), 0.2);

  // singular matrix determinant = 0
  Vector3<double> const c_deg0(1.0, 2.0, 3.0);
  Vector3<double> const c_deg1(2.0, 4.0, 6.0); // Linearly dependent
  Vector3<double> const c_deg2(0.0, 0.0, 1.0);
  Matrix3<double> const m_deg(c_deg0, c_deg1, c_deg2);

  EXPECT_NEAR(determinant(m_deg), 0.0, 1e-15);
  EXPECT_THROW((void)invert(m_deg), std::runtime_error);
}

TEST_F(LinearAlgebraTests, Transpose) {
  Vector3<double> const vec_col_0(1.0, 2.0, 3.0);
  Vector3<double> const vec_col_1(4.0, 5.0, 6.0);
  Vector3<double> const vec_col_2(7.0, 8.0, 9.0);
  Matrix3<double> const mat_m(vec_col_0, vec_col_1, vec_col_2);

  auto mat_t = transpose(mat_m);
  EXPECT_DOUBLE_EQ(mat_t(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(mat_t(1, 0), 4.0);
}

// --- Extreme / Edge-Case Tests ---

TEST_F(LinearAlgebraTests, Vector3DivisionByZero) {
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

TEST_F(LinearAlgebraTests, Matrix3IdentityInverse) {
  // Identity matrix inverse should be identity
  Vector3<double> const vec_col_0(1.0, 0.0, 0.0);
  Vector3<double> const vec_col_1(0.0, 1.0, 0.0);
  Vector3<double> const vec_col_2(0.0, 0.0, 1.0);
  Matrix3<double> const mat_identity(vec_col_0, vec_col_1, vec_col_2);

  EXPECT_DOUBLE_EQ(determinant(mat_identity), 1.0);

  auto inv = invert(mat_identity);
  EXPECT_DOUBLE_EQ(inv(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(inv(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(inv(2, 2), 1.0);
  EXPECT_DOUBLE_EQ(inv(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(inv(0, 2), 0.0);
  EXPECT_DOUBLE_EQ(inv(1, 0), 0.0);
}

TEST_F(LinearAlgebraTests, Matrix3NegativeDeterminant) {
  // Left-handed coordinate system: negative determinant
  Vector3<double> const vec_col_0(0.0, 1.0, 0.0);
  Vector3<double> const vec_col_1(1.0, 0.0, 0.0);
  Vector3<double> const vec_col_2(0.0, 0.0, 1.0);
  Matrix3<double> const mat_m(vec_col_0, vec_col_1, vec_col_2);

  EXPECT_DOUBLE_EQ(determinant(mat_m), -1.0);

  // Should still be invertible
  auto inv = invert(mat_m);
  // M * M^-1 should give identity
  auto product = mat_m * inv;
  EXPECT_NEAR(product(0, 0), 1.0, 1e-15);
  EXPECT_NEAR(product(1, 1), 1.0, 1e-15);
  EXPECT_NEAR(product(2, 2), 1.0, 1e-15);
  EXPECT_NEAR(product(0, 1), 0.0, 1e-15);
}

TEST_F(LinearAlgebraTests, Matrix3ExtremeValues) {
  // Very large values
  double const big = 1e10;
  Vector3<double> const vec_col_0(big, 0.0, 0.0);
  Vector3<double> const vec_col_1(0.0, big, 0.0);
  Vector3<double> const vec_col_2(0.0, 0.0, big);
  Matrix3<double> const mat_big(vec_col_0, vec_col_1, vec_col_2);

  EXPECT_NEAR(determinant(mat_big), big * big * big, big * big * 1e-6);

  auto inv_big = invert(mat_big);
  EXPECT_NEAR(inv_big(0, 0), 1.0 / big, 1e-25);

  // Very small values
  double const small = 1e-4;
  Vector3<double> const vec_small_0(small, 0.0, 0.0);
  Vector3<double> const vec_small_1(0.0, small, 0.0);
  Vector3<double> const vec_small_2(0.0, 0.0, small);
  Matrix3<double> const mat_small(vec_small_0, vec_small_1, vec_small_2);

  EXPECT_NEAR(determinant(mat_small), small * small * small, 1e-20);

  auto inv_small = invert(mat_small);
  EXPECT_NEAR(inv_small(0, 0), 1.0 / small, 1e-6);
}

// --- New Operation Tests ---

TEST_F(LinearAlgebraTests, Vector3UnaryNegation) {
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

TEST_F(LinearAlgebraTests, Vector3InPlaceScalarOps) {
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

TEST_F(LinearAlgebraTests, DistanceFunction) {
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

TEST_F(LinearAlgebraTests, Matrix3AdditionSubtraction) {
  Vector3<double> const vec_col_0(1.0, 2.0, 3.0);
  Vector3<double> const vec_col_1(4.0, 5.0, 6.0);
  Vector3<double> const vec_col_2(7.0, 8.0, 9.0);
  Matrix3<double> const mat_m1(vec_col_0, vec_col_1, vec_col_2);

  Vector3<double> const vec_col_3(9.0, 8.0, 7.0);
  Vector3<double> const vec_col_4(6.0, 5.0, 4.0);
  Vector3<double> const vec_col_5(3.0, 2.0, 1.0);
  Matrix3<double> const mat_m2(vec_col_3, vec_col_4, vec_col_5);

  // Addition
  auto mat_sum = mat_m1 + mat_m2;
  EXPECT_DOUBLE_EQ(mat_sum(0, 0), 10.0);
  EXPECT_DOUBLE_EQ(mat_sum(1, 1), 10.0);
  EXPECT_DOUBLE_EQ(mat_sum(2, 2), 10.0);

  // Subtraction
  auto diff = mat_m1 - mat_m2;
  EXPECT_DOUBLE_EQ(diff(0, 0), -8.0);
  EXPECT_DOUBLE_EQ(diff(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(diff(2, 2), 8.0);

  // In-place subtraction
  Matrix3<double> mat_m3 = mat_m1;
  mat_m3 -= mat_m2;
  EXPECT_DOUBLE_EQ(mat_m3(0, 0), diff(0, 0));
  EXPECT_DOUBLE_EQ(mat_m3(1, 1), diff(1, 1));
  EXPECT_DOUBLE_EQ(mat_m3(2, 2), diff(2, 2));

  // M - M should be zero
  Matrix3<double> const mat_m1_copy = mat_m1;
  auto mat_zero = mat_m1 - mat_m1_copy;
  EXPECT_DOUBLE_EQ(mat_zero(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(mat_zero(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(mat_zero(2, 2), 0.0);
}

TEST_F(LinearAlgebraTests, Matrix3Identity) {
  auto mat_ident = Matrix3<double>::identity();
  EXPECT_DOUBLE_EQ(mat_ident(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(mat_ident(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(mat_ident(2, 2), 1.0);
  EXPECT_DOUBLE_EQ(mat_ident(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(mat_ident(0, 2), 0.0);
  EXPECT_DOUBLE_EQ(mat_ident(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(mat_ident(1, 2), 0.0);
  EXPECT_DOUBLE_EQ(mat_ident(2, 0), 0.0);
  EXPECT_DOUBLE_EQ(mat_ident(2, 1), 0.0);
  EXPECT_DOUBLE_EQ(mat_ident.trace(), 3.0);
  EXPECT_DOUBLE_EQ(determinant(mat_ident), 1.0);

  // Identity * vector = vector
  Vector3<double> const vec_v(3.0, 7.0, -2.0);
  auto mat_result = mat_ident * vec_v;
  EXPECT_DOUBLE_EQ(mat_result.x(), vec_v.x());
  EXPECT_DOUBLE_EQ(mat_result.y(), vec_v.y());
  EXPECT_DOUBLE_EQ(mat_result.z(), vec_v.z());

  // Identity * matrix = matrix
  Vector3<double> const vec_col_0(1.0, 2.0, 3.0);
  Vector3<double> const vec_col_1(4.0, 5.0, 6.0);
  Vector3<double> const vec_col_2(7.0, 8.0, 9.0);
  Matrix3<double> const mat_m(vec_col_0, vec_col_1, vec_col_2);
  auto mat_product = mat_ident * mat_m;
  EXPECT_DOUBLE_EQ(mat_product(0, 0), mat_m(0, 0));
  EXPECT_DOUBLE_EQ(mat_product(1, 1), mat_m(1, 1));
  EXPECT_DOUBLE_EQ(mat_product(2, 2), mat_m(2, 2));
}

TEST_F(LinearAlgebraTests, Vector3Equality) {
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

TEST_F(LinearAlgebraTests, Matrix3Equality) {
  Vector3<double> const vec_col_0(1.0, 2.0, 3.0);
  Vector3<double> const vec_col_1(4.0, 5.0, 6.0);
  Vector3<double> const vec_col_2(7.0, 8.0, 9.0);
  Matrix3<double> const mat_m1(vec_col_0, vec_col_1, vec_col_2);
  Matrix3<double> const mat_m2(vec_col_0, vec_col_1, vec_col_2);

  EXPECT_TRUE(mat_m1 == mat_m2);
  EXPECT_FALSE(mat_m1 != mat_m2);

  Matrix3<double> mat_m3 = mat_m1;
  mat_m3(0, 0) = 99.0;
  EXPECT_FALSE(mat_m1 == mat_m3);
  EXPECT_TRUE(mat_m1 != mat_m3);

  // Identity comparisons
  Matrix3<double> const mat_ident = Matrix3<double>::identity();
  Matrix3<double> const mat_zero;
  EXPECT_FALSE(mat_ident == mat_zero);
  EXPECT_TRUE(mat_ident != mat_zero);
}

// --- Iterator Tests ---

TEST_F(LinearAlgebraTests, Vector3BeginEnd) {
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

// --- Commutative Scalar-Vector Multiplication ---

TEST_F(LinearAlgebraTests, ScalarTimesVector) {
  Vector3<double> const vec_v(2.0, 3.0, 4.0);

  // scalar * vector should equal vector * scalar
  auto lhs = 5.0 * vec_v;
  auto rhs = vec_v * 5.0;
  EXPECT_DOUBLE_EQ(lhs.x(), rhs.x());
  EXPECT_DOUBLE_EQ(lhs.y(), rhs.y());
  EXPECT_DOUBLE_EQ(lhs.z(), rhs.z());

  // Works with int scalar (cross-type via requires clause)
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

// --- Matrix3 In-Place Operations ---

TEST_F(LinearAlgebraTests, Matrix3InPlaceAddition) {
  Vector3<double> const vec_col_0(1.0, 2.0, 3.0);
  Vector3<double> const vec_col_1(4.0, 5.0, 6.0);
  Vector3<double> const vec_col_2(7.0, 8.0, 9.0);
  Matrix3<double> const mat_m1(vec_col_0, vec_col_1, vec_col_2);

  Vector3<double> const vec_col_3(10.0, 20.0, 30.0);
  Vector3<double> const vec_col_4(40.0, 50.0, 60.0);
  Vector3<double> const vec_col_5(70.0, 80.0, 90.0);
  Matrix3<double> const mat_m2(vec_col_3, vec_col_4, vec_col_5);

  Matrix3<double> mat_m3 = mat_m1;
  mat_m3 += mat_m2;
  EXPECT_DOUBLE_EQ(mat_m3(0, 0), 11.0);
  EXPECT_DOUBLE_EQ(mat_m3(1, 1), 55.0);
  EXPECT_DOUBLE_EQ(mat_m3(2, 2), 99.0);
  EXPECT_DOUBLE_EQ(mat_m3(0, 1), 44.0);
  EXPECT_DOUBLE_EQ(mat_m3(2, 0), 33.0);
}

TEST_F(LinearAlgebraTests, Matrix3InPlaceScalarMultiply) {
  Vector3<double> const vec_col_0(1.0, 0.0, 0.0);
  Vector3<double> const vec_col_1(0.0, 2.0, 0.0);
  Vector3<double> const vec_col_2(0.0, 0.0, 3.0);
  Matrix3<double> mat_m1(vec_col_0, vec_col_1, vec_col_2);

  Matrix3<double> mat_m2 = mat_m1;
  mat_m2 *= 4.0;
  EXPECT_DOUBLE_EQ(mat_m2(0, 0), 4.0);
  EXPECT_DOUBLE_EQ(mat_m2(1, 1), 8.0);
  EXPECT_DOUBLE_EQ(mat_m2(2, 2), 12.0);
  EXPECT_DOUBLE_EQ(mat_m2(0, 1), 0.0); // Off-diagonal stays zero

  // Multiply by zero
  mat_m2 *= 0.0;
  EXPECT_DOUBLE_EQ(mat_m2(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(mat_m2(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(mat_m2(2, 2), 0.0);
}

TEST_F(LinearAlgebraTests, Matrix3InPlaceMatrixMultiply) {
  // Diagonal matrices: multiplication is element-wise on diagonal
  Vector3<double> const vec_col_0(2.0, 0.0, 0.0);
  Vector3<double> const vec_col_1(0.0, 3.0, 0.0);
  Vector3<double> const vec_col_2(0.0, 0.0, 5.0);
  Matrix3<double> mat_m1(vec_col_0, vec_col_1, vec_col_2);

  Vector3<double> const vec_col_3(4.0, 0.0, 0.0);
  Vector3<double> const vec_col_4(0.0, 2.0, 0.0);
  Vector3<double> const vec_col_5(0.0, 0.0, 3.0);
  Matrix3<double> const mat_m2(vec_col_3, vec_col_4, vec_col_5);

  mat_m1 *= mat_m2;
  EXPECT_DOUBLE_EQ(mat_m1(0, 0), 8.0);
  EXPECT_DOUBLE_EQ(mat_m1(1, 1), 6.0);
  EXPECT_DOUBLE_EQ(mat_m1(2, 2), 15.0);
  EXPECT_DOUBLE_EQ(mat_m1(0, 1), 0.0);

  // Multiply by identity leaves matrix unchanged
  Vector3<double> const vec_col_6(1.0, 2.0, 3.0);
  Vector3<double> const vec_col_7(4.0, 5.0, 6.0);
  Vector3<double> const vec_col_8(7.0, 8.0, 9.0);
  Matrix3<double> const mat_m3(vec_col_6, vec_col_7, vec_col_8);
  Matrix3<double> mat_m3_copy = mat_m3;
  mat_m3_copy *= Matrix3<double>::identity();
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t col = 0; col < 3; ++col) {
      EXPECT_DOUBLE_EQ(mat_m3_copy(row, col), mat_m3(row, col));
    }
  }
}

// --- Matrix3::array() ---

TEST_F(LinearAlgebraTests, Matrix3ArrayConversion) {
  // Column-major storage to row-major array conversion
  Vector3<double> const vec_col_0(1.0, 4.0, 7.0); // Column 0
  Vector3<double> const vec_col_1(2.0, 5.0, 8.0); // Column 1
  Vector3<double> const vec_col_2(3.0, 6.0, 9.0); // Column 2
  Matrix3<double> const mat_m(vec_col_0, vec_col_1, vec_col_2);

  auto arr = mat_m.array();

  // Row-major layout: arr[row][col]
  // Row 0: m(0,0)=1, m(0,1)=2, m(0,2)=3
  EXPECT_DOUBLE_EQ(arr[0][0], 1.0);
  EXPECT_DOUBLE_EQ(arr[0][1], 2.0);
  EXPECT_DOUBLE_EQ(arr[0][2], 3.0);

  // Row 1: m(1,0)=4, m(1,1)=5, m(1,2)=6
  EXPECT_DOUBLE_EQ(arr[1][0], 4.0);
  EXPECT_DOUBLE_EQ(arr[1][1], 5.0);
  EXPECT_DOUBLE_EQ(arr[1][2], 6.0);

  // Row 2: m(2,0)=7, m(2,1)=8, m(2,2)=9
  EXPECT_DOUBLE_EQ(arr[2][0], 7.0);
  EXPECT_DOUBLE_EQ(arr[2][1], 8.0);
  EXPECT_DOUBLE_EQ(arr[2][2], 9.0);
}

// --- Cross Product Properties ---

TEST_F(LinearAlgebraTests, CrossProductProperties) {
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

// --- Transpose Properties ---

TEST_F(LinearAlgebraTests, TransposeProperties) {
  Vector3<double> const vec_col_0(1.0, 2.0, 3.0);
  Vector3<double> const vec_col_1(4.0, 5.0, 6.0);
  Vector3<double> const vec_col_2(7.0, 8.0, 9.0);
  Matrix3<double> const mat_m(vec_col_0, vec_col_1, vec_col_2);

  // Double transpose equals original
  auto mat_m_transposed_twice = transpose(transpose(mat_m));
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t col = 0; col < 3; ++col) {
      EXPECT_DOUBLE_EQ(mat_m_transposed_twice(row, col), mat_m(row, col));
    }
  }

  // Transpose of identity is identity
  auto mat_identity = Matrix3<double>::identity();
  auto mat_identity_t = transpose(mat_identity);
  EXPECT_TRUE(mat_identity == mat_identity_t);

  // Symmetric matrix: transpose(A) == A
  Vector3<double> const vec_sym_col_0(1.0, 2.0, 3.0);
  Vector3<double> const vec_sym_col_1(2.0, 5.0, 6.0);
  Vector3<double> const vec_sym_col_2(3.0, 6.0, 9.0);
  Matrix3<double> const mat_symmetric(vec_sym_col_0, vec_sym_col_1, vec_sym_col_2);
  auto mat_sym_t = transpose(mat_symmetric);
  EXPECT_TRUE(mat_symmetric == mat_sym_t);

  // Trace is invariant under transpose
  EXPECT_DOUBLE_EQ(mat_m.trace(), transpose(mat_m).trace());

  // Determinant is invariant under transpose
  Vector3<double> const vec_det_col_0(2.0, 3.0, 1.0);
  Vector3<double> const vec_det_col_1(4.0, 1.0, 3.0);
  Vector3<double> const vec_det_col_2(1.0, 2.0, 5.0);
  Matrix3<double> const non_singular(vec_det_col_0, vec_det_col_1, vec_det_col_2);
  EXPECT_DOUBLE_EQ(determinant(non_singular), determinant(transpose(non_singular)));
}

// --- Invert Round-Trip (General Non-Diagonal Matrix) ---

TEST_F(LinearAlgebraTests, InvertRoundTripGeneral) {
  // A general 3x3 matrix (not diagonal, not symmetric)
  Vector3<double> const vec_col_0(2.0, 1.0, 0.0);
  Vector3<double> const vec_col_1(-1.0, 3.0, 2.0);
  Vector3<double> const vec_col_2(4.0, -2.0, 5.0);
  Matrix3<double> const mat_m(vec_col_0, vec_col_1, vec_col_2);

  // Verify it's invertible
  double const det = determinant(mat_m);
  EXPECT_GT(std::abs(det), 1e-10);

  auto mat_m_inv = invert(mat_m);

  // M * M^-1 = I
  auto product = mat_m * mat_m_inv;
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t col = 0; col < 3; ++col) {
      double const expected = (row == col) ? 1.0 : 0.0;
      EXPECT_NEAR(product(row, col), expected, 1e-14);
    }
  }

  // M^-1 * M = I (inverse on both sides)
  auto product2 = mat_m_inv * mat_m;
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t col = 0; col < 3; ++col) {
      double const expected = (row == col) ? 1.0 : 0.0;
      EXPECT_NEAR(product2(row, col), expected, 1e-14);
    }
  }

  // det(M^-1) = 1/det(M)
  EXPECT_NEAR(determinant(mat_m_inv), 1.0 / det, 1e-14);
}

// --- Non-Double Type Instantiation ---

TEST_F(LinearAlgebraTests, Vector3IntType) {
  Vector3<int> const vec_a(1, 2, 3);
  Vector3<int> const vec_b(4, 5, 6);

  // Basic arithmetic
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

  // Scalar multiply
  auto scaled = vec_a * 3;
  EXPECT_EQ(scaled.x(), 3);
  EXPECT_EQ(scaled.y(), 6);
  EXPECT_EQ(scaled.z(), 9);

  // Dot product
  EXPECT_EQ(dot(vec_a, vec_b), 32);

  // Cross product
  auto vec_cross = cross(vec_a, vec_b);
  EXPECT_EQ(vec_cross.x(), 2 * 6 - 3 * 5); // -3
  EXPECT_EQ(vec_cross.y(), 3 * 4 - 1 * 6); // 6
  EXPECT_EQ(vec_cross.z(), 1 * 5 - 2 * 4); // -3

  // Equality
  Vector3<int> const vec_a2(1, 2, 3);
  EXPECT_TRUE(vec_a == vec_a2);
  EXPECT_FALSE(vec_a != vec_a2);
  EXPECT_FALSE(vec_a == vec_b);

  // empty
  Vector3<int> const zero;
  EXPECT_TRUE(zero.empty());
  EXPECT_FALSE(vec_a.empty());
}

TEST_F(LinearAlgebraTests, Vector3FloatType) {
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

TEST_F(LinearAlgebraTests, Matrix3IntType) {
  Vector3<int> const vec_col_0(1, 0, 0);
  Vector3<int> const vec_col_1(0, 2, 0);
  Vector3<int> const vec_col_2(0, 0, 3);
  Matrix3<int> const mat_m(vec_col_0, vec_col_1, vec_col_2);

  EXPECT_EQ(mat_m.trace(), 6);
  EXPECT_EQ(determinant(mat_m), 6);

  // Matrix * vector
  Vector3<int> const vec_v(1, 2, 3);
  auto mat_m_v = mat_m * vec_v;
  EXPECT_EQ(mat_m_v.x(), 1);
  EXPECT_EQ(mat_m_v.y(), 4);
  EXPECT_EQ(mat_m_v.z(), 9);

  // Identity
  auto ident = Matrix3<int>::identity();
  EXPECT_EQ(ident(0, 0), 1);
  EXPECT_EQ(ident(1, 1), 1);
  EXPECT_EQ(ident(2, 2), 1);
  EXPECT_EQ(ident(0, 1), 0);

  // Equality
  EXPECT_TRUE(mat_m == mat_m);
  EXPECT_FALSE(mat_m == ident);
}

// --- Norm Properties ---

TEST_F(LinearAlgebraTests, NormProperties) {
  Vector3<double> const vec_v(3.0, 4.0, 0.0);

  // norm_sq == norm^2
  EXPECT_DOUBLE_EQ(norm_sq(vec_v), norm(vec_v) * norm(vec_v));

  // norm of zero vector is zero
  Vector3<double> const zero;
  EXPECT_DOUBLE_EQ(norm(zero), 0.0);
  EXPECT_DOUBLE_EQ(norm_sq(zero), 0.0);

  // norm is non-negative
  Vector3<double> const neg(-3.0, -4.0, -5.0);
  EXPECT_GE(norm(neg), 0.0);

  // Triangle inequality: |a + b| <= |a| + |b|
  Vector3<double> const vec_a(1.0, 2.0, 3.0);
  Vector3<double> const vec_b(4.0, -1.0, 2.0);
  EXPECT_LE(norm(vec_a + vec_b), norm(vec_a) + norm(vec_b) + 1e-15);

  // Scaling: |s*v| = |s|*|v|
  double const scalar = -3.5;
  EXPECT_NEAR(norm(vec_v * scalar), std::abs(scalar) * norm(vec_v), 1e-15);

  // Normalized vector has unit length
  auto normalized = normalize(vec_v);
  EXPECT_NEAR(norm(normalized), 1.0, 1e-15);
}

// --- Matrix-Vector Interaction Properties ---

TEST_F(LinearAlgebraTests, MatrixVectorDistributive) {
  Vector3<double> const vec_col_0(2.0, 1.0, 0.0);
  Vector3<double> const vec_col_1(-1.0, 3.0, 2.0);
  Vector3<double> const vec_col_2(4.0, -2.0, 5.0);
  Matrix3<double> const mat_m(vec_col_0, vec_col_1, vec_col_2);

  Vector3<double> const vec_v0(1.0, 2.0, 3.0);
  Vector3<double> const vec_v1(4.0, -1.0, 2.0);

  // M * (v1 + v2) = M * v1 + M * v2
  auto lhs = mat_m * (vec_v0 + vec_v1);
  auto rhs = (mat_m * vec_v0) + (mat_m * vec_v1);
  EXPECT_NEAR(lhs.x(), rhs.x(), 1e-15);
  EXPECT_NEAR(lhs.y(), rhs.y(), 1e-15);
  EXPECT_NEAR(lhs.z(), rhs.z(), 1e-15);

  // (A * B) * v = A * (B * v)  (associativity)
  Vector3<double> const vec_d0(1.0, 3.0, -1.0);
  Vector3<double> const vec_d1(2.0, 0.0, 4.0);
  Vector3<double> const vec_d2(-1.0, 2.0, 1.0);
  Matrix3<double> const mat_n(vec_d0, vec_d1, vec_d2);

  auto lhs2 = (mat_m * mat_n) * vec_v0;
  auto rhs2 = mat_m * (mat_n * vec_v0);
  EXPECT_NEAR(lhs2.x(), rhs2.x(), 1e-13);
  EXPECT_NEAR(lhs2.y(), rhs2.y(), 1e-13);
  EXPECT_NEAR(lhs2.z(), rhs2.z(), 1e-13);
}

} // namespace correlation::testing
