// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "math/LinearAlgebra.hpp"

#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::math;

class LinearAlgebraTests : public ::testing::Test {};

// --- Vector3 Tests ---

TEST_F(LinearAlgebraTests, Vector3ConstructorsAndAccessors) {
  // Default constructor
  Vector3<double> v1;
  EXPECT_DOUBLE_EQ(v1.x(), 0.0);
  EXPECT_DOUBLE_EQ(v1.y(), 0.0);
  EXPECT_DOUBLE_EQ(v1.z(), 0.0);
  EXPECT_TRUE(v1.empty());

  // Parameterized constructor
  Vector3<double> v2(1.5, -2.5, 3.0);
  EXPECT_DOUBLE_EQ(v2.x(), 1.5);
  EXPECT_DOUBLE_EQ(v2.y(), -2.5);
  EXPECT_DOUBLE_EQ(v2.z(), 3.0);
  EXPECT_FALSE(v2.empty());

  // Array constructor
  std::array<double, 3> const arr = {10.0, 20.0, 30.0};
  Vector3<double> v3(arr);
  EXPECT_DOUBLE_EQ(v3[0], 10.0);
  EXPECT_DOUBLE_EQ(v3[1], 20.0);
  EXPECT_DOUBLE_EQ(v3[2], 30.0);
  EXPECT_EQ(v3.array(), arr);

  // Mutable access
  v3.x() = 1.0;
  v3[1] = 2.0;
  v3(2) = 3.0;
  EXPECT_DOUBLE_EQ(v3.x(), 1.0);
  EXPECT_DOUBLE_EQ(v3.y(), 2.0);
  EXPECT_DOUBLE_EQ(v3.z(), 3.0);
}

TEST_F(LinearAlgebraTests, Vector3Arithmetic) {
  Vector3<double> const a(1.0, 2.0, 3.0);
  Vector3<double> const b(4.0, 5.0, 6.0);

  // Addition
  auto c = a + b;
  EXPECT_DOUBLE_EQ(c.x(), 5.0);
  EXPECT_DOUBLE_EQ(c.y(), 7.0);
  EXPECT_DOUBLE_EQ(c.z(), 9.0);

  // Subtraction
  auto d = a - b;
  EXPECT_DOUBLE_EQ(d.x(), -3.0);
  EXPECT_DOUBLE_EQ(d.y(), -3.0);
  EXPECT_DOUBLE_EQ(d.z(), -3.0);

  // Scalar multiplication / division
  auto e = a * 2.0;
  EXPECT_DOUBLE_EQ(e.x(), 2.0);
  EXPECT_DOUBLE_EQ(e.y(), 4.0);
  EXPECT_DOUBLE_EQ(e.z(), 6.0);

  auto f = b / 2.0;
  EXPECT_DOUBLE_EQ(f.x(), 2.0);
  EXPECT_DOUBLE_EQ(f.y(), 2.5);
  EXPECT_DOUBLE_EQ(f.z(), 3.0);

  // In-place operators
  Vector3<double> g = a;
  g += b;
  EXPECT_DOUBLE_EQ(g.x(), 5.0);

  g -= b;
  EXPECT_DOUBLE_EQ(g.x(), 1.0);

  // Dot product operator
  double const dot_ab = a * b;
  EXPECT_DOUBLE_EQ(dot_ab, 32.0);
}

TEST_F(LinearAlgebraTests, Vector3FreeFunctions) {
  Vector3<double> const a(1.0, 2.0, 3.0);
  Vector3<double> const b(4.0, -5.0, 6.0);

  // dot
  EXPECT_DOUBLE_EQ(dot(a, b), 12.0);

  // cross
  auto cr = cross(a, b);
  EXPECT_DOUBLE_EQ(cr.x(), 2.0 * 6.0 - 3.0 * (-5.0)); // 27
  EXPECT_DOUBLE_EQ(cr.y(), 3.0 * 4.0 - 1.0 * 6.0);    // 6
  EXPECT_DOUBLE_EQ(cr.z(), 1.0 * (-5.0) - 2.0 * 4.0); // -13

  // norm_sq and norm
  Vector3<double> const v(3.0, 4.0, 0.0);
  EXPECT_DOUBLE_EQ(norm_sq(v), 25.0);
  EXPECT_DOUBLE_EQ(norm(v), 5.0);

  // normalize
  auto vn = normalize(v);
  EXPECT_DOUBLE_EQ(vn.x(), 0.6);
  EXPECT_DOUBLE_EQ(vn.y(), 0.8);
  EXPECT_DOUBLE_EQ(vn.z(), 0.0);

  // normalize singular vector throws
  Vector3<double> const zero_v;
  EXPECT_THROW(normalize(zero_v), std::domain_error);

  Vector3<double> const tiny_v(1e-301, 1e-301, 1e-301);
  EXPECT_THROW(normalize(tiny_v), std::domain_error);
}

// --- Matrix3 Tests ---

TEST_F(LinearAlgebraTests, Matrix3ConstructorsAndAccessors) {
  // Default constructor (zeros)
  Matrix3<double> m1;
  EXPECT_DOUBLE_EQ(m1(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(m1.trace(), 0.0);

  // Column constructor
  Vector3<double> const c0(1.0, 2.0, 3.0);
  Vector3<double> const c1(4.0, 5.0, 6.0);
  Vector3<double> const c2(7.0, 8.0, 9.0);
  Matrix3<double> m2(c0, c1, c2);

  EXPECT_DOUBLE_EQ(m2(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(m2(1, 0), 2.0); // Row 1, Col 0
  EXPECT_DOUBLE_EQ(m2(0, 1), 4.0); // Row 0, Col 1
  EXPECT_DOUBLE_EQ(m2(2, 2), 9.0);
  EXPECT_DOUBLE_EQ(m2.trace(), 15.0);

  // Mutable access
  m2[0][0] = -1.0;
  EXPECT_DOUBLE_EQ(m2(0, 0), -1.0);
}

TEST_F(LinearAlgebraTests, Matrix3Operations) {
  Vector3<double> const c0(1.0, 0.0, 0.0);
  Vector3<double> const c1(0.0, 2.0, 0.0);
  Vector3<double> const c2(0.0, 0.0, 3.0);
  Matrix3<double> const m(c0, c1, c2);

  // Matrix * Vector
  Vector3<double> const v(1.0, 2.0, 3.0);
  auto mv = m * v;
  EXPECT_DOUBLE_EQ(mv.x(), 1.0);
  EXPECT_DOUBLE_EQ(mv.y(), 4.0);
  EXPECT_DOUBLE_EQ(mv.z(), 9.0);

  // Matrix * Matrix
  auto mm = m * m;
  EXPECT_DOUBLE_EQ(mm(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(mm(1, 1), 4.0);
  EXPECT_DOUBLE_EQ(mm(2, 2), 9.0);

  // Scalar multiplication
  auto m_scaled = m * 2.0;
  EXPECT_DOUBLE_EQ(m_scaled(1, 1), 4.0);
}

TEST_F(LinearAlgebraTests, DeterminantAndInversion) {
  // Simple orthogonal matrix
  Vector3<double> const c0(2.0, 0.0, 0.0);
  Vector3<double> const c1(0.0, 4.0, 0.0);
  Vector3<double> const c2(0.0, 0.0, 5.0);
  Matrix3<double> const m(c0, c1, c2);

  EXPECT_DOUBLE_EQ(determinant(m), 40.0);

  auto minv = invert(m);
  EXPECT_DOUBLE_EQ(minv(0, 0), 0.5);
  EXPECT_DOUBLE_EQ(minv(1, 1), 0.25);
  EXPECT_DOUBLE_EQ(minv(2, 2), 0.2);

  // singular matrix determinant = 0
  Vector3<double> const c_deg0(1.0, 2.0, 3.0);
  Vector3<double> const c_deg1(2.0, 4.0, 6.0); // Linearly dependent
  Vector3<double> const c_deg2(0.0, 0.0, 1.0);
  Matrix3<double> const m_deg(c_deg0, c_deg1, c_deg2);

  EXPECT_NEAR(determinant(m_deg), 0.0, 1e-15);
  EXPECT_THROW(invert(m_deg), std::runtime_error);
}

TEST_F(LinearAlgebraTests, Transpose) {
  Vector3<double> const c0(1.0, 2.0, 3.0);
  Vector3<double> const c1(4.0, 5.0, 6.0);
  Vector3<double> const c2(7.0, 8.0, 9.0);
  Matrix3<double> const m(c0, c1, c2);

  auto mt = transpose(m);
  EXPECT_DOUBLE_EQ(mt(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(mt(1, 0), 4.0);
}

// --- Extreme / Edge-Case Tests ---

TEST_F(LinearAlgebraTests, Vector3DivisionByZero) {
  Vector3<double> const v(1.0, 2.0, 3.0);
  auto result = v / 0.0;

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
  Vector3<double> const c0(1.0, 0.0, 0.0);
  Vector3<double> const c1(0.0, 1.0, 0.0);
  Vector3<double> const c2(0.0, 0.0, 1.0);
  Matrix3<double> const identity(c0, c1, c2);

  EXPECT_DOUBLE_EQ(determinant(identity), 1.0);

  auto inv = invert(identity);
  EXPECT_DOUBLE_EQ(inv(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(inv(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(inv(2, 2), 1.0);
  EXPECT_DOUBLE_EQ(inv(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(inv(0, 2), 0.0);
  EXPECT_DOUBLE_EQ(inv(1, 0), 0.0);
}

TEST_F(LinearAlgebraTests, Matrix3NegativeDeterminant) {
  // Left-handed coordinate system: negative determinant
  Vector3<double> const c0(0.0, 1.0, 0.0);
  Vector3<double> const c1(1.0, 0.0, 0.0);
  Vector3<double> const c2(0.0, 0.0, 1.0);
  Matrix3<double> const m(c0, c1, c2);

  EXPECT_DOUBLE_EQ(determinant(m), -1.0);

  // Should still be invertible
  auto inv = invert(m);
  // M * M^-1 should give identity
  auto product = m * inv;
  EXPECT_NEAR(product(0, 0), 1.0, 1e-15);
  EXPECT_NEAR(product(1, 1), 1.0, 1e-15);
  EXPECT_NEAR(product(2, 2), 1.0, 1e-15);
  EXPECT_NEAR(product(0, 1), 0.0, 1e-15);
}

TEST_F(LinearAlgebraTests, Matrix3ExtremeValues) {
  // Very large values
  double const big = 1e10;
  Vector3<double> const c0(big, 0.0, 0.0);
  Vector3<double> const c1(0.0, big, 0.0);
  Vector3<double> const c2(0.0, 0.0, big);
  Matrix3<double> const m_big(c0, c1, c2);

  EXPECT_NEAR(determinant(m_big), big * big * big, big * big * 1e-6);

  auto inv_big = invert(m_big);
  EXPECT_NEAR(inv_big(0, 0), 1.0 / big, 1e-25);

  // Very small values
  double const small = 1e-4;
  Vector3<double> const s0(small, 0.0, 0.0);
  Vector3<double> const s1(0.0, small, 0.0);
  Vector3<double> const s2(0.0, 0.0, small);
  Matrix3<double> const m_small(s0, s1, s2);

  EXPECT_NEAR(determinant(m_small), small * small * small, 1e-20);

  auto inv_small = invert(m_small);
  EXPECT_NEAR(inv_small(0, 0), 1.0 / small, 1e-6);
}

// --- New Operation Tests ---

TEST_F(LinearAlgebraTests, Vector3UnaryNegation) {
  Vector3<double> const v(1.0, -2.5, 3.0);
  auto neg = -v;
  EXPECT_DOUBLE_EQ(neg.x(), -1.0);
  EXPECT_DOUBLE_EQ(neg.y(), 2.5);
  EXPECT_DOUBLE_EQ(neg.z(), -3.0);

  // Double negation should return original
  auto double_neg = -(-v);
  EXPECT_DOUBLE_EQ(double_neg.x(), v.x());
  EXPECT_DOUBLE_EQ(double_neg.y(), v.y());
  EXPECT_DOUBLE_EQ(double_neg.z(), v.z());

  // Negation of zero vector
  Vector3<double> const zero;
  auto neg_zero = -zero;
  EXPECT_TRUE(neg_zero.empty());
}

TEST_F(LinearAlgebraTests, Vector3InPlaceScalarOps) {
  Vector3<double> v1(2.0, 4.0, 6.0);
  v1 *= 3.0;
  EXPECT_DOUBLE_EQ(v1.x(), 6.0);
  EXPECT_DOUBLE_EQ(v1.y(), 12.0);
  EXPECT_DOUBLE_EQ(v1.z(), 18.0);

  v1 /= 2.0;
  EXPECT_DOUBLE_EQ(v1.x(), 3.0);
  EXPECT_DOUBLE_EQ(v1.y(), 6.0);
  EXPECT_DOUBLE_EQ(v1.z(), 9.0);

  // Multiply by zero
  v1 *= 0.0;
  EXPECT_TRUE(v1.empty());
}

TEST_F(LinearAlgebraTests, DistanceFunction) {
  Vector3<double> const a(1.0, 0.0, 0.0);
  Vector3<double> const b(4.0, 0.0, 0.0);
  EXPECT_DOUBLE_EQ(distance(a, b), 3.0);

  // Distance is symmetric
  EXPECT_DOUBLE_EQ(distance(a, b), distance(b, a));

  // Distance to self is zero
  EXPECT_DOUBLE_EQ(distance(a, a), 0.0);

  // 3-4-5 triangle
  Vector3<double> const origin(0.0, 0.0, 0.0);
  Vector3<double> const point(3.0, 4.0, 0.0);
  EXPECT_DOUBLE_EQ(distance(origin, point), 5.0);
}

TEST_F(LinearAlgebraTests, Matrix3AdditionSubtraction) {
  Vector3<double> const c0(1.0, 2.0, 3.0);
  Vector3<double> const c1(4.0, 5.0, 6.0);
  Vector3<double> const c2(7.0, 8.0, 9.0);
  Matrix3<double> const m1(c0, c1, c2);

  Vector3<double> const d0(9.0, 8.0, 7.0);
  Vector3<double> const d1(6.0, 5.0, 4.0);
  Vector3<double> const d2(3.0, 2.0, 1.0);
  Matrix3<double> const m2(d0, d1, d2);

  // Addition
  auto sum = m1 + m2;
  EXPECT_DOUBLE_EQ(sum(0, 0), 10.0);
  EXPECT_DOUBLE_EQ(sum(1, 1), 10.0);
  EXPECT_DOUBLE_EQ(sum(2, 2), 10.0);

  // Subtraction
  auto diff = m1 - m2;
  EXPECT_DOUBLE_EQ(diff(0, 0), -8.0);
  EXPECT_DOUBLE_EQ(diff(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(diff(2, 2), 8.0);

  // In-place subtraction
  Matrix3<double> m3 = m1;
  m3 -= m2;
  EXPECT_DOUBLE_EQ(m3(0, 0), diff(0, 0));
  EXPECT_DOUBLE_EQ(m3(1, 1), diff(1, 1));
  EXPECT_DOUBLE_EQ(m3(2, 2), diff(2, 2));

  // M - M should be zero
  auto zero_matrix = m1 - m1;
  EXPECT_DOUBLE_EQ(zero_matrix(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(zero_matrix(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(zero_matrix(2, 2), 0.0);
}

TEST_F(LinearAlgebraTests, Matrix3Identity) {
  auto ident = Matrix3<double>::identity();
  EXPECT_DOUBLE_EQ(ident(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(ident(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(ident(2, 2), 1.0);
  EXPECT_DOUBLE_EQ(ident(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(ident(0, 2), 0.0);
  EXPECT_DOUBLE_EQ(ident(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(ident(1, 2), 0.0);
  EXPECT_DOUBLE_EQ(ident(2, 0), 0.0);
  EXPECT_DOUBLE_EQ(ident(2, 1), 0.0);
  EXPECT_DOUBLE_EQ(ident.trace(), 3.0);
  EXPECT_DOUBLE_EQ(determinant(ident), 1.0);

  // Identity * vector = vector
  Vector3<double> const v(3.0, 7.0, -2.0);
  auto result = ident * v;
  EXPECT_DOUBLE_EQ(result.x(), v.x());
  EXPECT_DOUBLE_EQ(result.y(), v.y());
  EXPECT_DOUBLE_EQ(result.z(), v.z());

  // Identity * matrix = matrix
  Vector3<double> const c0(1.0, 2.0, 3.0);
  Vector3<double> const c1(4.0, 5.0, 6.0);
  Vector3<double> const c2(7.0, 8.0, 9.0);
  Matrix3<double> const m(c0, c1, c2);
  auto product = ident * m;
  EXPECT_DOUBLE_EQ(product(0, 0), m(0, 0));
  EXPECT_DOUBLE_EQ(product(1, 1), m(1, 1));
  EXPECT_DOUBLE_EQ(product(2, 2), m(2, 2));
}

TEST_F(LinearAlgebraTests, Vector3Equality) {
  Vector3<double> const a(1.0, 2.0, 3.0);
  Vector3<double> const b(1.0, 2.0, 3.0);
  Vector3<double> const c(1.0, 2.0, 3.1);

  EXPECT_TRUE(a == b);
  EXPECT_FALSE(a != b);

  EXPECT_FALSE(a == c);
  EXPECT_TRUE(a != c);

  // Zero vectors
  Vector3<double> const z1;
  Vector3<double> const z2;
  EXPECT_TRUE(z1 == z2);
}

TEST_F(LinearAlgebraTests, Matrix3Equality) {
  Vector3<double> const c0(1.0, 2.0, 3.0);
  Vector3<double> const c1(4.0, 5.0, 6.0);
  Vector3<double> const c2(7.0, 8.0, 9.0);
  Matrix3<double> const m1(c0, c1, c2);
  Matrix3<double> const m2(c0, c1, c2);

  EXPECT_TRUE(m1 == m2);
  EXPECT_FALSE(m1 != m2);

  Matrix3<double> m3 = m1;
  m3(0, 0) = 99.0;
  EXPECT_FALSE(m1 == m3);
  EXPECT_TRUE(m1 != m3);

  // Identity comparisons
  auto ident = Matrix3<double>::identity();
  Matrix3<double> const zero;
  EXPECT_FALSE(ident == zero);
  EXPECT_TRUE(ident != zero);
}

// --- Iterator Tests ---

TEST_F(LinearAlgebraTests, Vector3BeginEnd) {
  Vector3<double> const v(10.0, 20.0, 30.0);

  // Verify pointer-based iteration
  const double *iter = v.begin();
  EXPECT_DOUBLE_EQ(*iter, 10.0);
  ++iter;
  EXPECT_DOUBLE_EQ(*iter, 20.0);
  ++iter;
  EXPECT_DOUBLE_EQ(*iter, 30.0);
  ++iter;
  EXPECT_EQ(iter, v.end());

  // Range-for loop reads
  double sum = 0.0;
  for (double val : v) {
    sum += val;
  }
  EXPECT_DOUBLE_EQ(sum, 60.0);

  // Mutable iteration
  Vector3<double> v2(1.0, 2.0, 3.0);
  for (double &val : v2) {
    val *= 10.0;
  }
  EXPECT_DOUBLE_EQ(v2.x(), 10.0);
  EXPECT_DOUBLE_EQ(v2.y(), 20.0);
  EXPECT_DOUBLE_EQ(v2.z(), 30.0);

  // std::distance between begin and end
  EXPECT_EQ(v.end() - v.begin(), 3);
}

// --- Commutative Scalar-Vector Multiplication ---

TEST_F(LinearAlgebraTests, ScalarTimesVector) {
  Vector3<double> const v(2.0, 3.0, 4.0);

  // scalar * vector should equal vector * scalar
  auto lhs = 5.0 * v;
  auto rhs = v * 5.0;
  EXPECT_DOUBLE_EQ(lhs.x(), rhs.x());
  EXPECT_DOUBLE_EQ(lhs.y(), rhs.y());
  EXPECT_DOUBLE_EQ(lhs.z(), rhs.z());

  // Works with int scalar (cross-type via requires clause)
  auto int_result = 3 * v;
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
  Vector3<double> const c0(1.0, 2.0, 3.0);
  Vector3<double> const c1(4.0, 5.0, 6.0);
  Vector3<double> const c2(7.0, 8.0, 9.0);
  Matrix3<double> m1(c0, c1, c2);

  Vector3<double> const d0(10.0, 20.0, 30.0);
  Vector3<double> const d1(40.0, 50.0, 60.0);
  Vector3<double> const d2(70.0, 80.0, 90.0);
  Matrix3<double> const m2(d0, d1, d2);

  m1 += m2;
  EXPECT_DOUBLE_EQ(m1(0, 0), 11.0);
  EXPECT_DOUBLE_EQ(m1(1, 1), 55.0);
  EXPECT_DOUBLE_EQ(m1(2, 2), 99.0);
  EXPECT_DOUBLE_EQ(m1(0, 1), 44.0);
  EXPECT_DOUBLE_EQ(m1(2, 0), 33.0);
}

TEST_F(LinearAlgebraTests, Matrix3InPlaceScalarMultiply) {
  Vector3<double> const c0(1.0, 0.0, 0.0);
  Vector3<double> const c1(0.0, 2.0, 0.0);
  Vector3<double> const c2(0.0, 0.0, 3.0);
  Matrix3<double> m(c0, c1, c2);

  m *= 4.0;
  EXPECT_DOUBLE_EQ(m(0, 0), 4.0);
  EXPECT_DOUBLE_EQ(m(1, 1), 8.0);
  EXPECT_DOUBLE_EQ(m(2, 2), 12.0);
  EXPECT_DOUBLE_EQ(m(0, 1), 0.0); // Off-diagonal stays zero

  // Multiply by zero
  m *= 0.0;
  EXPECT_DOUBLE_EQ(m(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(m(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(m(2, 2), 0.0);
}

TEST_F(LinearAlgebraTests, Matrix3InPlaceMatrixMultiply) {
  // Diagonal matrices: multiplication is element-wise on diagonal
  Vector3<double> const c0(2.0, 0.0, 0.0);
  Vector3<double> const c1(0.0, 3.0, 0.0);
  Vector3<double> const c2(0.0, 0.0, 5.0);
  Matrix3<double> m(c0, c1, c2);

  Vector3<double> const d0(4.0, 0.0, 0.0);
  Vector3<double> const d1(0.0, 2.0, 0.0);
  Vector3<double> const d2(0.0, 0.0, 3.0);
  Matrix3<double> const m2(d0, d1, d2);

  m *= m2;
  EXPECT_DOUBLE_EQ(m(0, 0), 8.0);
  EXPECT_DOUBLE_EQ(m(1, 1), 6.0);
  EXPECT_DOUBLE_EQ(m(2, 2), 15.0);
  EXPECT_DOUBLE_EQ(m(0, 1), 0.0);

  // Multiply by identity leaves matrix unchanged
  Vector3<double> const e0(1.0, 2.0, 3.0);
  Vector3<double> const e1(4.0, 5.0, 6.0);
  Vector3<double> const e2(7.0, 8.0, 9.0);
  Matrix3<double> m3(e0, e1, e2);
  Matrix3<double> const orig = m3;
  m3 *= Matrix3<double>::identity();
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t col = 0; col < 3; ++col) {
      EXPECT_DOUBLE_EQ(m3(row, col), orig(row, col));
    }
  }
}

// --- Matrix3::array() ---

TEST_F(LinearAlgebraTests, Matrix3ArrayConversion) {
  // Column-major storage to row-major array conversion
  Vector3<double> const c0(1.0, 4.0, 7.0); // Column 0
  Vector3<double> const c1(2.0, 5.0, 8.0); // Column 1
  Vector3<double> const c2(3.0, 6.0, 9.0); // Column 2
  Matrix3<double> const m(c0, c1, c2);

  auto arr = m.array();

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
  Vector3<double> const a(1.0, 2.0, 3.0);
  Vector3<double> const b(4.0, -5.0, 6.0);

  // Anti-commutativity: a x b = -(b x a)
  auto axb = cross(a, b);
  auto bxa = cross(b, a);
  EXPECT_DOUBLE_EQ(axb.x(), -bxa.x());
  EXPECT_DOUBLE_EQ(axb.y(), -bxa.y());
  EXPECT_DOUBLE_EQ(axb.z(), -bxa.z());

  // Orthogonality: (a x b) . a = 0 and (a x b) . b = 0
  EXPECT_NEAR(dot(axb, a), 0.0, 1e-15);
  EXPECT_NEAR(dot(axb, b), 0.0, 1e-15);

  // Self cross product is zero
  auto axa = cross(a, a);
  EXPECT_NEAR(axa.x(), 0.0, 1e-15);
  EXPECT_NEAR(axa.y(), 0.0, 1e-15);
  EXPECT_NEAR(axa.z(), 0.0, 1e-15);

  // Basis vectors: x × y = z
  Vector3<double> const ex(1.0, 0.0, 0.0);
  Vector3<double> const ey(0.0, 1.0, 0.0);
  Vector3<double> const ez(0.0, 0.0, 1.0);
  auto xy = cross(ex, ey);
  EXPECT_DOUBLE_EQ(xy.x(), 0.0);
  EXPECT_DOUBLE_EQ(xy.y(), 0.0);
  EXPECT_DOUBLE_EQ(xy.z(), 1.0);

  // Lagrange's identity: |a x b|^2 = |a|^2*|b|^2 - (a.b)^2
  double const lhs = norm_sq(axb);
  double const rhs = norm_sq(a) * norm_sq(b) - dot(a, b) * dot(a, b);
  EXPECT_NEAR(lhs, rhs, 1e-10);
}

// --- Transpose Properties ---

TEST_F(LinearAlgebraTests, TransposeProperties) {
  Vector3<double> const c0(1.0, 2.0, 3.0);
  Vector3<double> const c1(4.0, 5.0, 6.0);
  Vector3<double> const c2(7.0, 8.0, 9.0);
  Matrix3<double> const m(c0, c1, c2);

  // Double transpose equals original
  auto mtt = transpose(transpose(m));
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t col = 0; col < 3; ++col) {
      EXPECT_DOUBLE_EQ(mtt(row, col), m(row, col));
    }
  }

  // Transpose of identity is identity
  auto ident = Matrix3<double>::identity();
  auto ident_t = transpose(ident);
  EXPECT_TRUE(ident == ident_t);

  // Symmetric matrix: transpose(A) == A
  Vector3<double> const s0(1.0, 2.0, 3.0);
  Vector3<double> const s1(2.0, 5.0, 6.0);
  Vector3<double> const s2(3.0, 6.0, 9.0);
  Matrix3<double> const sym(s0, s1, s2);
  auto sym_t = transpose(sym);
  EXPECT_TRUE(sym == sym_t);

  // Trace is invariant under transpose
  EXPECT_DOUBLE_EQ(m.trace(), transpose(m).trace());

  // Determinant is invariant under transpose
  Vector3<double> const r0(2.0, 3.0, 1.0);
  Vector3<double> const r1(4.0, 1.0, 3.0);
  Vector3<double> const r2(1.0, 2.0, 5.0);
  Matrix3<double> const non_singular(r0, r1, r2);
  EXPECT_DOUBLE_EQ(determinant(non_singular), determinant(transpose(non_singular)));
}

// --- Invert Round-Trip (General Non-Diagonal Matrix) ---

TEST_F(LinearAlgebraTests, InvertRoundTripGeneral) {
  // A general 3x3 matrix (not diagonal, not symmetric)
  Vector3<double> const c0(2.0, 1.0, 0.0);
  Vector3<double> const c1(-1.0, 3.0, 2.0);
  Vector3<double> const c2(4.0, -2.0, 5.0);
  Matrix3<double> const m(c0, c1, c2);

  // Verify it's invertible
  double const det = determinant(m);
  EXPECT_GT(std::abs(det), 1e-10);

  auto minv = invert(m);

  // M * M^-1 = I
  auto product = m * minv;
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t col = 0; col < 3; ++col) {
      double const expected = (row == col) ? 1.0 : 0.0;
      EXPECT_NEAR(product(row, col), expected, 1e-14);
    }
  }

  // M^-1 * M = I (inverse on both sides)
  auto product2 = minv * m;
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t col = 0; col < 3; ++col) {
      double const expected = (row == col) ? 1.0 : 0.0;
      EXPECT_NEAR(product2(row, col), expected, 1e-14);
    }
  }

  // det(M^-1) = 1/det(M)
  EXPECT_NEAR(determinant(minv), 1.0 / det, 1e-14);
}

// --- Non-Double Type Instantiation ---

TEST_F(LinearAlgebraTests, Vector3IntType) {
  Vector3<int> const a(1, 2, 3);
  Vector3<int> const b(4, 5, 6);

  // Basic arithmetic
  auto sum = a + b;
  EXPECT_EQ(sum.x(), 5);
  EXPECT_EQ(sum.y(), 7);
  EXPECT_EQ(sum.z(), 9);

  auto diff = a - b;
  EXPECT_EQ(diff.x(), -3);
  EXPECT_EQ(diff.y(), -3);
  EXPECT_EQ(diff.z(), -3);

  auto neg = -a;
  EXPECT_EQ(neg.x(), -1);
  EXPECT_EQ(neg.y(), -2);
  EXPECT_EQ(neg.z(), -3);

  // Scalar multiply
  auto scaled = a * 3;
  EXPECT_EQ(scaled.x(), 3);
  EXPECT_EQ(scaled.y(), 6);
  EXPECT_EQ(scaled.z(), 9);

  // Dot product
  EXPECT_EQ(dot(a, b), 32);

  // Cross product
  auto cr = cross(a, b);
  EXPECT_EQ(cr.x(), 2 * 6 - 3 * 5);  // -3
  EXPECT_EQ(cr.y(), 3 * 4 - 1 * 6);  // 6
  EXPECT_EQ(cr.z(), 1 * 5 - 2 * 4);  // -3

  // Equality
  Vector3<int> const a2(1, 2, 3);
  EXPECT_TRUE(a == a2);
  EXPECT_FALSE(a != a2);
  EXPECT_FALSE(a == b);

  // empty
  Vector3<int> const zero;
  EXPECT_TRUE(zero.empty());
  EXPECT_FALSE(a.empty());
}

TEST_F(LinearAlgebraTests, Vector3FloatType) {
  Vector3<float> const a(1.0F, 2.0F, 3.0F);
  Vector3<float> const b(4.0F, 5.0F, 6.0F);

  auto sum = a + b;
  EXPECT_FLOAT_EQ(sum.x(), 5.0F);
  EXPECT_FLOAT_EQ(sum.y(), 7.0F);
  EXPECT_FLOAT_EQ(sum.z(), 9.0F);

  EXPECT_FLOAT_EQ(dot(a, b), 32.0F);
  EXPECT_FLOAT_EQ(norm_sq(a), 14.0F);

  auto normalized = normalize(a);
  EXPECT_NEAR(norm(normalized), 1.0F, 1e-6F);
}

TEST_F(LinearAlgebraTests, Matrix3IntType) {
  Vector3<int> const c0(1, 0, 0);
  Vector3<int> const c1(0, 2, 0);
  Vector3<int> const c2(0, 0, 3);
  Matrix3<int> const m(c0, c1, c2);

  EXPECT_EQ(m.trace(), 6);
  EXPECT_EQ(determinant(m), 6);

  // Matrix * vector
  Vector3<int> const v(1, 2, 3);
  auto mv = m * v;
  EXPECT_EQ(mv.x(), 1);
  EXPECT_EQ(mv.y(), 4);
  EXPECT_EQ(mv.z(), 9);

  // Identity
  auto ident = Matrix3<int>::identity();
  EXPECT_EQ(ident(0, 0), 1);
  EXPECT_EQ(ident(1, 1), 1);
  EXPECT_EQ(ident(2, 2), 1);
  EXPECT_EQ(ident(0, 1), 0);

  // Equality
  EXPECT_TRUE(m == m);
  EXPECT_FALSE(m == ident);
}

// --- Norm Properties ---

TEST_F(LinearAlgebraTests, NormProperties) {
  Vector3<double> const v(3.0, 4.0, 0.0);

  // norm_sq == norm^2
  EXPECT_DOUBLE_EQ(norm_sq(v), norm(v) * norm(v));

  // norm of zero vector is zero
  Vector3<double> const zero;
  EXPECT_DOUBLE_EQ(norm(zero), 0.0);
  EXPECT_DOUBLE_EQ(norm_sq(zero), 0.0);

  // norm is non-negative
  Vector3<double> const neg(-3.0, -4.0, -5.0);
  EXPECT_GE(norm(neg), 0.0);

  // Triangle inequality: |a + b| <= |a| + |b|
  Vector3<double> const a(1.0, 2.0, 3.0);
  Vector3<double> const b(4.0, -1.0, 2.0);
  EXPECT_LE(norm(a + b), norm(a) + norm(b) + 1e-15);

  // Scaling: |s*v| = |s|*|v|
  double const scalar = -3.5;
  EXPECT_NEAR(norm(v * scalar), std::abs(scalar) * norm(v), 1e-15);

  // Normalized vector has unit length
  auto vn = normalize(v);
  EXPECT_NEAR(norm(vn), 1.0, 1e-15);
}

// --- Matrix-Vector Interaction Properties ---

TEST_F(LinearAlgebraTests, MatrixVectorDistributive) {
  Vector3<double> const c0(2.0, 1.0, 0.0);
  Vector3<double> const c1(-1.0, 3.0, 2.0);
  Vector3<double> const c2(4.0, -2.0, 5.0);
  Matrix3<double> const m(c0, c1, c2);

  Vector3<double> const v1(1.0, 2.0, 3.0);
  Vector3<double> const v2(4.0, -1.0, 2.0);

  // M * (v1 + v2) = M * v1 + M * v2
  auto lhs = m * (v1 + v2);
  auto rhs = (m * v1) + (m * v2);
  EXPECT_NEAR(lhs.x(), rhs.x(), 1e-15);
  EXPECT_NEAR(lhs.y(), rhs.y(), 1e-15);
  EXPECT_NEAR(lhs.z(), rhs.z(), 1e-15);

  // (A * B) * v = A * (B * v)  (associativity)
  Vector3<double> const d0(1.0, 3.0, -1.0);
  Vector3<double> const d1(2.0, 0.0, 4.0);
  Vector3<double> const d2(-1.0, 2.0, 1.0);
  Matrix3<double> const n(d0, d1, d2);

  auto lhs2 = (m * n) * v1;
  auto rhs2 = m * (n * v1);
  EXPECT_NEAR(lhs2.x(), rhs2.x(), 1e-13);
  EXPECT_NEAR(lhs2.y(), rhs2.y(), 1e-13);
  EXPECT_NEAR(lhs2.z(), rhs2.z(), 1e-13);
}

} // namespace correlation::testing

