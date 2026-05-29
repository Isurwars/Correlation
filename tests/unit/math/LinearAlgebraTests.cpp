// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
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
  std::array<double, 3> arr = {10.0, 20.0, 30.0};
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
  Vector3<double> a(1.0, 2.0, 3.0);
  Vector3<double> b(4.0, 5.0, 6.0);

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
  double dot_ab = a * b;
  EXPECT_DOUBLE_EQ(dot_ab, 32.0);
}

TEST_F(LinearAlgebraTests, Vector3FreeFunctions) {
  Vector3<double> a(1.0, 2.0, 3.0);
  Vector3<double> b(4.0, -5.0, 6.0);

  // dot
  EXPECT_DOUBLE_EQ(dot(a, b), 12.0);

  // cross
  auto cr = cross(a, b);
  EXPECT_DOUBLE_EQ(cr.x(), 2.0 * 6.0 - 3.0 * (-5.0)); // 27
  EXPECT_DOUBLE_EQ(cr.y(), 3.0 * 4.0 - 1.0 * 6.0);    // 6
  EXPECT_DOUBLE_EQ(cr.z(), 1.0 * (-5.0) - 2.0 * 4.0); // -13

  // norm_sq and norm
  Vector3<double> v(3.0, 4.0, 0.0);
  EXPECT_DOUBLE_EQ(norm_sq(v), 25.0);
  EXPECT_DOUBLE_EQ(norm(v), 5.0);

  // normalize
  auto vn = normalize(v);
  EXPECT_DOUBLE_EQ(vn.x(), 0.6);
  EXPECT_DOUBLE_EQ(vn.y(), 0.8);
  EXPECT_DOUBLE_EQ(vn.z(), 0.0);

  // normalize singular vector throws
  Vector3<double> zero_v;
  EXPECT_THROW(normalize(zero_v), std::domain_error);

  Vector3<double> tiny_v(1e-301, 1e-301, 1e-301);
  EXPECT_THROW(normalize(tiny_v), std::domain_error);
}

// --- Matrix3 Tests ---

TEST_F(LinearAlgebraTests, Matrix3ConstructorsAndAccessors) {
  // Default constructor (zeros)
  Matrix3<double> m1;
  EXPECT_DOUBLE_EQ(m1(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(m1.trace(), 0.0);

  // Column constructor
  Vector3<double> c0(1.0, 2.0, 3.0);
  Vector3<double> c1(4.0, 5.0, 6.0);
  Vector3<double> c2(7.0, 8.0, 9.0);
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
  Vector3<double> c0(1.0, 0.0, 0.0);
  Vector3<double> c1(0.0, 2.0, 0.0);
  Vector3<double> c2(0.0, 0.0, 3.0);
  Matrix3<double> m(c0, c1, c2);

  // Matrix * Vector
  Vector3<double> v(1.0, 2.0, 3.0);
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
  Vector3<double> c0(2.0, 0.0, 0.0);
  Vector3<double> c1(0.0, 4.0, 0.0);
  Vector3<double> c2(0.0, 0.0, 5.0);
  Matrix3<double> m(c0, c1, c2);

  EXPECT_DOUBLE_EQ(determinant(m), 40.0);

  auto minv = invert(m);
  EXPECT_DOUBLE_EQ(minv(0, 0), 0.5);
  EXPECT_DOUBLE_EQ(minv(1, 1), 0.25);
  EXPECT_DOUBLE_EQ(minv(2, 2), 0.2);

  // singular matrix determinant = 0
  Vector3<double> c_deg0(1.0, 2.0, 3.0);
  Vector3<double> c_deg1(2.0, 4.0, 6.0); // Linearly dependent
  Vector3<double> c_deg2(0.0, 0.0, 1.0);
  Matrix3<double> m_deg(c_deg0, c_deg1, c_deg2);

  EXPECT_NEAR(determinant(m_deg), 0.0, 1e-15);
  EXPECT_THROW(invert(m_deg), std::runtime_error);
}

TEST_F(LinearAlgebraTests, Transpose) {
  Vector3<double> c0(1.0, 2.0, 3.0);
  Vector3<double> c1(4.0, 5.0, 6.0);
  Vector3<double> c2(7.0, 8.0, 9.0);
  Matrix3<double> m(c0, c1, c2);

  auto mt = transpose(m);
  EXPECT_DOUBLE_EQ(mt(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(mt(1, 0), 4.0);
}

// --- Extreme / Edge-Case Tests ---

TEST_F(LinearAlgebraTests, Vector3DivisionByZero) {
  Vector3<double> v(1.0, 2.0, 3.0);
  auto result = v / 0.0;

  // Division by zero should produce infinity
  EXPECT_TRUE(std::isinf(result.x()));
  EXPECT_TRUE(std::isinf(result.y()));
  EXPECT_TRUE(std::isinf(result.z()));

  // Zero / zero should produce NaN
  Vector3<double> zero_v;
  auto nan_result = zero_v / 0.0;
  EXPECT_TRUE(std::isnan(nan_result.x()));
}

TEST_F(LinearAlgebraTests, Matrix3IdentityInverse) {
  // Identity matrix inverse should be identity
  Vector3<double> c0(1.0, 0.0, 0.0);
  Vector3<double> c1(0.0, 1.0, 0.0);
  Vector3<double> c2(0.0, 0.0, 1.0);
  Matrix3<double> identity(c0, c1, c2);

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
  Vector3<double> c0(0.0, 1.0, 0.0);
  Vector3<double> c1(1.0, 0.0, 0.0);
  Vector3<double> c2(0.0, 0.0, 1.0);
  Matrix3<double> m(c0, c1, c2);

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
  double big = 1e10;
  Vector3<double> c0(big, 0.0, 0.0);
  Vector3<double> c1(0.0, big, 0.0);
  Vector3<double> c2(0.0, 0.0, big);
  Matrix3<double> m_big(c0, c1, c2);

  EXPECT_NEAR(determinant(m_big), big * big * big, big * big * 1e-6);

  auto inv_big = invert(m_big);
  EXPECT_NEAR(inv_big(0, 0), 1.0 / big, 1e-25);

  // Very small values
  double small = 1e-4;
  Vector3<double> s0(small, 0.0, 0.0);
  Vector3<double> s1(0.0, small, 0.0);
  Vector3<double> s2(0.0, 0.0, small);
  Matrix3<double> m_small(s0, s1, s2);

  EXPECT_NEAR(determinant(m_small), small * small * small, 1e-20);

  auto inv_small = invert(m_small);
  EXPECT_NEAR(inv_small(0, 0), 1.0 / small, 1e-6);
}

} // namespace correlation::testing
