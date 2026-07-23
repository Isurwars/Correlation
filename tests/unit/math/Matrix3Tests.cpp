// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "math/LinearAlgebra.hpp"

#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::math;
namespace {
class Matrix3Tests : public ::testing::Test {};
} // namespace

TEST_F(Matrix3Tests, Matrix3ConstructorsAndAccessors) {
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

TEST_F(Matrix3Tests, Matrix3Operations) {
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

TEST_F(Matrix3Tests, DeterminantAndInversion) {
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

  // Singular matrix determinant = 0
  Vector3<double> const c_deg0(1.0, 2.0, 3.0);
  Vector3<double> const c_deg1(2.0, 4.0, 6.0); // Linearly dependent
  Vector3<double> const c_deg2(0.0, 0.0, 1.0);
  Matrix3<double> const m_deg(c_deg0, c_deg1, c_deg2);

  EXPECT_NEAR(determinant(m_deg), 0.0, 1e-15);
  EXPECT_THROW((void)invert(m_deg), std::runtime_error);
}

TEST_F(Matrix3Tests, Transpose) {
  Vector3<double> const vec_col_0(1.0, 2.0, 3.0);
  Vector3<double> const vec_col_1(4.0, 5.0, 6.0);
  Vector3<double> const vec_col_2(7.0, 8.0, 9.0);
  Matrix3<double> const mat_m(vec_col_0, vec_col_1, vec_col_2);

  auto mat_t = transpose(mat_m);
  EXPECT_DOUBLE_EQ(mat_t(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(mat_t(1, 0), 4.0);
}

TEST_F(Matrix3Tests, Matrix3IdentityInverse) {
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

TEST_F(Matrix3Tests, Matrix3NegativeDeterminant) {
  Vector3<double> const vec_col_0(0.0, 1.0, 0.0);
  Vector3<double> const vec_col_1(1.0, 0.0, 0.0);
  Vector3<double> const vec_col_2(0.0, 0.0, 1.0);
  Matrix3<double> const mat_m(vec_col_0, vec_col_1, vec_col_2);

  EXPECT_DOUBLE_EQ(determinant(mat_m), -1.0);

  auto inv = invert(mat_m);
  auto product = mat_m * inv;
  EXPECT_NEAR(product(0, 0), 1.0, 1e-15);
  EXPECT_NEAR(product(1, 1), 1.0, 1e-15);
  EXPECT_NEAR(product(2, 2), 1.0, 1e-15);
  EXPECT_NEAR(product(0, 1), 0.0, 1e-15);
}

TEST_F(Matrix3Tests, Matrix3ExtremeValues) {
  double const big = 1e10;
  Vector3<double> const vec_col_0(big, 0.0, 0.0);
  Vector3<double> const vec_col_1(0.0, big, 0.0);
  Vector3<double> const vec_col_2(0.0, 0.0, big);
  Matrix3<double> const mat_big(vec_col_0, vec_col_1, vec_col_2);

  EXPECT_NEAR(determinant(mat_big), big * big * big, big * big * 1e-6);

  auto inv_big = invert(mat_big);
  EXPECT_NEAR(inv_big(0, 0), 1.0 / big, 1e-25);

  double const small = 1e-4;
  Vector3<double> const vec_small_0(small, 0.0, 0.0);
  Vector3<double> const vec_small_1(0.0, small, 0.0);
  Vector3<double> const vec_small_2(0.0, 0.0, small);
  Matrix3<double> const mat_small(vec_small_0, vec_small_1, vec_small_2);

  EXPECT_NEAR(determinant(mat_small), small * small * small, 1e-20);

  auto inv_small = invert(mat_small);
  EXPECT_NEAR(inv_small(0, 0), 1.0 / small, 1e-6);
}

TEST_F(Matrix3Tests, Matrix3AdditionSubtraction) {
  Vector3<double> const vec_col_0(1.0, 2.0, 3.0);
  Vector3<double> const vec_col_1(4.0, 5.0, 6.0);
  Vector3<double> const vec_col_2(7.0, 8.0, 9.0);
  Matrix3<double> const mat_m1(vec_col_0, vec_col_1, vec_col_2);

  Vector3<double> const vec_col_3(9.0, 8.0, 7.0);
  Vector3<double> const vec_col_4(6.0, 5.0, 4.0);
  Vector3<double> const vec_col_5(3.0, 2.0, 1.0);
  Matrix3<double> const mat_m2(vec_col_3, vec_col_4, vec_col_5);

  auto mat_sum = mat_m1 + mat_m2;
  EXPECT_DOUBLE_EQ(mat_sum(0, 0), 10.0);
  EXPECT_DOUBLE_EQ(mat_sum(1, 1), 10.0);
  EXPECT_DOUBLE_EQ(mat_sum(2, 2), 10.0);

  auto diff = mat_m1 - mat_m2;
  EXPECT_DOUBLE_EQ(diff(0, 0), -8.0);
  EXPECT_DOUBLE_EQ(diff(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(diff(2, 2), 8.0);

  Matrix3<double> mat_m3 = mat_m1;
  mat_m3 -= mat_m2;
  EXPECT_DOUBLE_EQ(mat_m3(0, 0), diff(0, 0));
  EXPECT_DOUBLE_EQ(mat_m3(1, 1), diff(1, 1));
  EXPECT_DOUBLE_EQ(mat_m3(2, 2), diff(2, 2));

  Matrix3<double> const mat_m1_copy = mat_m1;
  auto mat_zero = mat_m1 - mat_m1_copy;
  EXPECT_DOUBLE_EQ(mat_zero(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(mat_zero(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(mat_zero(2, 2), 0.0);
}

TEST_F(Matrix3Tests, Matrix3Identity) {
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

  Vector3<double> const vec_v(3.0, 7.0, -2.0);
  auto mat_result = mat_ident * vec_v;
  EXPECT_DOUBLE_EQ(mat_result.x(), vec_v.x());
  EXPECT_DOUBLE_EQ(mat_result.y(), vec_v.y());
  EXPECT_DOUBLE_EQ(mat_result.z(), vec_v.z());

  Vector3<double> const vec_col_0(1.0, 2.0, 3.0);
  Vector3<double> const vec_col_1(4.0, 5.0, 6.0);
  Vector3<double> const vec_col_2(7.0, 8.0, 9.0);
  Matrix3<double> const mat_m(vec_col_0, vec_col_1, vec_col_2);
  auto mat_product = mat_ident * mat_m;
  EXPECT_DOUBLE_EQ(mat_product(0, 0), mat_m(0, 0));
  EXPECT_DOUBLE_EQ(mat_product(1, 1), mat_m(1, 1));
  EXPECT_DOUBLE_EQ(mat_product(2, 2), mat_m(2, 2));
}

TEST_F(Matrix3Tests, Matrix3Equality) {
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

  Matrix3<double> const mat_ident = Matrix3<double>::identity();
  Matrix3<double> const mat_zero;
  EXPECT_FALSE(mat_ident == mat_zero);
  EXPECT_TRUE(mat_ident != mat_zero);
}

TEST_F(Matrix3Tests, Matrix3InPlaceAddition) {
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

TEST_F(Matrix3Tests, Matrix3InPlaceScalarMultiply) {
  Vector3<double> const vec_col_0(1.0, 0.0, 0.0);
  Vector3<double> const vec_col_1(0.0, 2.0, 0.0);
  Vector3<double> const vec_col_2(0.0, 0.0, 3.0);
  Matrix3<double> mat_m1(vec_col_0, vec_col_1, vec_col_2);

  Matrix3<double> mat_m2 = mat_m1;
  mat_m2 *= 4.0;
  EXPECT_DOUBLE_EQ(mat_m2(0, 0), 4.0);
  EXPECT_DOUBLE_EQ(mat_m2(1, 1), 8.0);
  EXPECT_DOUBLE_EQ(mat_m2(2, 2), 12.0);
  EXPECT_DOUBLE_EQ(mat_m2(0, 1), 0.0);

  mat_m2 *= 0.0;
  EXPECT_DOUBLE_EQ(mat_m2(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(mat_m2(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(mat_m2(2, 2), 0.0);
}

TEST_F(Matrix3Tests, Matrix3InPlaceMatrixMultiply) {
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

TEST_F(Matrix3Tests, Matrix3ArrayConversion) {
  Vector3<double> const vec_col_0(1.0, 4.0, 7.0);
  Vector3<double> const vec_col_1(2.0, 5.0, 8.0);
  Vector3<double> const vec_col_2(3.0, 6.0, 9.0);
  Matrix3<double> const mat_m(vec_col_0, vec_col_1, vec_col_2);

  auto arr = mat_m.array();
  EXPECT_DOUBLE_EQ(arr[0][0], 1.0);
  EXPECT_DOUBLE_EQ(arr[0][1], 2.0);
  EXPECT_DOUBLE_EQ(arr[0][2], 3.0);
  EXPECT_DOUBLE_EQ(arr[1][0], 4.0);
  EXPECT_DOUBLE_EQ(arr[1][1], 5.0);
  EXPECT_DOUBLE_EQ(arr[1][2], 6.0);
  EXPECT_DOUBLE_EQ(arr[2][0], 7.0);
  EXPECT_DOUBLE_EQ(arr[2][1], 8.0);
  EXPECT_DOUBLE_EQ(arr[2][2], 9.0);
}

TEST_F(Matrix3Tests, TransposeProperties) {
  Vector3<double> const vec_col_0(1.0, 2.0, 3.0);
  Vector3<double> const vec_col_1(4.0, 5.0, 6.0);
  Vector3<double> const vec_col_2(7.0, 8.0, 9.0);
  Matrix3<double> const mat_m(vec_col_0, vec_col_1, vec_col_2);

  auto mat_m_transposed_twice = transpose(transpose(mat_m));
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t col = 0; col < 3; ++col) {
      EXPECT_DOUBLE_EQ(mat_m_transposed_twice(row, col), mat_m(row, col));
    }
  }

  auto mat_identity = Matrix3<double>::identity();
  auto mat_identity_t = transpose(mat_identity);
  EXPECT_TRUE(mat_identity == mat_identity_t);

  Vector3<double> const vec_sym_col_0(1.0, 2.0, 3.0);
  Vector3<double> const vec_sym_col_1(2.0, 5.0, 6.0);
  Vector3<double> const vec_sym_col_2(3.0, 6.0, 9.0);
  Matrix3<double> const mat_symmetric(vec_sym_col_0, vec_sym_col_1, vec_sym_col_2);
  auto mat_sym_t = transpose(mat_symmetric);
  EXPECT_TRUE(mat_symmetric == mat_sym_t);

  EXPECT_DOUBLE_EQ(mat_m.trace(), transpose(mat_m).trace());

  Vector3<double> const vec_det_col_0(2.0, 3.0, 1.0);
  Vector3<double> const vec_det_col_1(4.0, 1.0, 3.0);
  Vector3<double> const vec_det_col_2(1.0, 2.0, 5.0);
  Matrix3<double> const non_singular(vec_det_col_0, vec_det_col_1, vec_det_col_2);
  EXPECT_DOUBLE_EQ(determinant(non_singular), determinant(transpose(non_singular)));
}

TEST_F(Matrix3Tests, InvertRoundTripGeneral) {
  Vector3<double> const vec_col_0(2.0, 1.0, 0.0);
  Vector3<double> const vec_col_1(-1.0, 3.0, 2.0);
  Vector3<double> const vec_col_2(4.0, -2.0, 5.0);
  Matrix3<double> const mat_m(vec_col_0, vec_col_1, vec_col_2);

  double const det = determinant(mat_m);
  EXPECT_GT(std::abs(det), 1e-10);

  auto mat_m_inv = invert(mat_m);

  auto product = mat_m * mat_m_inv;
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t col = 0; col < 3; ++col) {
      double const expected = (row == col) ? 1.0 : 0.0;
      EXPECT_NEAR(product(row, col), expected, 1e-14);
    }
  }

  auto product2 = mat_m_inv * mat_m;
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t col = 0; col < 3; ++col) {
      double const expected = (row == col) ? 1.0 : 0.0;
      EXPECT_NEAR(product2(row, col), expected, 1e-14);
    }
  }

  EXPECT_NEAR(determinant(mat_m_inv), 1.0 / det, 1e-14);
}

TEST_F(Matrix3Tests, Matrix3IntType) {
  Vector3<int> const vec_col_0(1, 0, 0);
  Vector3<int> const vec_col_1(0, 2, 0);
  Vector3<int> const vec_col_2(0, 0, 3);
  Matrix3<int> const mat_m(vec_col_0, vec_col_1, vec_col_2);

  EXPECT_EQ(mat_m.trace(), 6);
  EXPECT_EQ(determinant(mat_m), 6);

  Vector3<int> const vec_v(1, 2, 3);
  auto mat_m_v = mat_m * vec_v;
  EXPECT_EQ(mat_m_v.x(), 1);
  EXPECT_EQ(mat_m_v.y(), 4);
  EXPECT_EQ(mat_m_v.z(), 9);

  auto ident = Matrix3<int>::identity();
  EXPECT_EQ(ident(0, 0), 1);
  EXPECT_EQ(ident(1, 1), 1);
  EXPECT_EQ(ident(2, 2), 1);
  EXPECT_EQ(ident(0, 1), 0);

  EXPECT_TRUE(mat_m == mat_m);
  EXPECT_FALSE(mat_m == ident);
}

TEST_F(Matrix3Tests, MatrixVectorDistributive) {
  Vector3<double> const vec_col_0(2.0, 1.0, 0.0);
  Vector3<double> const vec_col_1(-1.0, 3.0, 2.0);
  Vector3<double> const vec_col_2(4.0, -2.0, 5.0);
  Matrix3<double> const mat_m(vec_col_0, vec_col_1, vec_col_2);

  Vector3<double> const vec_v0(1.0, 2.0, 3.0);
  Vector3<double> const vec_v1(4.0, -1.0, 2.0);

  auto lhs = mat_m * (vec_v0 + vec_v1);
  auto rhs = (mat_m * vec_v0) + (mat_m * vec_v1);
  EXPECT_NEAR(lhs.x(), rhs.x(), 1e-15);
  EXPECT_NEAR(lhs.y(), rhs.y(), 1e-15);
  EXPECT_NEAR(lhs.z(), rhs.z(), 1e-15);

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
