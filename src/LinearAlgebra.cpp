// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/LinearAlgebra.hpp"

#include <cmath>

//---------------------------------------------------------------------------//
//----------------------------- Vector Operators ----------------------------//
//---------------------------------------------------------------------------//

Vector3D operator+(const Vector3D &a, const Vector3D &b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Vector3D operator-(const Vector3D &a, const Vector3D &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

double operator*(const Vector3D &a, const Vector3D &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double dot(const Vector3D &a, const Vector3D &b) { return a * b; }

double norm(const Vector3D &a) { return std::sqrt(dot(a, a)); }

Vector3D cross(const Vector3D &a, const Vector3D &b) {
  return {a[1] * b[2] - a[2] * b[1], a[0] * b[2] - a[2] * b[0],
          a[0] * b[1] - a[1] * b[0]};
}

//---------------------------------------------------------------------------//
//----------------------------- Matrix Operators ----------------------------//
//---------------------------------------------------------------------------//

double determinant(const Matrix3D &m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
} // determinant

Matrix3D invertMatrix(const Matrix3D &m) {
  Matrix3D inv;

  // Calculate determinant
  double det = determinant(m);

  if (std::fabs(det) < 1e-15) {
    return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
  }

  double inv_det = 1.0 / det;

  inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det;
  inv[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) * inv_det;
  inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det;

  inv[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) * inv_det;
  inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det;
  inv[1][2] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) * inv_det;

  inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det;
  inv[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) * inv_det;
  inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det;

  return inv;
} // invertMatrix

Vector3D matrixVectorMultiply(const Matrix3D &mat, const Vector3D &vec) {
  Vector3D result;

  for (int i = 0; i < 3; ++i) {
    result.data()[i] = 0;
    for (int j = 0; j < 3; ++j) {
      result[i] += mat[i][j] * vec[j];
    }
  }

  return result;
} // matrixVectorMultiply
