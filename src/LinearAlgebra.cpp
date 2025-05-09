/* ----------------------------------------------------------------------------
 * Correlation: An Analysis Tool for Liquids and for Amorphous Solids
 * Copyright (c) 2013-2025 Isaías Rodríguez <isurwars@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the MIT License version as published in:
 * https://github.com/Isurwars/Correlation/blob/main/LICENSE
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 * ----------------------------------------------------------------------------
 */
#include "../include/LinearAlgebra.hpp"
#include <cmath>

Matrix3D invertMatrix(const Matrix3D &mat) {
  Matrix3D inv;
  const auto &m = mat.data;

  // Calculate determinant
  double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
	       m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
	       m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

  if (std::fabs(det) < 1e-15) {
    return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}; // Identity matrix fallback
  }

  double inv_det = 1.0 / det;

  inv.data[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det;
  inv.data[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) * inv_det;
  inv.data[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det;

  inv.data[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) * inv_det;
  inv.data[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det;
  inv.data[1][2] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) * inv_det;

  inv.data[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det;
  inv.data[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) * inv_det;
  inv.data[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det;

  return inv;
}

Vector3D matrixVectorMultiply(const Matrix3D &mat, const Vector3D &vec) {
  Vector3D result;
  for (int i = 0; i < 3; ++i) {
    result.data[i] = 0;
    for (int j = 0; j < 3; ++j) {
      result.data[i] += mat.data[i][j] * vec.data[j];
    }
  }
  return result;
}

Vector3D vectorSubtract(const Vector3D &a, const Vector3D &b) {
  return {
      {a.data[0] - b.data[0], a.data[1] - b.data[1], a.data[2] - b.data[2]}};
}
