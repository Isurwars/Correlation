#ifndef INCLUDE_LINEAR_ALGEBRA_HPP_
#define INCLUDE_LINEAR_ALGEBRA_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include <array>

using Matrix3D = std::array<std::array<double, 3>, 3>;
using Vector3D = std::array<double, 3>;

//---------------------------------------------------------------------------//
//----------------------------- Vector Operators ----------------------------//
//---------------------------------------------------------------------------//

double dot(const Vector3D&, const Vector3D&);
double norm(const Vector3D&);

//---------------------------------------------------------------------------//
//----------------------------- Matrix Operators ----------------------------//
//---------------------------------------------------------------------------//

// Matrix Determinant for 3x3 matrices
double determinant(const Matrix3D &);
// Matrix inversion for 3x3 matrices
Matrix3D invertMatrix(const Matrix3D &mat);
// Matrix-vector multiplication
Vector3D matrixVectorMultiply(const Matrix3D &mat, const Vector3D &vec);

#endif // INCLUDE_LINEAR_ALGEBRA_HPP_
