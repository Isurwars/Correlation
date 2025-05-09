#ifndef INCLUDE_LINEAR_ALGEBRA_H_
#define INCLUDE_LINEAR_ALGEBRA_H_
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
#include <array>

struct Matrix3D {
  double data[3][3];
};

struct Vector3D {
  double data[3];
};

// Matrix inversion for 3x3 matrices
Matrix3D invertMatrix(const Matrix3D &mat);

// Matrix-vector multiplication
Vector3D matrixVectorMultiply(const Matrix3D &mat, const Vector3D &vec);

// Vector subtraction
Vector3D vectorSubtract(const Vector3D &a, const Vector3D &b);

#endif // INCLUDE_LINEAR_ALGEBRA_H_
