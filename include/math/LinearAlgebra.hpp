/**
 * @file LinearAlgebra.hpp
 * @brief Lightweight 3D vector and matrix classes for geometric calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @license SPDX-License-Identifier: MIT
 */

#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace correlation::math {

// -----------------------------------------------------------------------------
//  Vector3  –  lightweight, constexpr, stack-based
// -----------------------------------------------------------------------------
/**
 * @class Vector3
 * @brief A lightweight, constexpr, stack-based 3D vector.
 */
template <typename T> class Vector3 {
public:
  using value_type = T;

  constexpr Vector3() noexcept : data_{0, 0, 0} {}
  constexpr Vector3(T x, T y, T z) noexcept : data_{x, y, z} {}
  constexpr explicit Vector3(const std::array<T, 3> &a) noexcept : data_{a} {}

  // Access by index
  constexpr T operator[](std::size_t i) const noexcept { return data_[i]; }
  constexpr T &operator[](std::size_t i) noexcept { return data_[i]; }

  // Named component accessors
  constexpr T x() const noexcept { return data_[0]; }
  constexpr T &x() noexcept { return data_[0]; }
  constexpr T y() const noexcept { return data_[1]; }
  constexpr T &y() noexcept { return data_[1]; }
  constexpr T z() const noexcept { return data_[2]; }
  constexpr T &z() noexcept { return data_[2]; }

  constexpr T &operator()(std::size_t i) noexcept { return data_[i]; }
  constexpr T operator()(std::size_t i) const noexcept { return data_[i]; }

  // True if every component is exactly zero
  constexpr bool empty() const noexcept {
    return data_[0] == T{} && data_[1] == T{} && data_[2] == T{};
  }

  constexpr const T *begin() const noexcept { return data_.data(); }
  constexpr const T *end() const noexcept { return data_.data() + 3; }
  constexpr T *begin() noexcept { return data_.data(); }
  constexpr T *end() noexcept { return data_.data() + 3; }

  // arithmetic
  constexpr Vector3 operator+(const Vector3 &rhs) const noexcept {
    return {data_[0] + rhs[0], data_[1] + rhs[1], data_[2] + rhs[2]};
  }
  constexpr Vector3 operator-(const Vector3 &rhs) const noexcept {
    return {data_[0] - rhs[0], data_[1] - rhs[1], data_[2] - rhs[2]};
  }
  constexpr Vector3 operator*(T s) const noexcept {
    return {data_[0] * s, data_[1] * s, data_[2] * s};
  }
  constexpr Vector3 operator/(T s) const noexcept {
    return {data_[0] / s, data_[1] / s, data_[2] / s};
  }

  // compound
  constexpr Vector3 &operator+=(const Vector3 &rhs) noexcept {
    data_[0] += rhs[0];
    data_[1] += rhs[1];
    data_[2] += rhs[2];
    return *this;
  }
  constexpr Vector3 &operator-=(const Vector3 &rhs) noexcept {
    data_[0] -= rhs[0];
    data_[1] -= rhs[1];
    data_[2] -= rhs[2];
    return *this;
  }

  // dot product (shorthand for a*b)
  constexpr T operator*(const Vector3 &rhs) const noexcept {
    return data_[0] * rhs[0] + data_[1] * rhs[1] + data_[2] * rhs[2];
  }

  constexpr std::array<T, 3> array() const noexcept { return data_; }

protected:
  std::array<T, 3> data_;
};

// -----------------------------------------------------------------------------
//  Matrix3  –  column-major storage
// -----------------------------------------------------------------------------
/**
 * @class Matrix3
 * @brief A lightweight, constexpr, stack-based 3x3 matrix.
 * Uses column-major storage.
 */
template <typename T> class Matrix3 {
public:
  using value_type = T;
  // Constructors
  constexpr Matrix3() noexcept { data_.fill(Vector3<T>{0, 0, 0}); }
  constexpr Matrix3(const Vector3<T> &c0, const Vector3<T> &c1,
                    const Vector3<T> &c2) noexcept
      : data_{c0, c1, c2} {}
  // Column access
  constexpr const Vector3<T> &operator[](std::size_t c) const noexcept {
    return data_[c];
  }
  constexpr Vector3<T> &operator[](std::size_t c) noexcept { return data_[c]; }

  // Element access (row, col)
  constexpr T operator()(std::size_t r, std::size_t c) const noexcept {
    return data_[c][r];
  }
  constexpr T &operator()(std::size_t r, std::size_t c) noexcept {
    return data_[c][r];
  }

  constexpr std::array<std::array<T, 3>, 3> array() const noexcept {
    return {{{data_[0][0], data_[1][0], data_[2][0]},
             {data_[0][1], data_[1][1], data_[2][1]},
             {data_[0][2], data_[1][2], data_[2][2]}}};
  }

  // Trace: sum of diagonal elements
  constexpr T trace() const noexcept {
    return data_[0][0] + data_[1][1] + data_[2][2];
  }

  // Scalar multiplication
  constexpr Matrix3 operator*(T s) const noexcept {
    return Matrix3(data_[0] * s, data_[1] * s, data_[2] * s);
  }

  // In-place scalar multiplication
  constexpr Matrix3 &operator*=(T s) noexcept {
    data_[0] = data_[0] * s;
    data_[1] = data_[1] * s;
    data_[2] = data_[2] * s;
    return *this;
  }

  // In-place matrix addition
  constexpr Matrix3 &operator+=(const Matrix3 &rhs) noexcept {
    data_[0] += rhs.data_[0];
    data_[1] += rhs.data_[1];
    data_[2] += rhs.data_[2];
    return *this;
  }

  // In-place matrix multiplication
  constexpr Matrix3 &operator*=(const Matrix3 &rhs) noexcept {
    const Matrix3 lhs = *this;
    for (std::size_t j = 0; j < 3; ++j) {
      data_[j] = lhs[0] * rhs(0, j) + lhs[1] * rhs(1, j) + lhs[2] * rhs(2, j);
    }
    return *this;
  }

private:
  std::array<Vector3<T>, 3> data_;
};

// -----------------------------------------------------------------------------
//  Specializations for Double (SIMD Optimized)
// -----------------------------------------------------------------------------

#if defined(CORRELATION_SIMD_AVX2) || defined(CORRELATION_SIMD_AVX512)

/**
 * @brief SIMD-optimized specialization of Vector3 for double.
 */
template <> class CORRELATION_ALIGN(32) Vector3<double> {
public:
  using value_type = double;

  Vector3() noexcept { _mm256_store_pd(data_, _mm256_setzero_pd()); }
  Vector3(double x, double y, double z) noexcept {
    CORRELATION_ALIGN(32) double temp[4] = {x, y, z, 0.0};
    _mm256_store_pd(data_, _mm256_load_pd(temp));
  }
  explicit Vector3(const std::array<double, 3> &a) noexcept {
    CORRELATION_ALIGN(32) double temp[4] = {a[0], a[1], a[2], 0.0};
    _mm256_store_pd(data_, _mm256_load_pd(temp));
  }

  double operator[](std::size_t i) const noexcept { return data_[i]; }
  double &operator[](std::size_t i) noexcept { return data_[i]; }
  double x() const noexcept { return data_[0]; }
  double &x() noexcept { return data_[0]; }
  double y() const noexcept { return data_[1]; }
  double &y() noexcept { return data_[1]; }
  double z() const noexcept { return data_[2]; }
  double &z() noexcept { return data_[2]; }
  double &operator()(std::size_t i) noexcept { return data_[i]; }

  bool empty() const noexcept {
    return data_[0] == 0.0 && data_[1] == 0.0 && data_[2] == 0.0;
  }

  // Arithmetic with SIMD
  Vector3 operator+(const Vector3 &rhs) const noexcept {
    Vector3 res;
    _mm256_store_pd(res.data_, _mm256_add_pd(_mm256_load_pd(data_),
                                             _mm256_load_pd(rhs.data_)));
    return res;
  }
  Vector3 operator-(const Vector3 &rhs) const noexcept {
    Vector3 res;
    _mm256_store_pd(res.data_, _mm256_sub_pd(_mm256_load_pd(data_),
                                             _mm256_load_pd(rhs.data_)));
    return res;
  }
  Vector3 operator*(double s) const noexcept {
    Vector3 res;
    _mm256_store_pd(res.data_,
                    _mm256_mul_pd(_mm256_load_pd(data_), _mm256_set1_pd(s)));
    return res;
  }
  Vector3 operator/(double s) const noexcept {
    Vector3 res;
    _mm256_store_pd(res.data_,
                    _mm256_div_pd(_mm256_load_pd(data_), _mm256_set1_pd(s)));
    return res;
  }

  Vector3 &operator+=(const Vector3 &rhs) noexcept {
    _mm256_store_pd(
        data_, _mm256_add_pd(_mm256_load_pd(data_), _mm256_load_pd(rhs.data_)));
    return *this;
  }
  Vector3 &operator-=(const Vector3 &rhs) noexcept {
    _mm256_store_pd(
        data_, _mm256_sub_pd(_mm256_load_pd(data_), _mm256_load_pd(rhs.data_)));
    return *this;
  }

  double operator*(const Vector3 &rhs) const noexcept {
    __m256d mult =
        _mm256_mul_pd(_mm256_load_pd(data_), _mm256_load_pd(rhs.data_));
    // Mask out the 4th element (just in case) or just sum the first 3.
    // Actually, since 4th component is 0, we can sum all 4.
    __m128d lo = _mm256_castpd256_pd128(mult);
    __m128d hi = _mm256_extractf128_pd(mult, 1);
    __m128d res = _mm_add_pd(lo, hi);
    res = _mm_hadd_pd(res, res);
    return _mm_cvtsd_f64(res);
  }

  std::array<double, 3> array() const noexcept {
    return {data_[0], data_[1], data_[2]};
  }

  const double *begin() const noexcept { return data_; }
  const double *end() const noexcept { return data_ + 3; }
  double *begin() noexcept { return data_; }
  double *end() noexcept { return data_ + 3; }

private:
  CORRELATION_ALIGN(32) double data_[4]; // Padded to 4 doubles for AVX
};

#endif // SIMD Specialized Vector3<double>

// -----------------------------------------------------------------------------
//  Free Functions
// -----------------------------------------------------------------------------

/**
 * @brief Computes the dot product of two vectors.
 * 
 * @param a The first vector.
 * @param b The second vector.
 * @return The scalar dot product.
 */
template <typename T>
constexpr T dot(const Vector3<T> &a, const Vector3<T> &b) noexcept {
  return a * b;
}

/**
 * @brief Computes the cross product of two 3D vectors.
 * 
 * @param a The first vector.
 * @param b The second vector.
 * @return A new Vector3 representing the cross product a x b.
 */
template <typename T>
constexpr Vector3<T> cross(const Vector3<T> &a, const Vector3<T> &b) noexcept {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

// SIMD Cross Product for double
#if defined(CORRELATION_SIMD_AVX2) || defined(CORRELATION_SIMD_AVX512)
inline Vector3<double> cross(const Vector3<double> &a,
                             const Vector3<double> &b) noexcept {
  // Cross product (a1b2 - a2b1, a2b0 - a0b2, a0b1 - a1b0)
  // This can be done efficiently with shuffles, but is it worth it for a single
  // cross? Let's use the basic formula for now, but compiler will see the
  // alignment.
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}
#endif

/**
 * @brief Computes the squared vector norm (length).
 * 
 * @param v The vector.
 * @return The squared length of the vector.
 */
template <typename T> constexpr T norm_sq(const Vector3<T> &v) noexcept {
  return v * v;
}
/**
 * @brief Computes the vector norm (length).
 * 
 * @param v The vector.
 * @return The length of the vector.
 */
template <typename T> constexpr T norm(const Vector3<T> &v) noexcept {
  return std::sqrt(v * v);
}

/**
 * @brief Normalizes the given vector.
 * 
 * @param v The vector to normalize.
 * @return A normalized copy of the vector.
 * @throws std::domain_error if the vector length is too close to zero.
 */
template <typename T> inline Vector3<T> normalize(const Vector3<T> &v) {
  const T n = norm(v);
  if (n < static_cast<T>(1e-300))
    throw std::domain_error("normalize: zero-length vector");
  return v / n;
}

/**
 * @brief Computes the determinant of a 3x3 matrix.
 * 
 * @param m The matrix.
 * @return The determinant.
 */
template <typename T> constexpr T determinant(const Matrix3<T> &m) noexcept {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

/**
 * @brief Computes the inverse of a 3x3 matrix.
 * 
 * @param m The matrix to invert.
 * @return The inverted matrix.
 * @throws std::runtime_error if the matrix is singular (determinant near zero).
 */
template <typename T> constexpr Matrix3<T> invert(const Matrix3<T> &m) {
  T det = determinant(m);
  if (std::abs(det) < 1e-15)
    throw std::runtime_error("singular matrix");
  T inv = T{1} / det;

  Vector3<T> c0{m[1][1] * m[2][2] - m[1][2] * m[2][1],
                m[0][2] * m[2][1] - m[0][1] * m[2][2],
                m[0][1] * m[1][2] - m[0][2] * m[1][1]};

  Vector3<T> c1{m[1][2] * m[2][0] - m[1][0] * m[2][2],
                m[0][0] * m[2][2] - m[0][2] * m[2][0],
                m[0][2] * m[1][0] - m[0][0] * m[1][2]};

  Vector3<T> c2{m[1][0] * m[2][1] - m[1][1] * m[2][0],
                m[0][1] * m[2][0] - m[0][0] * m[2][1],
                m[0][0] * m[1][1] - m[0][1] * m[1][0]};

  return Matrix3<T>(c0 * inv, c1 * inv, c2 * inv);
}

/**
 * @brief Computes the transpose of a 3x3 matrix.
 * 
 * @param m The matrix to transpose.
 * @return The transposed matrix.
 */
template <typename T>
constexpr Matrix3<T> transpose(const Matrix3<T> &m) noexcept {
  return Matrix3<T>({m(0, 0), m(1, 0), m(2, 0)}, {m(0, 1), m(1, 1), m(2, 1)},
                    {m(0, 2), m(1, 2), m(2, 2)});
}

// Matrix-Vector product
template <typename T>
constexpr Vector3<T> operator*(const Matrix3<T> &m,
                               const Vector3<T> &v) noexcept {
  Vector3<T> res;
  res += m[0] * v.x();
  res += m[1] * v.y();
  res += m[2] * v.z();
  return res;
}

// Matrix-Matrix product
template <typename T>
constexpr Matrix3<T> operator*(const Matrix3<T> &a,
                               const Matrix3<T> &b) noexcept {
  return Matrix3<T>(a * b[0], a * b[1], a * b[2]);
}

// Scalar * vector product
template <typename Scalar, typename T>
constexpr std::enable_if_t<std::is_arithmetic<Scalar>::value, Vector3<T>>
operator*(Scalar s, const Vector3<T> &v) noexcept {
  return v * static_cast<T>(s);
}

} // namespace correlation::math
