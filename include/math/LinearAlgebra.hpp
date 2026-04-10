/**
 * @file LinearAlgebra.hpp
 * @brief Lightweight 3D vector and matrix classes for geometric calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
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
 * 
 * Provides basic linear algebra operations for 3D geometric calculations.
 * Supports SIMD optimization for specific types (e.g., double).
 * 
 * @tparam T The scalar type stored by the vector (e.g., float, double, int).
 */
template <typename T> class Vector3 {
public:
  using value_type = T; ///< Scalar type.

  /** @brief Default constructor. Initializes all components to zero. */
  constexpr Vector3() noexcept : data_{0, 0, 0} {}

  /**
   * @brief Parameterized constructor.
   * @param x X-component value.
   * @param y Y-component value.
   * @param z Z-component value.
   */
  constexpr Vector3(T x, T y, T z) noexcept : data_{x, y, z} {}

  /**
   * @brief Array-based constructor.
   * @param a Standard array containing {x, y, z}.
   */
  constexpr explicit Vector3(const std::array<T, 3> &a) noexcept : data_{a} {}

  /**
   * @brief Access component by index (0, 1, or 2).
   * @param i Index of the component.
   * @return The value at index i.
   */
  constexpr T operator[](std::size_t i) const noexcept { return data_[i]; }

  /**
   * @brief Access mutable component by index (0, 1, or 2).
   * @param i Index of the component.
   * @return Reference to the value at index i.
   */
  constexpr T &operator[](std::size_t i) noexcept { return data_[i]; }

  /** @return The X-component. */
  constexpr T x() const noexcept { return data_[0]; }
  /** @return Reference to the X-component. */
  constexpr T &x() noexcept { return data_[0]; }
  /** @return The Y-component. */
  constexpr T y() const noexcept { return data_[1]; }
  /** @return Reference to the Y-component. */
  constexpr T &y() noexcept { return data_[1]; }
  /** @return The Z-component. */
  constexpr T z() const noexcept { return data_[2]; }
  /** @return Reference to the Z-component. */
  constexpr T &z() noexcept { return data_[2]; }

  /**
   * @brief Access component by index (alternative syntax).
   * @param i Index.
   * @return Reference to component.
   */
  constexpr T &operator()(std::size_t i) noexcept { return data_[i]; }

  /**
   * @brief Access component by index (alternative syntax).
   * @param i Index.
   * @return Component value.
   */
  constexpr T operator()(std::size_t i) const noexcept { return data_[i]; }

  /**
   * @brief Checks if the vector is a zero vector.
   * @return True if every component is exactly zero.
   */
  constexpr bool empty() const noexcept {
    return data_[0] == T{} && data_[1] == T{} && data_[2] == T{};
  }

  /** @return Pointer to the beginning of the underlying data. */
  constexpr const T *begin() const noexcept { return data_.data(); }
  /** @return Pointer to the end of the underlying data. */
  constexpr const T *end() const noexcept { return data_.data() + 3; }
  /** @return Mutable pointer to the beginning of the underlying data. */
  constexpr T *begin() noexcept { return data_.data(); }
  /** @return Mutable pointer to the end of the underlying data. */
  constexpr T *end() noexcept { return data_.data() + 3; }

  /**
   * @brief Vector addition.
   * @param rhs The vector to add.
   * @return A new vector representing the sum.
   */
  constexpr Vector3 operator+(const Vector3 &rhs) const noexcept {
    return {data_[0] + rhs[0], data_[1] + rhs[1], data_[2] + rhs[2]};
  }

  /**
   * @brief Vector subtraction.
   * @param rhs The vector to subtract.
   * @return A new vector representing the difference.
   */
  constexpr Vector3 operator-(const Vector3 &rhs) const noexcept {
    return {data_[0] - rhs[0], data_[1] - rhs[1], data_[2] - rhs[2]};
  }

  /**
   * @brief Scalar multiplication.
   * @param s The scalar factor.
   * @return A new vector scaled by s.
   */
  constexpr Vector3 operator*(T s) const noexcept {
    return {data_[0] * s, data_[1] * s, data_[2] * s};
  }

  /**
   * @brief Scalar division.
   * @param s The scalar divisor.
   * @return A new vector scaled by 1/s.
   */
  constexpr Vector3 operator/(T s) const noexcept {
    return {data_[0] / s, data_[1] / s, data_[2] / s};
  }

  /**
   * @brief In-place addition.
   * @param rhs The vector to add.
   * @return Reference to this vector.
   */
  constexpr Vector3 &operator+=(const Vector3 &rhs) noexcept {
    data_[0] += rhs[0];
    data_[1] += rhs[1];
    data_[2] += rhs[2];
    return *this;
  }

  /**
   * @brief In-place subtraction.
   * @param rhs The vector to subtract.
   * @return Reference to this vector.
   */
  constexpr Vector3 &operator-=(const Vector3 &rhs) noexcept {
    data_[0] -= rhs[0];
    data_[1] -= rhs[1];
    data_[2] -= rhs[2];
    return *this;
  }

  /**
   * @brief Vector-vector multiplication (dot product shorthand).
   * @see dot()
   * @param rhs The vector to multiply with.
   * @return The scalar result of the dot product.
   */
  constexpr T operator*(const Vector3 &rhs) const noexcept {
    return data_[0] * rhs[0] + data_[1] * rhs[1] + data_[2] * rhs[2];
  }

  /**
   * @brief Converts the vector to a std::array.
   * @return Array containing {x, y, z}.
   */
  constexpr std::array<T, 3> array() const noexcept { return data_; }

protected:
  std::array<T, 3> data_; ///< Underlying vector components.
};

// -----------------------------------------------------------------------------
//  Matrix3  –  column-major storage
// -----------------------------------------------------------------------------
/**
 * @class Matrix3
 * @brief A lightweight, constexpr, stack-based 3x3 matrix.
 * 
 * Uses column-major storage format (standard for OpenGL and BLAS).
 * This means data is stored as an array of three column vectors.
 * 
 * @tparam T The scalar type stored by the matrix.
 */
template <typename T> class Matrix3 {
public:
  using value_type = T; ///< Scalar type.
  /** @brief Default constructor. Initializes an identity-like zero matrix. */
  constexpr Matrix3() noexcept { data_.fill(Vector3<T>{0, 0, 0}); }

  /**
   * @brief Parameterized constructor from column vectors.
   * @param c0 First column vector.
   * @param c1 Second column vector.
   * @param c2 Third column vector.
   */
  constexpr Matrix3(const Vector3<T> &c0, const Vector3<T> &c1,
                    const Vector3<T> &c2) noexcept
      : data_{c0, c1, c2} {}
  /**
   * @brief Access column vector by index.
   * @param c Column index (0, 1, or 2).
   * @return Constant reference to the column vector.
   */
  constexpr const Vector3<T> &operator[](std::size_t c) const noexcept {
    return data_[c];
  }

  /**
   * @brief Access mutable column vector by index.
   * @param c Column index (0, 1, or 2).
   * @return Reference to the column vector.
   */
  constexpr Vector3<T> &operator[](std::size_t c) noexcept { return data_[c]; }

  /**
   * @brief Access element by row and column.
   * @param r Row index.
   * @param c Column index.
   * @return The element value.
   */
  constexpr T operator()(std::size_t r, std::size_t c) const noexcept {
    return data_[c][r];
  }

  /**
   * @brief Access mutable element by row and column.
   * @param r Row index.
   * @param c Column index.
   * @return Reference to the element.
   */
  constexpr T &operator()(std::size_t r, std::size_t c) noexcept {
    return data_[c][r];
  }

  /**
   * @brief Converts the matrix to a nested std::array.
   * @return Row-major nested array representation.
   */
  constexpr std::array<std::array<T, 3>, 3> array() const noexcept {
    return {{{data_[0][0], data_[1][0], data_[2][0]},
             {data_[0][1], data_[1][1], data_[2][1]},
             {data_[0][2], data_[1][2], data_[2][2]}}};
  }

  /**
   * @brief Computes the matrix trace.
   * @return Sum of the diagonal elements.
   */
  constexpr T trace() const noexcept {
    return data_[0][0] + data_[1][1] + data_[2][2];
  }

  /**
   * @brief Matrix-scalar multiplication.
   * @param s Scalar factor.
   * @return A scaled matrix object.
   */
  constexpr Matrix3 operator*(T s) const noexcept {
    return Matrix3(data_[0] * s, data_[1] * s, data_[2] * s);
  }

  /**
   * @brief In-place scalar multiplication.
   * @param s Scalar factor.
   * @return Reference to this matrix.
   */
  constexpr Matrix3 &operator*=(T s) noexcept {
    data_[0] = data_[0] * s;
    data_[1] = data_[1] * s;
    data_[2] = data_[2] * s;
    return *this;
  }

  /**
   * @brief In-place matrix addition.
   * @param rhs The matrix to add.
   * @return Reference to this matrix.
   */
  constexpr Matrix3 &operator+=(const Matrix3 &rhs) noexcept {
    data_[0] += rhs.data_[0];
    data_[1] += rhs.data_[1];
    data_[2] += rhs.data_[2];
    return *this;
  }

  /**
   * @brief In-place matrix-matrix multiplication.
   * @param rhs The matrix to multiply with (post-multiplication).
   * @return Reference to this matrix.
   */
  constexpr Matrix3 &operator*=(const Matrix3 &rhs) noexcept {
    const Matrix3 lhs = *this;
    for (std::size_t j = 0; j < 3; ++j) {
      data_[j] = lhs[0] * rhs(0, j) + lhs[1] * rhs(1, j) + lhs[2] * rhs(2, j);
    }
    return *this;
  }

private:
  std::array<Vector3<T>, 3> data_; ///< Column vectors of the matrix.
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
  using value_type = double; ///< Double scalar type.

  /** @brief Default constructor. Initializes to zero and ensures 32-byte alignment. */
  Vector3() noexcept { _mm256_store_pd(data_, _mm256_setzero_pd()); }

  /**
   * @brief Parameterized constructor.
   * @param x X-component.
   * @param y Y-component.
   * @param z Z-component.
   */
  Vector3(double x, double y, double z) noexcept {
    CORRELATION_ALIGN(32) double temp[4] = {x, y, z, 0.0};
    _mm256_store_pd(data_, _mm256_load_pd(temp));
  }

  /**
   * @brief Array-based constructor.
   * @param a Input array {x, y, z}.
   */
  explicit Vector3(const std::array<double, 3> &a) noexcept {
    CORRELATION_ALIGN(32) double temp[4] = {a[0], a[1], a[2], 0.0};
    _mm256_store_pd(data_, _mm256_load_pd(temp));
  }

  /**
   * @brief Access component by index.
   * @param i Index (0, 1, or 2).
   * @return Component value.
   */
  double operator[](std::size_t i) const noexcept { return data_[i]; }

  /**
   * @brief Access mutable component by index.
   * @param i Index (0, 1, or 2).
   * @return Reference to component.
   */
  double &operator[](std::size_t i) noexcept { return data_[i]; }

  /** @return X-component. */
  double x() const noexcept { return data_[0]; }
  /** @return Reference to X-component. */
  double &x() noexcept { return data_[0]; }
  /** @return Y-component. */
  double y() const noexcept { return data_[1]; }
  /** @return Reference to Y-component. */
  double &y() noexcept { return data_[1]; }
  /** @return Z-component. */
  double z() const noexcept { return data_[2]; }
  /** @return Reference to Z-component. */
  double &z() noexcept { return data_[2]; }

  /**
   * @brief Access component by index.
   * @param i Index (0, 1, or 2).
   * @return Reference to component.
   */
  double &operator()(std::size_t i) noexcept { return data_[i]; }

  /** @return True if x, y, and z are all 0.0. */
  bool empty() const noexcept {
    return data_[0] == 0.0 && data_[1] == 0.0 && data_[2] == 0.0;
  }

  // Arithmetic with SIMD
  /**
   * @brief SIMD-accelerated vector addition.
   * @param rhs Vector to add.
   * @return Sum vector.
   */
  Vector3 operator+(const Vector3 &rhs) const noexcept {
    Vector3 res;
    _mm256_store_pd(res.data_, _mm256_add_pd(_mm256_load_pd(data_),
                                             _mm256_load_pd(rhs.data_)));
    return res;
  }

  /**
   * @brief SIMD-accelerated vector subtraction.
   * @param rhs Vector to subtract.
   * @return Difference vector.
   */
  Vector3 operator-(const Vector3 &rhs) const noexcept {
    Vector3 res;
    _mm256_store_pd(res.data_, _mm256_sub_pd(_mm256_load_pd(data_),
                                             _mm256_load_pd(rhs.data_)));
    return res;
  }

  /**
   * @brief SIMD-accelerated scalar multiplication.
   * @param s Scalar factor.
   * @return Scaled vector.
   */
  Vector3 operator*(double s) const noexcept {
    Vector3 res;
    _mm256_store_pd(res.data_,
                    _mm256_mul_pd(_mm256_load_pd(data_), _mm256_set1_pd(s)));
    return res;
  }

  /**
   * @brief SIMD-accelerated scalar division.
   * @param s Scalar divisor.
   * @return Vector scaled by 1/s.
   */
  Vector3 operator/(double s) const noexcept {
    Vector3 res;
    _mm256_store_pd(res.data_,
                    _mm256_div_pd(_mm256_load_pd(data_), _mm256_set1_pd(s)));
    return res;
  }

  /**
   * @brief SIMD-accelerated in-place addition.
   * @param rhs Vector to add.
   * @return Reference to this vector.
   */
  Vector3 &operator+=(const Vector3 &rhs) noexcept {
    _mm256_store_pd(
        data_, _mm256_add_pd(_mm256_load_pd(data_), _mm256_load_pd(rhs.data_)));
    return *this;
  }

  /**
   * @brief SIMD-accelerated in-place subtraction.
   * @param rhs Vector to subtract.
   * @return Reference to this vector.
   */
  Vector3 &operator-=(const Vector3 &rhs) noexcept {
    _mm256_store_pd(
        data_, _mm256_sub_pd(_mm256_load_pd(data_), _mm256_load_pd(rhs.data_)));
    return *this;
  }

  /**
   * @brief SIMD-accelerated dot product.
   * @param rhs Vector to multiply with.
   * @return Scalar result.
   */
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

  /**
   * @brief Converts the vector to a std::array.
   * @return Array containing {x, y, z}.
   */
  std::array<double, 3> array() const noexcept {
    return {data_[0], data_[1], data_[2]};
  }

  /** @return Constant pointer to the beginning of the SIMD-aligned data. */
  const double *begin() const noexcept { return data_; }
  /** @return Constant pointer to the end of the data. */
  const double *end() const noexcept { return data_ + 3; }
  /** @return Mutable pointer to the beginning of the SIMD-aligned data. */
  double *begin() noexcept { return data_; }
  /** @return Mutable pointer to the end of the data. */
  double *end() noexcept { return data_ + 3; }

private:
  CORRELATION_ALIGN(32) double data_[4]; ///< Padded to 4 doubles for 256-bit SIMD alignment.
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

/**
 * @brief Computes the cross product of two 3D vectors.
 * 
 * Optimized specialization for double using manual calculations to help the
 * compiler generate aligned SIMD code.
 *
 * @param a The first vector.
 * @param b The second vector.
 * @return A new Vector3 representing the cross product a x b.
 */
#if defined(CORRELATION_SIMD_AVX2) || defined(CORRELATION_SIMD_AVX512)
inline Vector3<double> cross(const Vector3<double> &a,
                             const Vector3<double> &b) noexcept {
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

/**
 * @brief Matrix-vector multiplication.
 * @tparam T Coordinate type.
 * @param m 3x3 matrix.
 * @param v 3D vector.
 * @return Transformed vector m*v.
 */
template <typename T>
constexpr Vector3<T> operator*(const Matrix3<T> &m,
                               const Vector3<T> &v) noexcept {
  Vector3<T> res;
  res += m[0] * v.x();
  res += m[1] * v.y();
  res += m[2] * v.z();
  return res;
}

/**
 * @brief Matrix-matrix multiplication.
 * @tparam T Scalar type.
 * @param a Left matrix.
 * @param b Right matrix.
 * @return Resulting matrix product a*b.
 */
template <typename T>
constexpr Matrix3<T> operator*(const Matrix3<T> &a,
                               const Matrix3<T> &b) noexcept {
  return Matrix3<T>(a * b[0], a * b[1], a * b[2]);
}

/**
 * @brief Scalar-vector multiplication.
 * @tparam Scalar Arithmetic type.
 * @tparam T Vector coordinate type.
 * @param s Scaling factor.
 * @param v Input vector.
 * @return Scaled vector s*v.
 */
template <typename Scalar, typename T>
constexpr std::enable_if_t<std::is_arithmetic<Scalar>::value, Vector3<T>>
operator*(Scalar s, const Vector3<T> &v) noexcept {
  return v * static_cast<T>(s);
}

} // namespace correlation::math
