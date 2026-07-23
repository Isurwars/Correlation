/**
 * @file LinearAlgebra.hpp
 * @brief Lightweight 3D vector and matrix classes for geometric calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/Precision.hpp"
#include "math/SIMDConfig.hpp"
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
   * @param x_val X-component value.
   * @param y_val Y-component value.
   * @param z_val Z-component value.
   */
  constexpr Vector3(T x_val, T y_val, T z_val) noexcept : data_{x_val, y_val, z_val} {}

  template <typename U>
  constexpr Vector3(U x_val, U y_val, U z_val) noexcept
      : data_{static_cast<T>(x_val), static_cast<T>(y_val), static_cast<T>(z_val)} {}

  template <typename U>
  constexpr Vector3(const Vector3<U> &other) noexcept
      : data_{static_cast<T>(other.x()), static_cast<T>(other.y()), static_cast<T>(other.z())} {}

  /**
   * @brief Array-based constructor.
   * @param arr Standard array containing {x, y, z}.
   */
  constexpr explicit Vector3(const std::array<T, 3> &arr) noexcept : data_{arr} {}

  /**
   * @brief Access component by index (0, 1, or 2).
   * @param idx Index of the component.
   * @return The value at index i.
   */
  [[nodiscard]] constexpr T operator[](std::size_t idx) const noexcept { return data_.at(idx); }

  /**
   * @brief Access mutable component by index (0, 1, or 2).
   * @param idx Index of the component.
   * @return Reference to the value at index i.
   */
  constexpr T &operator[](std::size_t idx) noexcept { return data_.at(idx); }

  /** @return The X-component. */
  [[nodiscard]] constexpr T x() const noexcept { return data_[0]; }
  /** @return Reference to the X-component. */
  constexpr T &x() noexcept { return data_[0]; }
  /** @return The Y-component. */
  [[nodiscard]] constexpr T y() const noexcept { return data_[1]; }
  /** @return Reference to the Y-component. */
  constexpr T &y() noexcept { return data_[1]; }
  /** @return The Z-component. */
  [[nodiscard]] constexpr T z() const noexcept { return data_[2]; }
  /** @return Reference to the Z-component. */
  constexpr T &z() noexcept { return data_[2]; }

  /**
   * @brief Access component by index (alternative syntax).
   * @param idx Index.
   * @return Reference to component.
   */
  constexpr T &operator()(std::size_t idx) noexcept { return data_.at(idx); }

  /**
   * @brief Access component by index (alternative syntax).
   * @param idx Index.
   * @return Component value.
   */
  [[nodiscard]] constexpr T operator()(std::size_t idx) const noexcept { return data_.at(idx); }

  /**
   * @brief Checks if the vector is a zero vector.
   * @return True if every component is exactly zero.
   */
  [[nodiscard]] constexpr bool empty() const noexcept { return data_[0] == T{} && data_[1] == T{} && data_[2] == T{}; }

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
  [[nodiscard]] constexpr Vector3 operator+(const Vector3 &rhs) const noexcept {
    return {data_[0] + rhs[0], data_[1] + rhs[1], data_[2] + rhs[2]};
  }

  /**
   * @brief Vector subtraction.
   * @param rhs The vector to subtract.
   * @return A new vector representing the difference.
   */
  [[nodiscard]] constexpr Vector3 operator-(const Vector3 &rhs) const noexcept {
    return {data_[0] - rhs[0], data_[1] - rhs[1], data_[2] - rhs[2]};
  }

  /**
   * @brief Unary negation.
   * @return A new vector with all components negated.
   */
  [[nodiscard]] constexpr Vector3 operator-() const noexcept { return {-data_[0], -data_[1], -data_[2]}; }

  /**
   * @brief Scalar multiplication.
   * @param scalar The scalar factor.
   * @return A new vector scaled by s.
   */
  [[nodiscard]] constexpr Vector3 operator*(T scalar) const noexcept {
    return {data_[0] * scalar, data_[1] * scalar, data_[2] * scalar};
  }

  /**
   * @brief Scalar division.
   * @param scalar The scalar divisor.
   * @return A new vector scaled by 1/s.
   */
  [[nodiscard]] constexpr Vector3 operator/(T scalar) const noexcept {
    return {data_[0] / scalar, data_[1] / scalar, data_[2] / scalar};
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
   * @brief In-place scalar multiplication.
   * @param scalar The scalar factor.
   * @return Reference to this vector.
   */
  constexpr Vector3 &operator*=(T scalar) noexcept {
    data_[0] *= scalar;
    data_[1] *= scalar;
    data_[2] *= scalar;
    return *this;
  }

  /**
   * @brief In-place scalar division.
   * @param scalar The scalar divisor.
   * @return Reference to this vector.
   */
  constexpr Vector3 &operator/=(T scalar) noexcept {
    data_[0] /= scalar;
    data_[1] /= scalar;
    data_[2] /= scalar;
    return *this;
  }

  /**
   * @brief Vector-vector multiplication (dot product shorthand).
   * @see dot()
   * @param rhs The vector to multiply with.
   * @return The scalar result of the dot product.
   */
  [[nodiscard]] constexpr T operator*(const Vector3 &rhs) const noexcept {
    return data_[0] * rhs[0] + data_[1] * rhs[1] + data_[2] * rhs[2];
  }

  /**
   * @brief Converts the vector to a std::array.
   * @return Array containing {x, y, z}.
   */
  [[nodiscard]] constexpr std::array<T, 3> array() const noexcept { return data_; }

  /**
   * @brief Equality comparison.
   * @param rhs The vector to compare with.
   * @return True if all components are exactly equal.
   */
  [[nodiscard]] constexpr bool operator==(const Vector3 &rhs) const noexcept {
    return data_[0] == rhs[0] && data_[1] == rhs[1] && data_[2] == rhs[2];
  }

  /**
   * @brief Inequality comparison.
   * @param rhs The vector to compare with.
   * @return True if any component differs.
   */
  [[nodiscard]] constexpr bool operator!=(const Vector3 &rhs) const noexcept { return !(*this == rhs); }

private:
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
  /** @brief Default constructor. Initializes all elements to zero. */
  constexpr Matrix3() noexcept { data_.fill(Vector3<T>{0, 0, 0}); }

  /**
   * @brief Returns the 3×3 identity matrix.
   * @return Identity matrix with ones on the diagonal.
   */
  [[nodiscard]] static constexpr Matrix3 identity() noexcept {
    return Matrix3({T{1}, T{0}, T{0}}, {T{0}, T{1}, T{0}}, {T{0}, T{0}, T{1}});
  }

  /**
   * @brief Parameterized constructor from column vectors.
   * @param c_zero First column vector.
   * @param c_one Second column vector.
   * @param c_two Third column vector.
   */
  constexpr Matrix3(const Vector3<T> &c_zero, const Vector3<T> &c_one, const Vector3<T> &c_two) noexcept
      : data_{c_zero, c_one, c_two} {}

  template <typename U>
  constexpr Matrix3(const Matrix3<U> &other) noexcept
      : data_{Vector3<T>(other[0]), Vector3<T>(other[1]), Vector3<T>(other[2])} {}
  /**
   * @brief Access column vector by index.
   * @param col_idx Column index (0, 1, or 2).
   * @return Constant reference to the column vector.
   */
  [[nodiscard]] constexpr const Vector3<T> &operator[](std::size_t col_idx) const noexcept { return data_.at(col_idx); }

  /**
   * @brief Access mutable column vector by index.
   * @param col_idx Column index (0, 1, or 2).
   * @return Reference to the column vector.
   */
  constexpr Vector3<T> &operator[](std::size_t col_idx) noexcept { return data_.at(col_idx); }

  /**
   * @brief Access element by row and column.
   * @param row_idx Row index.
   * @param col_idx Column index.
   * @return The element value.
   */
  [[nodiscard]] constexpr T operator()(std::size_t row_idx, std::size_t col_idx) const noexcept {
    return data_.at(col_idx)[row_idx];
  }

  /**
   * @brief Access mutable element by row and column.
   * @param row_idx Row index.
   * @param col_idx Column index.
   * @return Reference to the element.
   */
  constexpr T &operator()(std::size_t row_idx, std::size_t col_idx) noexcept { return data_.at(col_idx)[row_idx]; }

  /**
   * @brief Converts the matrix to a nested std::array.
   * @return Row-major nested array representation.
   */
  [[nodiscard]] constexpr std::array<std::array<T, 3>, 3> array() const noexcept {
    return {{{data_[0][0], data_[1][0], data_[2][0]},
             {data_[0][1], data_[1][1], data_[2][1]},
             {data_[0][2], data_[1][2], data_[2][2]}}};
  }

  /**
   * @brief Computes the matrix trace.
   * @return Sum of the diagonal elements.
   */
  [[nodiscard]] constexpr T trace() const noexcept { return data_[0][0] + data_[1][1] + data_[2][2]; }

  /**
   * @brief Matrix-scalar multiplication.
   * @param scalar Scalar factor.
   * @return A scaled matrix object.
   */
  [[nodiscard]] constexpr Matrix3 operator*(T scalar) const noexcept {
    return Matrix3(data_[0] * scalar, data_[1] * scalar, data_[2] * scalar);
  }

  /**
   * @brief In-place scalar multiplication.
   * @param scalar Scalar factor.
   * @return Reference to this matrix.
   */
  constexpr Matrix3 &operator*=(T scalar) noexcept {
    data_[0] *= scalar;
    data_[1] *= scalar;
    data_[2] *= scalar;
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
   * @brief Matrix addition.
   * @param rhs The matrix to add.
   * @return A new matrix representing the sum.
   */
  [[nodiscard]] constexpr Matrix3 operator+(const Matrix3 &rhs) const noexcept {
    return Matrix3(data_[0] + rhs.data_[0], data_[1] + rhs.data_[1], data_[2] + rhs.data_[2]);
  }

  /**
   * @brief In-place matrix subtraction.
   * @param rhs The matrix to subtract.
   * @return Reference to this matrix.
   */
  constexpr Matrix3 &operator-=(const Matrix3 &rhs) noexcept {
    data_[0] -= rhs.data_[0];
    data_[1] -= rhs.data_[1];
    data_[2] -= rhs.data_[2];
    return *this;
  }

  /**
   * @brief Matrix subtraction.
   * @param rhs The matrix to subtract.
   * @return A new matrix representing the difference.
   */
  [[nodiscard]] constexpr Matrix3 operator-(const Matrix3 &rhs) const noexcept {
    return Matrix3(data_[0] - rhs.data_[0], data_[1] - rhs.data_[1], data_[2] - rhs.data_[2]);
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

  /**
   * @brief Equality comparison.
   * @param rhs The matrix to compare with.
   * @return True if all elements are exactly equal.
   */
  [[nodiscard]] constexpr bool operator==(const Matrix3 &rhs) const noexcept {
    return data_[0] == rhs.data_[0] && data_[1] == rhs.data_[1] && data_[2] == rhs.data_[2];
  }

  /**
   * @brief Inequality comparison.
   * @param rhs The matrix to compare with.
   * @return True if any element differs.
   */
  [[nodiscard]] constexpr bool operator!=(const Matrix3 &rhs) const noexcept { return !(*this == rhs); }

private:
  std::array<Vector3<T>, 3> data_; ///< Column vectors of the matrix.
};

// -----------------------------------------------------------------------------
//  Specializations for Double and Float (SIMD Optimized)
// -----------------------------------------------------------------------------

#include "math/detail/Vector3SIMD.hpp" // IWYU pragma: export

// -----------------------------------------------------------------------------
//  Free Functions
// -----------------------------------------------------------------------------

/**
 * @brief Computes the dot product of two vectors.
 *
 * @param vec_a The first vector.
 * @param vec_b The second vector.
 * @return The scalar dot product.
 */
template <typename T> [[nodiscard]] constexpr T dot(const Vector3<T> &vec_a, const Vector3<T> &vec_b) noexcept {
  return vec_a * vec_b;
}

/**
 * @brief Computes the cross product of two 3D vectors.
 *
 * @param vec_a The first vector.
 * @param vec_b The second vector.
 * @return A new Vector3 representing the cross product vec_a x vec_b.
 */
template <typename T>
[[nodiscard]] constexpr Vector3<T> cross(const Vector3<T> &vec_a, const Vector3<T> &vec_b) noexcept {
  return {vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1], vec_a[2] * vec_b[0] - vec_a[0] * vec_b[2],
          vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]};
}

/**
 * @brief Computes the cross product of two 3D vectors.
 *
 * Optimized specialization for double using manual calculations to help the
 * compiler generate aligned SIMD code.
 *
 * @param vec_a The first vector.
 * @param vec_b The second vector.
 * @return A new Vector3 representing the cross product vec_a x vec_b.
 */
#if defined(CORRELATION_SIMD_AVX2) || defined(CORRELATION_SIMD_AVX512)
[[nodiscard]] inline Vector3<double> cross(const Vector3<double> &vec_a, const Vector3<double> &vec_b) noexcept {
  return {vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1], vec_a[2] * vec_b[0] - vec_a[0] * vec_b[2],
          vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]};
}
#endif

/**
 * @brief Computes the squared vector norm (length).
 *
 * @param vec_a The vector.
 * @return The squared length of the vector.
 */
template <typename T> [[nodiscard]] constexpr T norm_sq(const Vector3<T> &vec_a) noexcept { return vec_a * vec_a; }
/**
 * @brief Computes the vector norm (length).
 *
 * @param vec_a The vector.
 * @return The length of the vector.
 */
template <typename T> [[nodiscard]] constexpr T norm(const Vector3<T> &vec_a) noexcept {
  return std::sqrt(vec_a * vec_a);
}

/**
 * @brief Computes the distance between two vectors.
 *
 * @param vec_a The first vector.
 * @param vec_b The second vector.
 * @return The Euclidean distance between vec_a and vec_b.
 */
template <typename T> [[nodiscard]] constexpr T distance(const Vector3<T> &vec_a, const Vector3<T> &vec_b) noexcept {
  return norm(vec_a - vec_b);
}

/**
 * @brief Normalizes the given vector.
 *
 * @param vec_a The vector to normalize.
 * @return A normalized copy of the vector.
 * @throws std::domain_error if the vector length is too close to zero.
 */
template <typename T> [[nodiscard]] constexpr Vector3<T> normalize(const Vector3<T> &vec_a) {
  const T norm_val = norm(vec_a);
  if (norm_val < static_cast<T>(1e-300)) {
    throw std::domain_error("normalize: zero-length vector");
  }
  return vec_a / norm_val;
}

/**
 * @brief Computes the determinant of a 3x3 matrix.
 *
 * @param matrix The matrix.
 * @return The determinant.
 */
template <typename T> [[nodiscard]] constexpr T determinant(const Matrix3<T> &matrix) noexcept {
  return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
         matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
         matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

/**
 * @brief Computes the inverse of a 3x3 matrix.
 *
 * @param matrix The matrix to invert.
 * @return The inverted matrix.
 * @throws std::runtime_error if the matrix is singular (determinant near zero).
 */
template <typename T> [[nodiscard]] constexpr Matrix3<T> invert(const Matrix3<T> &matrix) {
  const T det = determinant(matrix);
  if (std::abs(det) < 1e-15) {
    throw std::runtime_error("singular matrix");
  }
  const T inv = T{1} / det;

  const Vector3<T> c_zero{matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1],
                          matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2],
                          matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]};

  const Vector3<T> c_one{matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2],
                         matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0],
                         matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]};

  const Vector3<T> c_two{matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0],
                         matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1],
                         matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]};

  return Matrix3<T>(c_zero * inv, c_one * inv, c_two * inv);
}

/**
 * @brief Computes the transpose of a 3x3 matrix.
 *
 * @param matrix The matrix to transpose.
 * @return The transposed matrix.
 */
template <typename T> [[nodiscard]] constexpr Matrix3<T> transpose(const Matrix3<T> &matrix) noexcept {
  return Matrix3<T>({matrix(0, 0), matrix(0, 1), matrix(0, 2)}, {matrix(1, 0), matrix(1, 1), matrix(1, 2)},
                    {matrix(2, 0), matrix(2, 1), matrix(2, 2)});
}

/**
 * @brief Matrix-vector multiplication.
 * @tparam T Coordinate type.
 * @param matrix 3x3 matrix.
 * @param vector 3D vector.
 * @return Transformed vector matrix * vector.
 */
template <typename T, typename U>
[[nodiscard]] constexpr auto operator*(const Matrix3<T> &matrix, const Vector3<U> &vector) noexcept {
  using CommonT = std::common_type_t<T, U>;
  return Vector3<CommonT>(std::fma(static_cast<CommonT>(vector.z()), static_cast<CommonT>(matrix[2].x()),
                                   std::fma(static_cast<CommonT>(vector.y()), static_cast<CommonT>(matrix[1].x()),
                                            static_cast<CommonT>(vector.x()) * static_cast<CommonT>(matrix[0].x()))),
                          std::fma(static_cast<CommonT>(vector.z()), static_cast<CommonT>(matrix[2].y()),
                                   std::fma(static_cast<CommonT>(vector.y()), static_cast<CommonT>(matrix[1].y()),
                                            static_cast<CommonT>(vector.x()) * static_cast<CommonT>(matrix[0].y()))),
                          std::fma(static_cast<CommonT>(vector.z()), static_cast<CommonT>(matrix[2].z()),
                                   std::fma(static_cast<CommonT>(vector.y()), static_cast<CommonT>(matrix[1].z()),
                                            static_cast<CommonT>(vector.x()) * static_cast<CommonT>(matrix[0].z()))));
}

/**
 * @brief Matrix-matrix multiplication.
 * @tparam T Scalar type.
 * @param mat_a Left matrix.
 * @param mat_b Right matrix.
 * @return Resulting matrix product mat_a * mat_b.
 */
template <typename T>
[[nodiscard]] constexpr Matrix3<T> operator*(const Matrix3<T> &mat_a, const Matrix3<T> &mat_b) noexcept {
  return Matrix3<T>(mat_a * mat_b[0], mat_a * mat_b[1], mat_a * mat_b[2]);
}

/**
 * @brief Scalar-vector multiplication.
 * @tparam Scalar Arithmetic type.
 * @tparam T Vector coordinate type.
 * @param scalar Scaling factor.
 * @param vec Input vector.
 * @return Scaled vector scalar * vec.
 */
template <typename Scalar, typename T>
[[nodiscard]] constexpr Vector3<T> operator*(Scalar scalar, const Vector3<T> &vec) noexcept
  requires std::is_arithmetic_v<Scalar>
{
  return vec * static_cast<T>(scalar);
}

/// Precision-configurable aliases for 3D vectors and 3x3 matrices.
using Vector3R = Vector3<real_t>;
using Matrix3R = Matrix3<real_t>;

} // namespace correlation::math
