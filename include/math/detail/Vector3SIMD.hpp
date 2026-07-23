/**
 * @file Vector3SIMD.hpp
 * @brief Explicit SIMD specializations for Vector3<double> and Vector3<float>.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/Precision.hpp" // IWYU pragma: export
#include "math/SIMDConfig.hpp"

#include <array>
#include <cstddef>

namespace correlation::math {

template <typename T> class Vector3;

#if defined(CORRELATION_SIMD_AVX2) || defined(CORRELATION_SIMD_AVX512)

/**
 * @brief SIMD-optimized specialization of Vector3 for double.
 */
template <> class CORRELATION_ALIGN(32) Vector3<double> {
public:
  using value_type = double; ///< Double scalar type.

  /** @brief Default constructor. Initializes to zero and ensures 32-byte alignment. */
  Vector3() noexcept { _mm256_store_pd(data_.data(), _mm256_setzero_pd()); }

  /**
   * @brief Parameterized constructor.
   * @param x_val X-component.
   * @param y_val Y-component.
   * @param z_val Z-component.
   */
  Vector3(double x_val, double y_val, double z_val) noexcept {
    CORRELATION_ALIGN(32) std::array<double, 4> temp = {x_val, y_val, z_val, 0.0};
    _mm256_store_pd(data_.data(), _mm256_load_pd(temp.data()));
  }

  /**
   * @brief Array-based constructor.
   * @param arr Input array {x, y, z}.
   */
  explicit Vector3(const std::array<double, 3> &arr) noexcept {
    CORRELATION_ALIGN(32) std::array<double, 4> temp = {arr[0], arr[1], arr[2], 0.0};
    _mm256_store_pd(data_.data(), _mm256_load_pd(temp.data()));
  }

  /**
   * @brief Access component by index.
   * @param idx Index (0, 1, or 2).
   * @return Component value.
   */
  [[nodiscard]] double operator[](std::size_t idx) const noexcept {
    return data_[idx]; // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
  }

  /**
   * @brief Access mutable component by index.
   * @param idx Index (0, 1, or 2).
   * @return Reference to component.
   */
  double &operator[](std::size_t idx) noexcept {
    return data_[idx]; // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
  }

  /** @return X-component. */
  [[nodiscard]] double x() const noexcept { return data_[0]; }
  /** @return Reference to X-component. */
  double &x() noexcept { return data_[0]; }
  /** @return Y-component. */
  [[nodiscard]] double y() const noexcept { return data_[1]; }
  /** @return Reference to Y-component. */
  double &y() noexcept { return data_[1]; }
  /** @return Z-component. */
  [[nodiscard]] double z() const noexcept { return data_[2]; }
  /** @return Reference to Z-component. */
  double &z() noexcept { return data_[2]; }

  /**
   * @brief Access component by index.
   * @param idx Index (0, 1, or 2).
   * @return Reference to component.
   */
  double &operator()(std::size_t idx) noexcept {
    return data_[idx]; // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
  }

  /** @return True if x, y, and z are all 0.0. */
  [[nodiscard]] bool empty() const noexcept { return data_[0] == 0.0 && data_[1] == 0.0 && data_[2] == 0.0; }

  // Arithmetic with SIMD
  /**
   * @brief SIMD-accelerated vector addition.
   * @param rhs Vector to add.
   * @return Sum vector.
   */
  [[nodiscard]] Vector3 operator+(const Vector3 &rhs) const noexcept {
    Vector3 res;
    _mm256_store_pd(res.data_.data(), _mm256_add_pd(_mm256_load_pd(data_.data()), _mm256_load_pd(rhs.data_.data())));
    return res;
  }

  /**
   * @brief SIMD-accelerated vector subtraction.
   * @param rhs Vector to subtract.
   * @return Difference vector.
   */
  [[nodiscard]] Vector3 operator-(const Vector3 &rhs) const noexcept {
    Vector3 res;
    _mm256_store_pd(res.data_.data(), _mm256_sub_pd(_mm256_load_pd(data_.data()), _mm256_load_pd(rhs.data_.data())));
    return res;
  }

  /**
   * @brief SIMD-accelerated unary negation.
   * @return Negated vector.
   */
  [[nodiscard]] Vector3 operator-() const noexcept {
    Vector3 res;
    _mm256_store_pd(res.data_.data(), _mm256_sub_pd(_mm256_setzero_pd(), _mm256_load_pd(data_.data())));
    return res;
  }

  /**
   * @brief SIMD-accelerated scalar multiplication.
   * @param scalar Scalar factor.
   * @return Scaled vector.
   */
  [[nodiscard]] Vector3 operator*(double scalar) const noexcept {
    Vector3 res;
    _mm256_store_pd(res.data_.data(), _mm256_mul_pd(_mm256_load_pd(data_.data()), _mm256_set1_pd(scalar)));
    return res;
  }

  /**
   * @brief SIMD-accelerated scalar division.
   * @param scalar Scalar divisor.
   * @return Vector scaled by 1/scalar.
   */
  [[nodiscard]] Vector3 operator/(double scalar) const noexcept {
    Vector3 res;
    _mm256_store_pd(res.data_.data(), _mm256_div_pd(_mm256_load_pd(data_.data()), _mm256_set1_pd(scalar)));
    return res;
  }

  /**
   * @brief SIMD-accelerated in-place addition.
   * @param rhs Vector to add.
   * @return Reference to this vector.
   */
  Vector3 &operator+=(const Vector3 &rhs) noexcept {
    _mm256_store_pd(data_.data(), _mm256_add_pd(_mm256_load_pd(data_.data()), _mm256_load_pd(rhs.data_.data())));
    return *this;
  }

  /**
   * @brief SIMD-accelerated in-place subtraction.
   * @param rhs Vector to subtract.
   * @return Reference to this vector.
   */
  Vector3 &operator-=(const Vector3 &rhs) noexcept {
    _mm256_store_pd(data_.data(), _mm256_sub_pd(_mm256_load_pd(data_.data()), _mm256_load_pd(rhs.data_.data())));
    return *this;
  }

  /**
   * @brief SIMD-accelerated in-place scalar multiplication.
   * @param scalar Scalar factor.
   * @return Reference to this vector.
   */
  Vector3 &operator*=(double scalar) noexcept {
    _mm256_store_pd(data_.data(), _mm256_mul_pd(_mm256_load_pd(data_.data()), _mm256_set1_pd(scalar)));
    return *this;
  }

  /**
   * @brief SIMD-accelerated in-place scalar division.
   * @param scalar Scalar divisor.
   * @return Reference to this vector.
   */
  Vector3 &operator/=(double scalar) noexcept {
    _mm256_store_pd(data_.data(), _mm256_div_pd(_mm256_load_pd(data_.data()), _mm256_set1_pd(scalar)));
    return *this;
  }

  /**
   * @brief SIMD-accelerated dot product.
   * @param rhs Vector to multiply with.
   * @return Scalar result.
   */
  [[nodiscard]] double operator*(const Vector3 &rhs) const noexcept {
    __m256d mult = _mm256_mul_pd(_mm256_load_pd(data_.data()), _mm256_load_pd(rhs.data_.data()));
    __m128d low = _mm256_castpd256_pd128(mult);
    __m128d high = _mm256_extractf128_pd(mult, 1);
    __m128d res = _mm_add_pd(low, high);
    res = _mm_hadd_pd(res, res);
    return _mm_cvtsd_f64(res);
  }

  /**
   * @brief Converts the vector to a std::array.
   * @return Array containing {x, y, z}.
   */
  [[nodiscard]] std::array<double, 3> array() const noexcept { return {data_[0], data_[1], data_[2]}; }

  /**
   * @brief Equality comparison.
   * @param rhs The vector to compare with.
   * @return True if all components are exactly equal.
   */
  [[nodiscard]] bool operator==(const Vector3 &rhs) const noexcept {
    return data_[0] == rhs.data_[0] && data_[1] == rhs.data_[1] && data_[2] == rhs.data_[2];
  }

  /**
   * @brief Inequality comparison.
   * @param rhs The vector to compare with.
   * @return True if any component differs.
   */
  [[nodiscard]] bool operator!=(const Vector3 &rhs) const noexcept { return !(*this == rhs); }

  /** @return Constant pointer to the beginning of the SIMD-aligned data. */
  [[nodiscard]] const double *begin() const noexcept { return data_.data(); }
  /** @return Constant pointer to the end of the data. */
  [[nodiscard]] const double *end() const noexcept { return data_.data() + 3; }
  /** @return Mutable pointer to the beginning of the SIMD-aligned data. */
  [[nodiscard]] double *begin() noexcept { return data_.data(); }
  /** @return Mutable pointer to the end of the data. */
  [[nodiscard]] double *end() noexcept { return data_.data() + 3; }

private:
  CORRELATION_ALIGN(32) std::array<double, 4> data_ {}; ///< Padded to 4 doubles for 256-bit SIMD alignment.
};

/**
 * @brief SIMD-optimized specialization of Vector3 for float.
 */
template <> class CORRELATION_ALIGN(16) Vector3<float> {
public:
  using value_type = float; ///< Float scalar type.

  /** @brief Default constructor. Initializes to zero and ensures 16-byte alignment. */
  Vector3() noexcept { _mm_store_ps(data_.data(), _mm_setzero_ps()); }

  /**
   * @brief Parameterized constructor.
   * @param x_val X-component.
   * @param y_val Y-component.
   * @param z_val Z-component.
   */
  Vector3(float x_val, float y_val, float z_val) noexcept {
    CORRELATION_ALIGN(16) std::array<float, 4> temp = {x_val, y_val, z_val, 0.0F};
    _mm_store_ps(data_.data(), _mm_load_ps(temp.data()));
  }

  /**
   * @brief Array-based constructor.
   * @param arr Input array {x, y, z}.
   */
  explicit Vector3(const std::array<float, 3> &arr) noexcept {
    CORRELATION_ALIGN(16) std::array<float, 4> temp = {arr[0], arr[1], arr[2], 0.0F};
    _mm_store_ps(data_.data(), _mm_load_ps(temp.data()));
  }

  /**
   * @brief Converting constructor from Vector3<double>.
   */
  explicit Vector3(const Vector3<double> &other) noexcept {
    data_[0] = static_cast<float>(other.x());
    data_[1] = static_cast<float>(other.y());
    data_[2] = static_cast<float>(other.z());
    data_[3] = 0.0F;
  }

  /**
   * @brief Access component by index.
   * @param idx Index (0, 1, or 2).
   * @return Component value.
   */
  [[nodiscard]] float operator[](std::size_t idx) const noexcept {
    return data_[idx]; // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
  }

  /**
   * @brief Mutable access component by index.
   * @param idx Index (0, 1, or 2).
   * @return Reference to component value.
   */
  float &operator[](std::size_t idx) noexcept {
    return data_[idx]; // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
  }

  /** @return The X-component. */
  [[nodiscard]] float x() const noexcept { return data_[0]; }
  /** @return Reference to the X-component. */
  float &x() noexcept { return data_[0]; }
  /** @return The Y-component. */
  [[nodiscard]] float y() const noexcept { return data_[1]; }
  /** @return Reference to the Y-component. */
  float &y() noexcept { return data_[1]; }
  /** @return The Z-component. */
  [[nodiscard]] float z() const noexcept { return data_[2]; }
  /** @return Reference to the Z-component. */
  float &z() noexcept { return data_[2]; }

  /**
   * @brief Access component by index (alternative syntax).
   * @param idx Index.
   * @return Reference to component.
   */
  float &operator()(std::size_t idx) noexcept {
    return data_[idx]; // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
  }

  /**
   * @brief Vector addition.
   * @param rhs The vector to add.
   * @return A new vector representing the sum.
   */
  [[nodiscard]] Vector3 operator+(const Vector3 &rhs) const noexcept {
    Vector3 result;
    __m128 vec_a = _mm_load_ps(data_.data());
    __m128 vec_b = _mm_load_ps(rhs.data_.data());
    _mm_store_ps(result.data_.data(), _mm_add_ps(vec_a, vec_b));
    return result;
  }

  /**
   * @brief Vector subtraction.
   * @param rhs The vector to subtract.
   * @return A new vector representing the difference.
   */
  [[nodiscard]] Vector3 operator-(const Vector3 &rhs) const noexcept {
    Vector3 result;
    __m128 vec_a = _mm_load_ps(data_.data());
    __m128 vec_b = _mm_load_ps(rhs.data_.data());
    _mm_store_ps(result.data_.data(), _mm_sub_ps(vec_a, vec_b));
    return result;
  }

  /**
   * @brief Unary negation.
   * @return A new vector with all components negated.
   */
  [[nodiscard]] Vector3 operator-() const noexcept {
    Vector3 result;
    __m128 vec_a = _mm_load_ps(data_.data());
    _mm_store_ps(result.data_.data(), _mm_sub_ps(_mm_setzero_ps(), vec_a));
    return result;
  }

  /**
   * @brief Scalar multiplication.
   * @param scalar The scalar factor.
   * @return A new vector scaled by scalar.
   */
  [[nodiscard]] Vector3 operator*(float scalar) const noexcept {
    Vector3 result;
    __m128 vec_a = _mm_load_ps(data_.data());
    __m128 sclr = _mm_set1_ps(scalar);
    _mm_store_ps(result.data_.data(), _mm_mul_ps(vec_a, sclr));
    return result;
  }

  /**
   * @brief Scalar division.
   * @param scalar The scalar divisor.
   * @return A new vector scaled by 1/scalar.
   */
  [[nodiscard]] Vector3 operator/(float scalar) const noexcept {
    Vector3 result;
    __m128 vec_a = _mm_load_ps(data_.data());
    __m128 sclr = _mm_set1_ps(scalar);
    _mm_store_ps(result.data_.data(), _mm_div_ps(vec_a, sclr));
    return result;
  }

  /** @brief In-place addition. */
  Vector3 &operator+=(const Vector3 &rhs) noexcept {
    _mm_store_ps(data_.data(), _mm_add_ps(_mm_load_ps(data_.data()), _mm_load_ps(rhs.data_.data())));
    return *this;
  }

  /** @brief In-place subtraction. */
  Vector3 &operator-=(const Vector3 &rhs) noexcept {
    _mm_store_ps(data_.data(), _mm_sub_ps(_mm_load_ps(data_.data()), _mm_load_ps(rhs.data_.data())));
    return *this;
  }

  /** @brief In-place scalar multiplication. */
  Vector3 &operator*=(float scalar) noexcept {
    _mm_store_ps(data_.data(), _mm_mul_ps(_mm_load_ps(data_.data()), _mm_set1_ps(scalar)));
    return *this;
  }

  /** @brief In-place scalar division. */
  Vector3 &operator/=(float scalar) noexcept {
    _mm_store_ps(data_.data(), _mm_div_ps(_mm_load_ps(data_.data()), _mm_set1_ps(scalar)));
    return *this;
  }

  /**
   * @brief Dot product.
   * @param rhs The vector to dot with.
   * @return The scalar dot product.
   */
  [[nodiscard]] float operator*(const Vector3 &rhs) const noexcept {
    return data_[0] * rhs.data_[0] + data_[1] * rhs.data_[1] + data_[2] * rhs.data_[2];
  }

  /** @return Array containing {x, y, z}. */
  [[nodiscard]] std::array<float, 3> array() const noexcept { return {data_[0], data_[1], data_[2]}; }

  /** @brief Checks if the vector is a zero vector. */
  [[nodiscard]] bool empty() const noexcept { return data_[0] == 0.0F && data_[1] == 0.0F && data_[2] == 0.0F; }

  /** @brief Equality comparison. */
  [[nodiscard]] bool operator==(const Vector3 &rhs) const noexcept {
    return data_[0] == rhs[0] && data_[1] == rhs[1] && data_[2] == rhs[2];
  }

  /** @brief Inequality comparison. */
  [[nodiscard]] bool operator!=(const Vector3 &rhs) const noexcept { return !(*this == rhs); }

  /** @return Pointer to the beginning of the data. */
  [[nodiscard]] const float *begin() const noexcept { return data_.data(); }
  /** @return Pointer past the end of the data. */
  [[nodiscard]] const float *end() const noexcept { return data_.data() + 3; }
  /** @return Mutable pointer to the beginning of the data. */
  [[nodiscard]] float *begin() noexcept { return data_.data(); }
  /** @return Mutable pointer to the end of the data. */
  [[nodiscard]] float *end() noexcept { return data_.data() + 3; }

private:
  CORRELATION_ALIGN(16) std::array<float, 4> data_ {}; ///< Padded to 4 floats for 128-bit SSE alignment.
};

#endif // SIMD Specialized Vector3<float>

} // namespace correlation::math
