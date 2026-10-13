/**
 * @file Math3D.hpp
 * @brief Self-contained 3D transformation mathematics (Vector4, Matrix4, Quaternion).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/LinearAlgebra.hpp"
#include <array>
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace correlation::render3d {

// -----------------------------------------------------------------------------
// Vector4
// -----------------------------------------------------------------------------

/**
 * @class Vector4
 * @brief Stack-allocated 4D homogeneous vector.
 * @tparam T Floating-point scalar type.
 */
template <typename T = float>
  requires std::is_floating_point_v<T>
class Vector4 {
public:
  using value_type = T;

  /** @brief Default constructor. Initializes all components to zero. */
  constexpr Vector4() noexcept : data_{0, 0, 0, 0} {}

  /**
   * @brief Component-wise constructor.
   * @param coord_x X coordinate.
   * @param coord_y Y coordinate.
   * @param coord_z Z coordinate.
   * @param coord_w W coordinate (homogeneous scale).
   */
  constexpr Vector4(T coord_x, T coord_y, T coord_z, T coord_w) noexcept
      : data_{coord_x, coord_y, coord_z, coord_w} {}

  /**
   * @brief Constructs Vector4 from 3D vector and w component.
   * @param vec3 3D vector.
   * @param coord_w Homogeneous coordinate (default 1.0).
   */
  constexpr explicit Vector4(const math::Vector3<T> &vec3, T coord_w = static_cast<T>(1)) noexcept
      : data_{vec3.x(), vec3.y(), vec3.z(), coord_w} {}

  /**
   * @brief Constructs Vector4 from std::array.
   * @param arr Array containing {x, y, z, w}.
   */
  constexpr explicit Vector4(const std::array<T, 4> &arr) noexcept : data_{arr} {}

  [[nodiscard]] constexpr T x() const noexcept { return data_.at(0); }
  [[nodiscard]] constexpr T y() const noexcept { return data_.at(1); }
  [[nodiscard]] constexpr T z() const noexcept { return data_.at(2); }
  [[nodiscard]] constexpr T w() const noexcept { return data_.at(3); }

  [[nodiscard]] constexpr T &operator[](std::size_t idx) noexcept { return data_.at(idx); }
  [[nodiscard]] constexpr const T &operator[](std::size_t idx) const noexcept {
    return data_.at(idx);
  }

  [[nodiscard]] constexpr const std::array<T, 4> &data() const noexcept { return data_; }

  [[nodiscard]] constexpr math::Vector3<T> toVector3() const noexcept {
    if (std::abs(data_.at(3)) > static_cast<T>(1e-7)) {
      return math::Vector3<T>(data_.at(0) / data_.at(3), data_.at(1) / data_.at(3),
                              data_.at(2) / data_.at(3));
    }
    return math::Vector3<T>(data_.at(0), data_.at(1), data_.at(2));
  }

  [[nodiscard]] constexpr Vector4 operator+(const Vector4 &rhs) const noexcept {
    return Vector4(data_.at(0) + rhs.data_.at(0), data_.at(1) + rhs.data_.at(1),
                   data_.at(2) + rhs.data_.at(2), data_.at(3) + rhs.data_.at(3));
  }

  [[nodiscard]] constexpr Vector4 operator-(const Vector4 &rhs) const noexcept {
    return Vector4(data_.at(0) - rhs.data_.at(0), data_.at(1) - rhs.data_.at(1),
                   data_.at(2) - rhs.data_.at(2), data_.at(3) - rhs.data_.at(3));
  }

  [[nodiscard]] constexpr Vector4 operator*(T scalar) const noexcept {
    return Vector4(data_.at(0) * scalar, data_.at(1) * scalar, data_.at(2) * scalar,
                   data_.at(3) * scalar);
  }

  [[nodiscard]] constexpr Vector4 operator/(T scalar) const {
    if (std::abs(scalar) < static_cast<T>(1e-9)) {
      throw std::runtime_error("Division by zero in Vector4.");
    }
    const T inv = static_cast<T>(1) / scalar;
    return *this * inv;
  }

  [[nodiscard]] constexpr T dot(const Vector4 &rhs) const noexcept {
    return data_.at(0) * rhs.data_.at(0) + data_.at(1) * rhs.data_.at(1) +
           data_.at(2) * rhs.data_.at(2) + data_.at(3) * rhs.data_.at(3);
  }

  [[nodiscard]] T length() const noexcept { return std::sqrt(dot(*this)); }

  [[nodiscard]] Vector4 normalized() const {
    const T len = length();
    if (len < static_cast<T>(1e-9)) {
      return Vector4();
    }
    return *this * (static_cast<T>(1) / len);
  }

  [[nodiscard]] constexpr bool operator==(const Vector4 &rhs) const noexcept {
    return data_ == rhs.data_;
  }

private:
  std::array<T, 4> data_{0, 0, 0, 0};
};

// -----------------------------------------------------------------------------
// Matrix4 (Column-Major)
// -----------------------------------------------------------------------------

/**
 * @class Matrix4
 * @brief 4x4 column-major transformation matrix.
 * @tparam T Floating-point scalar type.
 */
template <typename T = float>
  requires std::is_floating_point_v<T>
class Matrix4 {
public:
  using value_type = T;

  /** @brief Default constructor. Initializes to zero matrix. */
  constexpr Matrix4() noexcept = default;

  /**
   * @brief Constructs matrix from column-major elements array.
   * @param raw 16-element array in column-major order.
   */
  constexpr explicit Matrix4(const std::array<T, 16> &raw) noexcept : data_{raw} {}

  /**
   * @brief Creates a 4x4 identity matrix.
   * @return Identity matrix.
   */
  [[nodiscard]] static constexpr Matrix4 identity() noexcept {
    Matrix4 result_mat;
    result_mat.data_.at(0) = static_cast<T>(1);
    result_mat.data_.at(5) = static_cast<T>(1);
    result_mat.data_.at(10) = static_cast<T>(1);
    result_mat.data_.at(15) = static_cast<T>(1);
    return result_mat;
  }

  /**
   * @brief Element accessor by (row, column).
   * @param row_idx Row index in [0, 3].
   * @param col_idx Column index in [0, 3].
   * @return Element reference.
   */
  [[nodiscard]] constexpr T &operator()(std::size_t row_idx, std::size_t col_idx) noexcept {
    return data_.at(col_idx * 4 + row_idx);
  }

  [[nodiscard]] constexpr const T &operator()(std::size_t row_idx,
                                              std::size_t col_idx) const noexcept {
    return data_.at(col_idx * 4 + row_idx);
  }

  [[nodiscard]] constexpr const std::array<T, 16> &data() const noexcept { return data_; }
  [[nodiscard]] constexpr std::array<T, 16> &data() noexcept { return data_; }

  /**
   * @brief Matrix-vector multiplication.
   * @param vec 4D vector.
   * @return Transformed vector.
   */
  [[nodiscard]] constexpr Vector4<T> operator*(const Vector4<T> &vec) const noexcept {
    return Vector4<T>((*this)(0, 0) * vec.x() + (*this)(0, 1) * vec.y() + (*this)(0, 2) * vec.z() +
                          (*this)(0, 3) * vec.w(),
                      (*this)(1, 0) * vec.x() + (*this)(1, 1) * vec.y() + (*this)(1, 2) * vec.z() +
                          (*this)(1, 3) * vec.w(),
                      (*this)(2, 0) * vec.x() + (*this)(2, 1) * vec.y() + (*this)(2, 2) * vec.z() +
                          (*this)(2, 3) * vec.w(),
                      (*this)(3, 0) * vec.x() + (*this)(3, 1) * vec.y() + (*this)(3, 2) * vec.z() +
                          (*this)(3, 3) * vec.w());
  }

  /**
   * @brief Matrix-matrix multiplication.
   * @param rhs Right-hand side matrix.
   * @return Product matrix.
   */
  [[nodiscard]] constexpr Matrix4 operator*(const Matrix4 &rhs) const noexcept {
    Matrix4 result_mat;
    for (std::size_t col_idx = 0; col_idx < 4; ++col_idx) {
      for (std::size_t row_idx = 0; row_idx < 4; ++row_idx) {
        T sum = 0;
        for (std::size_t elem_idx = 0; elem_idx < 4; ++elem_idx) {
          sum += (*this)(row_idx, elem_idx) * rhs(elem_idx, col_idx);
        }
        result_mat(row_idx, col_idx) = sum;
      }
    }
    return result_mat;
  }

  /**
   * @brief Transforms a 3D point (w = 1) with perspective division.
   * @param point 3D point.
   * @return Projected 3D point.
   */
  [[nodiscard]] constexpr math::Vector3<T>
  transformPoint(const math::Vector3<T> &point) const noexcept {
    const Vector4<T> vec(point, static_cast<T>(1));
    const Vector4<T> transformed = *this * vec;
    return transformed.toVector3();
  }

  /**
   * @brief Transforms a 3D direction vector (w = 0).
   * @param direction 3D direction.
   * @return Transformed direction.
   */
  [[nodiscard]] constexpr math::Vector3<T>
  transformDirection(const math::Vector3<T> &direction) const noexcept {
    const Vector4<T> vec(direction, static_cast<T>(0));
    const Vector4<T> transformed = *this * vec;
    return math::Vector3<T>(transformed.x(), transformed.y(), transformed.z());
  }

  /**
   * @brief Creates a translation matrix.
   * @param offset_x Translation in X.
   * @param offset_y Translation in Y.
   * @param offset_z Translation in Z.
   * @return Translation matrix.
   */
  [[nodiscard]] static constexpr Matrix4 translation(T offset_x, T offset_y, T offset_z) noexcept {
    Matrix4 result_mat = Matrix4::identity();
    result_mat(0, 3) = offset_x;
    result_mat(1, 3) = offset_y;
    result_mat(2, 3) = offset_z;
    return result_mat;
  }

  [[nodiscard]] static constexpr Matrix4 translation(const math::Vector3<T> &offset) noexcept {
    return translation(offset.x(), offset.y(), offset.z());
  }

  /**
   * @brief Creates a non-uniform scaling matrix.
   * @param scale_x Scale in X.
   * @param scale_y Scale in Y.
   * @param scale_z Scale in Z.
   * @return Scaling matrix.
   */
  [[nodiscard]] static constexpr Matrix4 scaling(T scale_x, T scale_y, T scale_z) noexcept {
    Matrix4 result_mat;
    result_mat(0, 0) = scale_x;
    result_mat(1, 1) = scale_y;
    result_mat(2, 2) = scale_z;
    result_mat(3, 3) = static_cast<T>(1);
    return result_mat;
  }

  [[nodiscard]] static constexpr Matrix4 scaling(T uniform_scale) noexcept {
    return scaling(uniform_scale, uniform_scale, uniform_scale);
  }

  /**
   * @brief Transposes the matrix.
   * @return Transposed matrix.
   */
  [[nodiscard]] constexpr Matrix4 transpose() const noexcept {
    Matrix4 result_mat;
    for (std::size_t row_idx = 0; row_idx < 4; ++row_idx) {
      for (std::size_t col_idx = 0; col_idx < 4; ++col_idx) {
        result_mat(col_idx, row_idx) = (*this)(row_idx, col_idx);
      }
    }
    return result_mat;
  }

  /**
   * @brief Generates a Look-At view matrix (Right-Handed).
   * @param eye Camera eye position.
   * @param target Target look-at point.
   * @param up_direction World up vector.
   * @return View matrix.
   */
  [[nodiscard]] static Matrix4 lookAt(const math::Vector3<T> &eye, const math::Vector3<T> &target,
                                      const math::Vector3<T> &up_direction) {
    const math::Vector3<T> forward = math::normalize(target - eye);
    const math::Vector3<T> right = math::normalize(math::cross(forward, up_direction));
    const math::Vector3<T> actual_up = math::cross(right, forward);

    Matrix4 result_mat = Matrix4::identity();
    result_mat(0, 0) = right.x();
    result_mat(0, 1) = right.y();
    result_mat(0, 2) = right.z();
    result_mat(0, 3) = -math::dot(right, eye);

    result_mat(1, 0) = actual_up.x();
    result_mat(1, 1) = actual_up.y();
    result_mat(1, 2) = actual_up.z();
    result_mat(1, 3) = -math::dot(actual_up, eye);

    result_mat(2, 0) = -forward.x();
    result_mat(2, 1) = -forward.y();
    result_mat(2, 2) = -forward.z();
    result_mat(2, 3) = math::dot(forward, eye);

    return result_mat;
  }

  /**
   * @brief Generates perspective projection matrix (RH, depth in [0, 1] suitable for
   * Vulkan/Software).
   * @param fov_y_rad Vertical field of view in radians.
   * @param aspect Viewport aspect ratio (width / height).
   * @param z_near Near clipping plane distance (> 0).
   * @param z_far Far clipping plane distance (> z_near).
   * @return Perspective projection matrix.
   */
  [[nodiscard]] static Matrix4 perspective(T fov_y_rad, T aspect, T z_near, T z_far) {
    if (std::abs(aspect) < static_cast<T>(1e-6) || z_near <= 0 || z_far <= z_near) {
      throw std::invalid_argument("Invalid projection parameters in perspective().");
    }
    const T tan_half_fov = std::tan(fov_y_rad / static_cast<T>(2));
    Matrix4 result_mat;
    result_mat(0, 0) = static_cast<T>(1) / (aspect * tan_half_fov);
    result_mat(1, 1) = static_cast<T>(1) / tan_half_fov;
    result_mat(2, 2) = z_far / (z_near - z_far);
    result_mat(2, 3) = (z_near * z_far) / (z_near - z_far);
    result_mat(3, 2) = static_cast<T>(-1);
    return result_mat;
  }

  /**
   * @brief Generates orthographic projection matrix (RH, depth in [0, 1]).
   * @param left Left boundary.
   * @param right Right boundary.
   * @param bottom Bottom boundary.
   * @param top Top boundary.
   * @param z_near Near plane.
   * @param z_far Far plane.
   * @return Orthographic matrix.
   */
  [[nodiscard]] static Matrix4 orthographic(T left, T right, T bottom, T top, T z_near, T z_far) {
    Matrix4 result_mat = Matrix4::identity();
    result_mat(0, 0) = static_cast<T>(2) / (right - left);
    result_mat(1, 1) = static_cast<T>(2) / (top - bottom);
    result_mat(2, 2) = static_cast<T>(1) / (z_near - z_far);
    result_mat(0, 3) = -(right + left) / (right - left);
    result_mat(1, 3) = -(top + bottom) / (top - bottom);
    result_mat(2, 3) = z_near / (z_near - z_far);
    return result_mat;
  }

private:
  std::array<T, 16> data_{};
};

// -----------------------------------------------------------------------------
// Quaternion
// -----------------------------------------------------------------------------

/**
 * @class Quaternion
 * @brief Unit quaternion for robust 3D rotation representation without gimbal lock.
 * @tparam T Floating-point scalar type.
 */
template <typename T = float>
  requires std::is_floating_point_v<T>
class Quaternion {
public:
  using value_type = T;

  /**
   * @struct Components
   * @brief Explicit bundle of quaternion scalar and imaginary components.
   */
  struct Components {
    T real_w{static_cast<T>(1)};
    T imag_x{static_cast<T>(0)};
    T imag_y{static_cast<T>(0)};
    T imag_z{static_cast<T>(0)};
  };

  /** @brief Default constructor. Initializes to identity quaternion (w = 1). */
  constexpr Quaternion() noexcept = default;

  /**
   * @brief Constructs quaternion from explicit components bundle.
   * @param components Struct holding real and imaginary parts.
   */
  constexpr explicit Quaternion(const Components &components) noexcept
      : w_{components.real_w}, x_{components.imag_x}, y_{components.imag_y}, z_{components.imag_z} {
  }

  /**
   * @brief Constructs quaternion from scalar real part and 3D imaginary vector.
   * @param real_val Real scalar component.
   * @param imag_vec Imaginary 3D vector.
   */
  constexpr Quaternion(T real_val, const math::Vector3<T> &imag_vec) noexcept
      : w_{real_val}, x_{imag_vec.x()}, y_{imag_vec.y()}, z_{imag_vec.z()} {}

  [[nodiscard]] constexpr T w() const noexcept { return w_; }
  [[nodiscard]] constexpr T x() const noexcept { return x_; }
  [[nodiscard]] constexpr T y() const noexcept { return y_; }
  [[nodiscard]] constexpr T z() const noexcept { return z_; }

  /**
   * @brief Creates quaternion from rotation axis and angle.
   * @param axis Normalized rotation axis.
   * @param angle_rad Rotation angle in radians.
   * @return Rotation quaternion.
   */
  [[nodiscard]] static Quaternion fromAxisAngle(const math::Vector3<T> &axis, T angle_rad) {
    const T half_angle = angle_rad / static_cast<T>(2);
    const T sin_half = std::sin(half_angle);
    const math::Vector3<T> norm_axis = math::normalize(axis);
    return Quaternion(std::cos(half_angle), norm_axis * sin_half);
  }

  [[nodiscard]] T length() const noexcept {
    return std::sqrt(w_ * w_ + x_ * x_ + y_ * y_ + z_ * z_);
  }

  [[nodiscard]] Quaternion normalized() const {
    const T len = length();
    if (len < static_cast<T>(1e-9)) {
      return Quaternion();
    }
    const T inv_len = static_cast<T>(1) / len;
    return Quaternion(w_ * inv_len, math::Vector3<T>(x_ * inv_len, y_ * inv_len, z_ * inv_len));
  }

  [[nodiscard]] constexpr Quaternion conjugate() const noexcept {
    return Quaternion(w_, math::Vector3<T>(-x_, -y_, -z_));
  }

  /**
   * @brief Hamilton quaternion multiplication (q1 * q2).
   * @param rhs Right-hand side quaternion.
   * @return Combined rotation.
   */
  [[nodiscard]] constexpr Quaternion operator*(const Quaternion &rhs) const noexcept {
    const T new_w = w_ * rhs.w_ - x_ * rhs.x_ - y_ * rhs.y_ - z_ * rhs.z_;
    const math::Vector3<T> new_imag(w_ * rhs.x_ + x_ * rhs.w_ + y_ * rhs.z_ - z_ * rhs.y_,
                                    w_ * rhs.y_ - x_ * rhs.z_ + y_ * rhs.w_ + z_ * rhs.x_,
                                    w_ * rhs.z_ + x_ * rhs.y_ - y_ * rhs.x_ + z_ * rhs.w_);
    return Quaternion(new_w, new_imag);
  }

  /**
   * @brief Rotates a 3D vector by this quaternion.
   * @param vector3 Input 3D vector.
   * @return Rotated vector.
   */
  [[nodiscard]] math::Vector3<T> rotate(const math::Vector3<T> &vector3) const {
    const Quaternion pure_quat(static_cast<T>(0), vector3);
    const Quaternion quat_norm = normalized();
    const Quaternion quat_res = quat_norm * pure_quat * quat_norm.conjugate();
    return math::Vector3<T>(quat_res.x(), quat_res.y(), quat_res.z());
  }

  /**
   * @brief Converts unit quaternion to a 4x4 rotation matrix.
   * @return Rotation matrix.
   */
  [[nodiscard]] Matrix4<T> toMatrix4() const {
    const Quaternion quat = normalized();
    Matrix4<T> result_mat = Matrix4<T>::identity();

    const T x_sq = quat.x_ * quat.x_;
    const T y_sq = quat.y_ * quat.y_;
    const T z_sq = quat.z_ * quat.z_;
    const T xy_prod = quat.x_ * quat.y_;
    const T xz_prod = quat.x_ * quat.z_;
    const T yz_prod = quat.y_ * quat.z_;
    const T wx_prod = quat.w_ * quat.x_;
    const T wy_prod = quat.w_ * quat.y_;
    const T wz_prod = quat.w_ * quat.z_;

    result_mat(0, 0) = static_cast<T>(1) - static_cast<T>(2) * (y_sq + z_sq);
    result_mat(0, 1) = static_cast<T>(2) * (xy_prod - wz_prod);
    result_mat(0, 2) = static_cast<T>(2) * (xz_prod + wy_prod);

    result_mat(1, 0) = static_cast<T>(2) * (xy_prod + wz_prod);
    result_mat(1, 1) = static_cast<T>(1) - static_cast<T>(2) * (x_sq + z_sq);
    result_mat(1, 2) = static_cast<T>(2) * (yz_prod - wx_prod);

    result_mat(2, 0) = static_cast<T>(2) * (xz_prod - wy_prod);
    result_mat(2, 1) = static_cast<T>(2) * (yz_prod + wx_prod);
    result_mat(2, 2) = static_cast<T>(1) - static_cast<T>(2) * (x_sq + y_sq);

    return result_mat;
  }

private:
  T w_{static_cast<T>(1)};
  T x_{static_cast<T>(0)};
  T y_{static_cast<T>(0)};
  T z_{static_cast<T>(0)};
};

} // namespace correlation::render3d
