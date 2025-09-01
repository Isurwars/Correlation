#ifndef INCLUDE_LINEAR_ALGEBRA_HPP_
#define INCLUDE_LINEAR_ALGEBRA_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include <array>
#include <cmath>
#include <stdexcept>

namespace linalg {

// -----------------------------------------------------------------------------
//  Vector3  –  lightweight, constexpr, stack-based
// -----------------------------------------------------------------------------
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

  // dot product
  constexpr T operator*(const Vector3 &rhs) const noexcept {
    return data_[0] * rhs[0] + data_[1] * rhs[1] + data_[2] * rhs[2];
  }

  constexpr std::array<T, 3> array() const noexcept { return data_; }

private:
  std::array<T, 3> data_;
};

// free scalar * vector
template <typename Scalar, typename T>
constexpr std::enable_if_t<std::is_arithmetic<Scalar>::value, Vector3<T>>
operator*(Scalar s, const Vector3<T> &v) noexcept {
  return v * static_cast<T>(s);
}

// -----------------------------------------------------------------------------
//  Matrix3  –  column-major storage
// -----------------------------------------------------------------------------
template <typename T> class Matrix3 {
public:
  using value_type = T;

  constexpr Matrix3() noexcept { data_.fill(Vector3<T>{0, 0, 0}); }
  constexpr Matrix3(const Vector3<T> &c0, const Vector3<T> &c1,
                    const Vector3<T> &c2) noexcept
      : data_{c0, c1, c2} {}

  constexpr const Vector3<T> &operator[](std::size_t c) const noexcept {
    return data_[c];
  }
  constexpr Vector3<T> &operator[](std::size_t c) noexcept { return data_[c]; }

  constexpr std::array<std::array<T, 3>, 3> array() const noexcept {
    return {{{data_[0][0], data_[1][0], data_[2][0]},
             {data_[0][1], data_[1][1], data_[2][1]},
             {data_[0][2], data_[1][2], data_[2][2]}}};
  }

private:
  std::array<Vector3<T>, 3> data_;
};

// -----------------------------------------------------------------------------
//  Free functions  –  drop-in API
// -----------------------------------------------------------------------------
template <typename T>
constexpr T dot(const Vector3<T> &a, const Vector3<T> &b) noexcept {
  return a * b;
}

template <typename T>
constexpr Vector3<T> cross(const Vector3<T> &a, const Vector3<T> &b) noexcept {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

template <typename T> constexpr T norm(const Vector3<T> &v) noexcept {
  return std::sqrt(v * v);
}

template <typename T> constexpr T determinant(const Matrix3<T> &m) noexcept {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

template <typename T> constexpr Matrix3<T> invertMatrix(const Matrix3<T> &m) {
  T det = determinant(m);
  if (det == T{0})
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

template <typename T>
constexpr Vector3<T> matrixVectorMultiply(const Matrix3<T> &m,
                                          const Vector3<T> &v) noexcept {
  return {m[0] * v, m[1] * v, m[2] * v};
}
} // namespace linalg
#endif // INCLUDE_LINEAR_ALGEBRA_HPP_
