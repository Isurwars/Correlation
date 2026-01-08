// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <algorithm>
#include <cstdint>
#include <string>

#include "LinearAlgebra.hpp"

using AtomID = std::uint32_t;

struct ElementID {
  int value;
  constexpr bool operator==(const ElementID other) const {
    return value == other.value;
  };
};

struct Element {
  std::string symbol;
  ElementID id{-1};
};

struct Neighbor {
  AtomID index;
  double distance;
  linalg::Vector3<double> r_ij;
};

class Atom {
public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//
  explicit Atom() = default;

  explicit Atom(Element element, linalg::Vector3<double> pos,
                AtomID id) noexcept
      : element_(std::move(element)), position_(pos), id_(id) {}

  //-------------------------------------------------------------------------//
  //------------------------------- Accessors -------------------------------//
  //-------------------------------------------------------------------------//

  AtomID id() const noexcept { return id_; }
  void setID(std::uint32_t num) { id_ = num; }

  const linalg::Vector3<double> &position() const noexcept { return position_; }
  void setPosition(linalg::Vector3<double> pos) { position_ = pos; }

  const Element &element() const { return element_; }
  void setElement(const Element &ele) { element_ = ele; }

  int element_id() const { return element_.id.value; }

private:
  AtomID id_;
  linalg::Vector3<double> position_;
  Element element_;
};

/**
 * @brief Calculates the Euclidean distance between two atoms.
 */
inline double distance(const Atom &a, const Atom &b) noexcept {
  return linalg::norm(a.position() - b.position());
}

/**
 * @brief Calculates the angle (in radians) formed by three atoms.
 * @param center The atom at the vertex of the angle.
 * @param a One of the outer atoms.
 * @param b The other outer atom.
 * @return The angle in radians, or 0.0 if vectors are collinear or zero.
 */
inline double angle(const Atom &center, const Atom &a, const Atom &b) noexcept {
  const linalg::Vector3<double> vA = a.position() - center.position();
  const linalg::Vector3<double> vB = b.position() - center.position();

  const double norm_sq_A = linalg::dot(vA, vA);
  const double norm_sq_B = linalg::dot(vB, vB);

  if (norm_sq_A == 0.0 || norm_sq_B == 0.0) {
    return 0.0;
  }

  double cos_theta = linalg::dot(vA, vB) / std::sqrt(norm_sq_A * norm_sq_B);

  // Clamp for numerical stability
  cos_theta = std::clamp(cos_theta, -1.0, 1.0);

  return std::acos(cos_theta);
}
