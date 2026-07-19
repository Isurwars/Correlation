/**
 * @file Atom.hpp
 * @brief Atom data structure and AtomID type definitions.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"

#include <algorithm>
#include <cstdint>
#include <string>

namespace correlation::core {

/**
 * @brief Unsigned integer type used for unique atom identification.
 */
using AtomID = std::uint32_t;

/**
 * @brief Represents a unique integer ID for an element type.
 */
struct ElementID {
  int value; ///< Unique integer value representing the element type.

  /**
   * @brief Equality operator for ElementID.
   * @param other The other ID to compare against.
   * @return True if values are equal.
   */
  constexpr bool operator==(const ElementID&) const = default;
};

/**
 * @brief Represents a chemical element with its symbol and unique ID.
 */
struct Element {
  std::string symbol; ///< Chemical symbol (e.g., "Si").
  ElementID id{-1};   ///< Assigned unique integer ID.

  /**
   * @brief Equality operator for Element.
   * @param other The other element to compare against.
   * @return True if symbols are identical.
   */
  constexpr bool operator==(const Element &other) const { return symbol == other.symbol; };
};

/**
 * @brief Represents an atom in the simulation cell.
 *
 * Stores the element type, position, and unique ID of the atom.
 */
class Atom {
public:
  /** @name Constructors */
  ///@{
  explicit Atom() = default;

  /**
   * @brief Parameterized constructor.
   * @param element The element type of the atom.
   * @param pos The position vector of the atom.
   * @param atom_id The unique ID of the atom.
   */
  explicit Atom(Element element, const math::Vector3<real_t> &pos, AtomID atom_id) noexcept
      : element_(std::move(element)), position_(pos), id_(atom_id) {}

  ///@}

  /** @name Accessors */
  ///@{

  /**
   * @brief Gets the unique ID of the atom.
   * @return The atom ID.
   */
  [[nodiscard]] AtomID id() const noexcept { return id_; }

  /**
   * @brief Sets the unique ID of the atom.
   * @param num The new atom ID.
   */
  void setID(std::uint32_t num) { id_ = num; }

  /**
   * @brief Gets the position of the atom.
   * @return A const reference to the position vector.
   */
  [[nodiscard]] const math::Vector3<real_t> &position() const noexcept { return position_; }

  /**
   * @brief Sets the position of the atom.
   * @param pos The new position vector.
   */
  void setPosition(const math::Vector3<real_t> &pos) { position_ = pos; }

  /**
   * @brief Gets the velocity of the atom.
   * @return A const reference to the velocity vector.
   */
  [[nodiscard]] const math::Vector3<real_t> &velocity() const noexcept { return velocity_; }

  /**
   * @brief Sets the velocity of the atom.
   * @param vel The new velocity vector.
   */
  void setVelocity(const math::Vector3<real_t> &vel) { velocity_ = vel; }

  /**
   * @brief Gets the element type of the atom.
   * @return A const reference to the Element struct.
   */
  [[nodiscard]] const Element &element() const { return element_; }

  /**
   * @brief Sets the element type of the atom.
   * @param ele The new Element struct.
   */
  void setElement(const Element &ele) { element_ = ele; }

  /**
   * @brief Gets the integer ID of the element type.
   * @return The element ID value.
   */
  [[nodiscard]] int element_id() const { return element_.id.value; }

  ///@}

private:
  AtomID id_{0};                   ///< Unique identification number.
  math::Vector3<real_t> position_; ///< Cartesian coordinates in Angstroms.
  math::Vector3<real_t> velocity_; ///< Velocity in Angstroms/fs.
  Element element_;                ///< Chemical element properties.
};

/**
 * @brief Calculates the Euclidean distance between two atoms.
 * @param atom_a The first atom.
 * @param atom_b The second atom.
 * @return The straight-line distance between atoms (not accounting for PBC).
 */
[[nodiscard]] inline real_t distance(const Atom &atom_a, const Atom &atom_b) noexcept {
  return math::norm(atom_a.position() - atom_b.position());
}

/**
 * @brief Calculates the angle (in radians) formed by three atoms.
 * @param center The atom at the vertex of the angle.
 * @param atom_a One of the outer atoms.
 * @param atom_b The other outer atom.
 * @return The angle in radians, or 0.0 if vectors are collinear or zero.
 */
[[nodiscard]] inline real_t angle(const Atom &center, const Atom &atom_a, const Atom &atom_b) noexcept {
  const math::Vector3<real_t> vec_A = atom_a.position() - center.position();
  const math::Vector3<real_t> vec_B = atom_b.position() - center.position();

  const real_t norm_sq_A = math::dot(vec_A, vec_A);
  const real_t norm_sq_B = math::dot(vec_B, vec_B);

  // Guard against near-zero norms to prevent division by zero or NaN underflows
  if (norm_sq_A < 1e-12 || norm_sq_B < 1e-12) {
    return 0.0;
  }

  real_t cos_theta = math::dot(vec_A, vec_B) / std::sqrt(norm_sq_A * norm_sq_B);

  // Clamp for numerical stability
  cos_theta = std::clamp(cos_theta, static_cast<real_t>(-1.0), static_cast<real_t>(1.0));

  return std::acos(cos_theta);
}

} // namespace correlation::core
