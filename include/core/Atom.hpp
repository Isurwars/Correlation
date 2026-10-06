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
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <unordered_map>

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
  constexpr bool operator==(const ElementID &other) const = default;
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
  constexpr bool operator==(const Element &other) const { return symbol == other.symbol; }
};

/**
 * @brief Thread-safe flyweight pool for chemical Element instances.
 */
class ElementPool {
public:
  /**
   * @brief Returns a shared sentinel empty element instance.
   * @return Pointer to default element.
   */
  [[nodiscard]] static const Element *defaultElement() noexcept {
    static const Element DEFAULT_ELEM{.symbol = "", .id = ElementID{-1}};
    return &DEFAULT_ELEM;
  }

  /**
   * @brief Interns an Element into the pool or returns an existing instance.
   * @param elem Element struct to intern.
   * @return Immutable pointer to shared Element.
   */
  [[nodiscard]] static const Element *intern(const Element &elem) {
    if (elem.symbol.empty() && elem.id.value == -1) {
      return defaultElement();
    }
    return intern(elem.symbol, elem.id);
  }

  /**
   * @brief Interns an element by symbol and ID into the pool.
   * @param symbol Chemical symbol.
   * @param elem_id Element identifier.
   * @return Immutable pointer to shared Element.
   */
  [[nodiscard]] static const Element *intern(std::string_view symbol, ElementID elem_id) {
    struct Key {
      std::string symbol;
      int id_val;
      bool operator==(const Key &other) const noexcept {
        return id_val == other.id_val && symbol == other.symbol;
      }
    };
    struct KeyHash {
      size_t operator()(const Key &key_obj) const noexcept {
        return std::hash<std::string_view>{}(key_obj.symbol) ^
               (std::hash<int>{}(key_obj.id_val) << 1);
      }
    };

    static std::mutex pool_mutex;
    static std::unordered_map<Key, std::unique_ptr<Element>, KeyHash> pool;

    std::scoped_lock lock(pool_mutex);
    Key const key{std::string(symbol), elem_id.value};
    auto iter = pool.find(key);
    if (iter != pool.end()) {
      return iter->second.get();
    }

    auto new_element = std::make_unique<Element>(std::string(symbol), elem_id);
    const Element *ptr = new_element.get();
    pool.emplace(key, std::move(new_element));
    return ptr;
  }
};

/**
 * @brief Represents an atom in the simulation cell.
 *
 * Stores the element type as a flyweight pointer, position, and unique ID of the atom.
 */
class Atom {
public:
  /** @name Constructors */
  ///@{
  explicit Atom() noexcept : element_(ElementPool::defaultElement()) {}

  /**
   * @brief Parameterized constructor.
   * @param element The element type of the atom.
   * @param pos The position vector of the atom.
   * @param atom_id The unique ID of the atom.
   */
  explicit Atom(const Element &element, const math::Vector3<real_t> &pos, AtomID atom_id) noexcept
      : element_(ElementPool::intern(element)), position_(pos), id_(atom_id) {}

  /**
   * @brief Pointer-based flyweight constructor.
   * @param element_ptr Pointer to the interned element.
   * @param pos The position vector of the atom.
   * @param atom_id The unique ID of the atom.
   */
  explicit Atom(const Element *element_ptr, const math::Vector3<real_t> &pos,
                AtomID atom_id) noexcept
      : element_(element_ptr != nullptr ? element_ptr : ElementPool::defaultElement()),
        position_(pos), id_(atom_id) {}

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
  void setID(std::uint32_t num) noexcept { id_ = num; }

  /**
   * @brief Gets the position of the atom.
   * @return A const reference to the position vector.
   */
  [[nodiscard]] const math::Vector3<real_t> &position() const noexcept { return position_; }

  /**
   * @brief Sets the position of the atom.
   * @param pos The new position vector.
   */
  void setPosition(const math::Vector3<real_t> &pos) noexcept { position_ = pos; }

  /**
   * @brief Gets the velocity of the atom.
   * @return A const reference to the velocity vector.
   */
  [[nodiscard]] const math::Vector3<real_t> &velocity() const noexcept { return velocity_; }

  /**
   * @brief Sets the velocity of the atom.
   * @param vel The new velocity vector.
   */
  void setVelocity(const math::Vector3<real_t> &vel) noexcept { velocity_ = vel; }

  /**
   * @brief Gets the element type of the atom.
   * @return A const reference to the Element struct.
   */
  [[nodiscard]] const Element &element() const noexcept {
    return element_ != nullptr ? *element_ : *ElementPool::defaultElement();
  }

  /**
   * @brief Sets the element type of the atom.
   * @param ele The new Element struct.
   */
  void setElement(const Element &ele) { element_ = ElementPool::intern(ele); }

  /**
   * @brief Gets the integer ID of the element type.
   * @return The element ID value.
   */
  [[nodiscard]] int elementId() const noexcept { return element().id.value; }

  ///@}

private:
  AtomID id_{0};                   ///< Unique identification number.
  math::Vector3<real_t> position_; ///< Cartesian coordinates in Angstroms.
  math::Vector3<real_t> velocity_; ///< Velocity in Angstroms/fs.
  const Element *element_{
      ElementPool::defaultElement()}; ///< Flyweight pointer to element metadata.
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
[[nodiscard]] inline real_t angle(const Atom &center, const Atom &atom_a,
                                  const Atom &atom_b) noexcept {
  const math::Vector3<real_t> vec_a = atom_a.position() - center.position();
  const math::Vector3<real_t> vec_b = atom_b.position() - center.position();

  const real_t norm_sq_a = math::dot(vec_a, vec_a);
  const real_t norm_sq_b = math::dot(vec_b, vec_b);

  // Guard against near-zero norms to prevent division by zero or NaN underflows
  if (norm_sq_a < 1e-12 || norm_sq_b < 1e-12) {
    return 0.0;
  }

  using std::acos;
  using std::sqrt;
  real_t cos_theta = math::dot(vec_a, vec_b) / sqrt(norm_sq_a * norm_sq_b);

  // Clamp for numerical stability
  cos_theta = std::clamp(cos_theta, static_cast<real_t>(-1.0), static_cast<real_t>(1.0));

  return acos(cos_theta);
}

} // namespace correlation::core
