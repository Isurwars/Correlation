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
#include <array>
#include <atomic>
#include <cstdint>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace correlation::core {

/**
 * @brief Unsigned integer type used for unique atom identification.
 */
using AtomID = std::uint32_t;

/**
 * @brief Represents a unique integer ID for an element type.
 */
struct ElementID {
  std::int32_t value{
      -1,
  }; ///< Unique integer value representing the element type (-1 for unassigned).

  /**
   * @brief Equality operator for ElementID.
   * @param[in] other The other ID to compare against.
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
   * @param[in] other The other element to compare against.
   * @return True if symbols are identical.
   */
  constexpr bool operator==(const Element &other) const { return symbol == other.symbol; }
};

/**
 * @brief Thread-safe flyweight pool for chemical Element instances.
 */
class ElementPool {
public:
  static constexpr std::size_t kMaxElements = 65536;

  /**
   * @brief Returns a shared sentinel empty element instance.
   * @return Pointer to default element.
   */
  [[nodiscard]] static const Element *defaultElement() noexcept {
    static const Element DEFAULT_ELEM{.symbol = "", .id = ElementID{-1}};
    return &DEFAULT_ELEM;
  }

  /**
   * @brief Lock-free O(1) lookup of an interned Element by its compact 16-bit pool index.
   * @param[in] index Compact index returned by internIndex().
   * @return Const reference to the interned Element.
   */
  [[nodiscard]] static const Element &get(std::uint16_t index) noexcept {
    const Element *ptr = table()[index].load(std::memory_order_acquire);
    return ptr != nullptr ? *ptr : *defaultElement();
  }

  /**
   * @brief Interns an Element and returns its compact 16-bit pool index.
   * @param[in] elem Element struct to intern.
   * @return Compact pool index (0 for default empty element).
   */
  [[nodiscard]] static std::uint16_t internIndex(const Element &elem) {
    if (elem.symbol.empty() && elem.id.value == -1) {
      return 0;
    }
    return internIndex(elem.symbol, elem.id);
  }

  /**
   * @brief Interns an element by symbol and ID and returns its compact 16-bit pool index.
   * @param[in] symbol Chemical symbol.
   * @param[in] elem_id Element identifier.
   * @return Compact pool index (0 for default empty element).
   */
  [[nodiscard]] static std::uint16_t internIndex(std::string_view symbol, ElementID elem_id) {
    if (symbol.empty() && elem_id.value == -1) {
      return 0;
    }
    struct Key {
      std::string symbol;
      std::int32_t id_val;
      bool operator==(const Key &other) const noexcept {
        return id_val == other.id_val && symbol == other.symbol;
      }
    };
    struct KeyHash {
      size_t operator()(const Key &key_obj) const noexcept {
        return std::hash<std::string_view>{}(key_obj.symbol) ^
               (std::hash<std::int32_t>{}(key_obj.id_val) << 1);
      }
    };

    static std::mutex pool_mutex;
    static std::unordered_map<Key, std::uint16_t, KeyHash> index_map;
    static std::vector<std::unique_ptr<Element>> storage;

    std::scoped_lock lock(pool_mutex);
    Key const key{std::string(symbol), elem_id.value};
    if (const auto iter = index_map.find(key); iter != index_map.end()) {
      return iter->second;
    }

    const std::size_t next_slot = storage.size() + 1;
    if (next_slot >= kMaxElements) {
      throw std::overflow_error("ElementPool capacity exceeded (max 65535 distinct elements).");
    }

    auto new_element =
        std::unique_ptr<Element>(new Element{.symbol = std::string(symbol), .id = elem_id});
    const Element *ptr = new_element.get();
    storage.push_back(std::move(new_element));
    const auto idx = static_cast<std::uint16_t>(next_slot);
    table()[idx].store(ptr, std::memory_order_release);
    index_map.emplace(key, idx);
    return idx;
  }

  /**
   * @brief Interns an Element into the pool or returns an existing instance.
   * @param[in] elem Element struct to intern.
   * @return Immutable pointer to shared Element.
   */
  [[nodiscard]] static const Element *intern(const Element &elem) {
    return &get(internIndex(elem));
  }

  /**
   * @brief Interns an element by symbol and ID into the pool.
   * @param[in] symbol Chemical symbol.
   * @param[in] elem_id Element identifier.
   * @return Immutable pointer to shared Element.
   */
  [[nodiscard]] static const Element *intern(std::string_view symbol, ElementID elem_id) {
    return &get(internIndex(symbol, elem_id));
  }

private:
  static std::array<std::atomic<const Element *>, kMaxElements> &table() noexcept {
    static std::array<std::atomic<const Element *>, kMaxElements> entries{};
    return entries;
  }
};

/**
 * @brief Represents an atom in the simulation cell.
 *
 * Compact 20-byte (float) / 32-byte (double) layout storing Cartesian position, unique ID,
 * and a lock-free 16-bit ElementPool index. Velocities are stored lazily in Cell::velocities_.
 */
class Atom {
public:
  /** @name Constructors */
  ///@{
  explicit Atom() noexcept = default;

  /**
   * @brief Parameterized constructor.
   * @param[in] element The element type of the atom.
   * @param[in] pos The position vector of the atom.
   * @param[in] atom_id The unique ID of the atom.
   */
  explicit Atom(const Element &element, const math::Vector3<real_t> &pos, AtomID atom_id)
      : position_(pos), id_(atom_id), element_pool_idx_(ElementPool::internIndex(element)),
        element_id_(static_cast<std::int16_t>(element.id.value)) {}

  /**
   * @brief Pointer-based flyweight constructor.
   * @param[in] element_ptr Pointer to the interned element.
   * @param[in] pos The position vector of the atom.
   * @param[in] atom_id The unique ID of the atom.
   */
  explicit Atom(const Element *element_ptr, const math::Vector3<real_t> &pos, AtomID atom_id)
      : position_(pos), id_(atom_id),
        element_pool_idx_(element_ptr != nullptr ? ElementPool::internIndex(*element_ptr) : 0),
        element_id_(element_ptr != nullptr ? static_cast<std::int16_t>(element_ptr->id.value)
                                           : static_cast<std::int16_t>(-1)) {}

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
   * @param[in] num The new atom ID.
   */
  void setID(std::uint32_t num) noexcept { id_ = num; }

  /**
   * @brief Gets the position of the atom.
   * @return A const reference to the position vector.
   */
  [[nodiscard]] const math::Vector3<real_t> &position() const noexcept { return position_; }

  /**
   * @brief Sets the position of the atom.
   * @param[in] pos The new position vector.
   */
  void setPosition(const math::Vector3<real_t> &pos) noexcept { position_ = pos; }

  /**
   * @brief Gets the element type of the atom.
   * @return A const reference to the Element struct.
   */
  [[nodiscard]] const Element &element() const noexcept {
    return ElementPool::get(element_pool_idx_);
  }

  /**
   * @brief Sets the element type of the atom.
   * @param[in] ele The new Element struct.
   */
  void setElement(const Element &ele) {
    element_pool_idx_ = ElementPool::internIndex(ele);
    element_id_ = static_cast<std::int16_t>(ele.id.value);
  }

  /**
   * @brief Gets the integer ID of the element type.
   * @return The element ID value.
   */
  [[nodiscard]] std::int32_t elementId() const noexcept {
    return static_cast<std::int32_t>(element_id_);
  }

  ///@}

private:
  math::Vector3<real_t> position_;    ///< Cartesian coordinates in Angstroms.
  AtomID id_{0};                      ///< Unique identification number.
  std::uint16_t element_pool_idx_{0}; ///< Compact index into ElementPool.
  std::int16_t element_id_{-1};       ///< Cached integer element ID (-1 if unassigned).
};

/**
 * @brief Calculates the Euclidean distance between two atoms.
 * @param[in] atom_a The first atom.
 * @param[in] atom_b The second atom.
 * @return The straight-line distance between atoms (not accounting for PBC).
 */
[[nodiscard]] inline real_t distance(const Atom &atom_a, const Atom &atom_b) noexcept {
  return math::norm(atom_a.position() - atom_b.position());
}

/**
 * @brief Calculates the angle (in radians) formed by three atoms.
 * @param[in] center The atom at the vertex of the angle.
 * @param[in] atom_a One of the outer atoms.
 * @param[in] atom_b The other outer atom.
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
