/**
 * @file Cell.hpp
 * @brief Simulation cell structure with periodic boundary conditions.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Atom.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"

#include <array>
#include <optional>
#include <string_view>
#include <vector>

namespace correlation::core {

/**
 * @brief Represents the simulation cell (periodic box).
 *
 * This class handles the lattice vectors, periodic boundary conditions, and
 * stores the atoms contained within the cell.
 */
class Cell {

public:
  /** @name Constructors */
  ///@{
  explicit Cell() = default;

  /**
   * @brief Constructs a Cell from three lattice vectors.
   * @param vec_a The first lattice vector.
   * @param vec_b The second lattice vector.
   * @param vec_c The third lattice vector.
   */
  explicit Cell(const math::Vector3<real_t> &vec_a, const math::Vector3<real_t> &vec_b,
                const math::Vector3<real_t> &vec_c);

  /**
   * @brief Constructs a Cell from lattice parameters {a, b, c, alpha, beta,
   * gamma}.
   * @param params An array containing the six lattice parameters:
   *               - params[0]: a (length of vector a in Angstroms)
   *               - params[1]: b (length of vector b in Angstroms)
   *               - params[2]: c (length of vector c in Angstroms)
   *               - params[3]: alpha (angle between b and c in degrees)
   *               - params[4]: beta (angle between a and c in degrees)
   *               - params[5]: gamma (angle between a and b in degrees)
   */
  explicit Cell(const std::array<real_t, 6> &params);

  /**
   * @brief Move constructor.
   * Transfers ownership of lattice vectors and atom data.
   * @param other Cell object to move from.
   */
  Cell(Cell &&other) noexcept = default;

  /**
   * @brief Move assignment operator.
   * @param other Cell object to move from.
   * @return Reference to this cell.
   */
  Cell &operator=(Cell &&other) noexcept = default;

  /**
   * @brief Copy constructor.
   * @param other Cell object to copy from.
   */
  Cell(const Cell &other) = default;

  /**
   * @brief Copy assignment operator.
   * @param other Cell object to copy from.
   * @return Reference to this cell.
   */
  Cell &operator=(const Cell &other) = default;

  /**
   * @brief Destructor.
   */
  ~Cell() = default;

  ///@}

  /** @name Accessors */
  ///@{

  // Lattice Parameters
  /**
   * @brief Gets the lattice parameters (a, b, c, alpha, beta, gamma).
   * @return Array of 6 real_t containing the parameters.
   */
  [[nodiscard]] const std::array<real_t, 6> &latticeParameters() const noexcept {
    return lattice_parameters_;
  }

  void setLatticeParameters(const std::array<real_t, 6> &params);

  // Lattice Vectors
  /**
   * @brief Gets the lattice vectors as a 3x3 matrix.
   * @return Constant reference to the lattice vectors matrix.
   */
  [[nodiscard]] const math::Matrix3<real_t> &latticeVectors() const noexcept {
    return lattice_vectors_;
  }

  /**
   * @brief Gets the inverse lattice vectors as a 3x3 matrix.
   * Useful for converting Cartesian coordinates to fractional coordinates.
   * @return Constant reference to the inverse lattice vectors matrix.
   */
  [[nodiscard]] const math::Matrix3<real_t> &inverseLatticeVectors() const noexcept {
    return inverse_lattice_vectors_;
  }

  // Volume
  /**
   * @brief Gets the volume of the simulation cell.
   * @return The volume value in cubic Angstroms (typically).
   */
  [[nodiscard]] const real_t &volume() const noexcept { return volume_; }

  // Atoms
  /**
   * @brief Gets the list of atoms in the cell.
   * @return Constant reference to the vector of atoms.
   */
  [[nodiscard]] const std::vector<Atom> &atoms() const noexcept { return atoms_; }

  /**
   * @brief Gets a mutable list of atoms in the cell.
   * @return Reference to the vector of atoms.
   */
  [[nodiscard]] std::vector<Atom> &atoms() noexcept { return atoms_; }

  // Elements
  /**
   * @brief Gets the list of unique chemical elements in the cell.
   * @return Constant reference to the vector of elements.
   */
  [[nodiscard]] const std::vector<Element> &elements() const noexcept { return elements_; }

  /**
   * @brief Gets the total number of atoms in the cell.
   * @return The number of atoms.
   */
  [[nodiscard]] size_t atomCount() const noexcept { return atoms_.size(); }

  /**
   * @brief Checks whether the cell contains no atoms.
   * @return True if empty, false otherwise.
   */
  [[nodiscard]] bool isEmpty() const noexcept { return atoms_.empty(); }

  /**
   * @brief Finds the Element properties for a given element symbol.
   * @param symbol The element symbol (e.g., "Si").
   * @return An optional containing the Element struct if found, otherwise
   * std::nullopt.
   */
  [[nodiscard]] std::optional<Element> findElement(std::string_view symbol) const;

  ///@}

  /** @name Methods */
  ///@{

  /**
   * @brief Applies the minimum image convention to a distance vector.
   *
   * Finds the shortest distance vector between two points under periodic
   * boundary conditions.
   *
   * @param distance The Cartesian distance vector to wrap.
   * @return The minimum image Cartesian distance vector.
   */
  [[nodiscard]] math::Vector3<real_t> minimumImage(const math::Vector3<real_t> &distance) const;

  /**
   * @brief Adds a new atom to the cell.
   * The atom's element type is automatically registered if it's the first
   * mention of this symbol.
   * @param symbol The chemical element symbol (e.g. "Fe").
   * @param position The Cartesian position [x, y, z] in Angstroms.
   * @return A reference to the newly created Atom.
   */
  Atom &addAtom(std::string_view symbol, const math::Vector3<real_t> &position);

  /**
   * @brief Applies periodic boundary conditions to all atom positions.
   * Wraps all atoms back into the primary simulation cell [0, 1) in fractional
   * coordinates.
   */
  void wrapPositions();

  /**
   * @brief Computes the perpendicular widths (interplanar spacing) between opposing cell faces.
   *
   * For lattice vectors a, b, c with volume V, returns {V / ||b x c||, V / ||a x c||, V / ||a x
   * b||}. If volume is non-positive or vectors are degenerate, returns {0, 0, 0}.
   *
   * @return Array containing {w_a, w_b, w_c} in Angstroms.
   */
  [[nodiscard]] std::array<real_t, 3> perpendicularWidths() const noexcept;

  /**
   * @brief Replicates the simulation cell na x nb x nc times along lattice vectors.
   *
   * @param rep_a Replication factor along lattice vector a (>= 1).
   * @param rep_b Replication factor along lattice vector b (>= 1).
   * @param rep_c Replication factor along lattice vector c (>= 1).
   * @return Replicated supercell with expanded lattice vectors and translated atoms.
   * @throws std::invalid_argument If any factor is < 1 or total atom count exceeds 50,000.
   */
  [[nodiscard]] Cell replicate(int rep_a, int rep_b, int rep_c) const;

  /**
   * @brief Automatically expands the cell to satisfy the Minimum Image Convention for cutoff r_cut.
   *
   * Computes minimal replication factors along each axis such that every perpendicular width
   * satisfies w_k >= 2 * r_cut, up to safety limits.
   *
   * @param r_cut Interaction or evaluation cutoff radius in Angstroms.
   * @param max_replication Maximum allowable replication factor per axis (default: 10).
   * @param max_radius Maximum allowable cutoff radius ceiling (default: 50.0 Angstroms).
   * @return A minimally replicated supercell satisfying w_perp >= 2 * r_cut, or *this if already
   * large enough.
   * @throws std::invalid_argument If r_cut > max_radius or required replication exceeds
   * max_replication.
   */
  [[nodiscard]] Cell autoSupercell(real_t r_cut, int max_replication = 10,
                                   real_t max_radius = correlation::math::MAX_CUTOFF_RADIUS) const;

  /**
   * @brief Sets the energy of the cell frame.
   * @param energy The energy value.
   */
  void setEnergy(real_t energy) noexcept { energy_ = energy; }

  /**
   * @brief Gets the energy of the cell frame.
   * @return The energy value.
   */
  [[nodiscard]] real_t getEnergy() const noexcept { return energy_; }

  /**
   * @brief Reserves capacity for atom storage to minimize reallocations.
   * @param count Expected number of atoms.
   */
  void reserveAtoms(std::size_t count) { atoms_.reserve(count); }

  /**
   * @brief Updates lattice vectors and recomputes volume, inverse matrix, and scalar parameters.
   * @param new_lattice New 3x3 lattice matrix.
   */
  void updateLattice(const math::Matrix3<real_t> &new_lattice);

  ///@}

private:
  /**
   * @brief Internal helper to synchronize scalar parameters with vector matrix.
   */
  void updateLatticeParametersFromVectors();

  /**
   * @brief Registers an element symbol if not already present and returns its Element reference.
   * @param symbol Element symbol (e.g. "O").
   * @return Reference to the registered Element.
   */
  const Element &getOrRegisterElement(std::string_view symbol);

  math::Matrix3<real_t> lattice_vectors_;         ///< Basis vectors of the box.
  math::Matrix3<real_t> inverse_lattice_vectors_; ///< Inverse matrix for fractional mapping.
  std::array<real_t, 6> lattice_parameters_{};    ///< {a, b, c, alpha, beta, gamma}.
  real_t volume_{0.0};                            ///< Cached volume in Angstroms^3.
  real_t energy_{0.0};            ///< Potential energy of this specific coordinate set.
  std::vector<Atom> atoms_;       ///< Collection of atoms in the cell.
  std::vector<Element> elements_; ///< Unique elements present in the system.
};

} // namespace correlation::core
