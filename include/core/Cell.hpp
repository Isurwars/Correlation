/**
 * @file Cell.hpp
 * @brief Simulation cell structure with periodic boundary conditions.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "core/Atom.hpp"
#include "math/LinearAlgebra.hpp"

#include <array>
#include <optional>
#include <string>
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
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//
  explicit Cell() = default;

  /**
   * @brief Constructs a Cell from three lattice vectors.
   * @param a The first lattice vector.
   * @param b The second lattice vector.
   * @param c The third lattice vector.
   */
  explicit Cell(const math::Vector3<double> &a, const math::Vector3<double> &b,
                const math::Vector3<double> &c);

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
  explicit Cell(const std::array<double, 6> &params);

  /**
   * @brief Move constructor.
   * Transfers ownership of lattice vectors and atom data.
   * @param other Cell object to move from.
   */
  Cell(Cell &&other) noexcept;

  /**
   * @brief Move assignment operator.
   * @param other Cell object to move from.
   * @return Reference to this cell.
   */
  Cell &operator=(Cell &&other) noexcept;

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

  //-------------------------------------------------------------------------//
  //------------------------------- Accessors -------------------------------//
  //-------------------------------------------------------------------------//

  // Lattice Parameters
  /**
   * @brief Gets the lattice parameters (a, b, c, alpha, beta, gamma).
   * @return Array of 6 doubles containing the parameters.
   */
  [[nodiscard]] const std::array<double, 6> &lattice_parameters() const {
    return lattice_parameters_;
  }

  /**
   * @brief Sets the lattice parameters and updates lattice vectors.
   * @param params Array of 6 doubles (a, b, c, alpha, beta, gamma).
   */
  void setLatticeParameters(std::array<double, 6>);

  // Lattice Vectors
  /**
   * @brief Gets the lattice vectors as a 3x3 matrix.
   * @return Constant reference to the lattice vectors matrix.
   */
  [[nodiscard]] const math::Matrix3<double> &latticeVectors() const noexcept {
    return lattice_vectors_;
  }

  /**
   * @brief Gets the inverse lattice vectors as a 3x3 matrix.
   * Useful for converting Cartesian coordinates to fractional coordinates.
   * @return Constant reference to the inverse lattice vectors matrix.
   */
  [[nodiscard]] const math::Matrix3<double> &
  inverseLatticeVectors() const noexcept {
    return inverse_lattice_vectors_;
  }

  // Volume
  /**
   * @brief Gets the volume of the simulation cell.
   * @return The volume value in cubic Angstroms (typically).
   */
  [[nodiscard]] const double &volume() const noexcept { return volume_; }

  // Atoms
  /**
   * @brief Gets the list of atoms in the cell.
   * @return Constant reference to the vector of atoms.
   */
  [[nodiscard]] const std::vector<Atom> &atoms() const noexcept {
    return atoms_;
  }

  // Elements
  /**
   * @brief Gets the list of unique chemical elements in the cell.
   * @return Constant reference to the vector of elements.
   */
  [[nodiscard]] const std::vector<Element> &elements() const {
    return elements_;
  }

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
  [[nodiscard]] std::optional<Element>
  findElement(const std::string &symbol) const;

  //-------------------------------------------------------------------------//
  //-------------------------------- Methods --------------------------------//
  //-------------------------------------------------------------------------//

  /**
   * @brief Applies the minimum image convention to a distance vector.
   *
   * Finds the shortest distance vector between two points under periodic
   * boundary conditions.
   *
   * @param distance The Cartesian distance vector to wrap.
   * @return The minimum image Cartesian distance vector.
   */
  [[nodiscard]] math::Vector3<double>
  minimumImage(const math::Vector3<double> &distance) const;

  /**
   * @brief Adds a new atom to the cell.
   * The atom's element type is automatically registered if it's the first
   * mention of this symbol.
   * @param symbol The chemical element symbol (e.g. "Fe").
   * @param position The Cartesian position [x, y, z] in Angstroms.
   * @return A reference to the newly created Atom.
   */
  Atom &addAtom(const std::string &symbol,
                const math::Vector3<double> &position);

  /**
   * @brief Applies periodic boundary conditions to all atom positions.
   * Wraps all atoms back into the primary simulation cell [0, 1) in fractional
   * coordinates.
   */
  void wrapPositions();

  /**
   * @brief Sets the energy of the cell frame.
   * @param energy The energy value.
   */
  void setEnergy(double energy) { energy_ = energy; }

  /**
   * @brief Gets the energy of the cell frame.
   * @return The energy value.
   */
  [[nodiscard]] double getEnergy() const { return energy_; }

  /**
   * @brief Internal helper to update lattice vectors and recompute volume/inverse.
   * @param new_lattice New 3x3 lattice matrix.
   */
  void updateLattice(const math::Matrix3<double> &new_lattice);

  /**
   * @brief Internal helper to synchronize scalar parameters with vector matrix.
   */
  void updateLatticeParametersFromVectors();

  /**
   * @brief Registers an element symbol and returns its ID.
   * @param symbol Element symbol (e.g. "O").
   * @return The existing or newly assigned ElementID.
   */
  ElementID getOrRegisterElement(const std::string &symbol);

  math::Matrix3<double> lattice_vectors_; ///< Basis vectors of the box.
  math::Matrix3<double> inverse_lattice_vectors_; ///< Inverse matrix for fractional mapping.
  std::array<double, 6> lattice_parameters_; ///< {a, b, c, alpha, beta, gamma}.
  double volume_{0.0}; ///< Cached volume in Angstroms^3.
  double energy_{0.0}; ///< Potential energy of this specific coordinate set.
  std::vector<Atom> atoms_; ///< Collection of atoms in the cell.
  std::vector<Element> elements_; ///< Unique elements present in the system.
};

} // namespace correlation::core
