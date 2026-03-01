// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <array>
#include <optional>
#include <string>
#include <vector>

#include "Atom.hpp"
#include "LinearAlgebra.hpp"

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
  explicit Cell(const linalg::Vector3<double> &a,
                const linalg::Vector3<double> &b,
                const linalg::Vector3<double> &c);

  /**
   * @brief Constructs a Cell from lattice parameters {a, b, c, alpha, beta,
   * gamma}.
   * @param params An array containing the six lattice parameters.
   */
  explicit Cell(const std::array<double, 6> &params);

  // Move constructor
  Cell(Cell &&other) noexcept;

  // Move assignment operator
  Cell &operator=(Cell &&other) noexcept;

  // The copy constructor and assignment operator can be explicitly defaulted,
  // or implicitly generated if no other special member functions are declared.
  Cell(const Cell &) = default;
  Cell &operator=(const Cell &) = default;

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
  [[nodiscard]] const linalg::Matrix3<double> &latticeVectors() const noexcept {
    return lattice_vectors_;
  }

  /**
   * @brief Gets the inverse lattice vectors as a 3x3 matrix.
   * Useful for converting Cartesian coordinates to fractional coordinates.
   * @return Constant reference to the inverse lattice vectors matrix.
   */
  [[nodiscard]] const linalg::Matrix3<double> &
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
   * @brief Adds a new atom to the cell.
   * The atom's element type is automatically registered.
   * @param symbol The element symbol of the atom.
   * @param position The Cartesian position of the atom.
   * @return The newly created Atom object.
   */
  Atom &addAtom(const std::string &symbol,
                const linalg::Vector3<double> &position);

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

private:
  void updateLattice(const linalg::Matrix3<double> &new_lattice);
  void updateLatticeParametersFromVectors();
  ElementID getOrRegisterElement(const std::string &symbol);

  linalg::Matrix3<double> lattice_vectors_;
  linalg::Matrix3<double> inverse_lattice_vectors_;
  std::array<double, 6> lattice_parameters_;
  double volume_{0.0};
  double energy_{0.0};
  std::vector<Atom> atoms_;
  std::vector<Element> elements_;
};
