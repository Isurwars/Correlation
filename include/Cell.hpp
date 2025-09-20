// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <array>
#include <optional>
#include <string>
#include <vector>

#include "Atom.hpp"
#include "LinearAlgebra.hpp"

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
  const std::array<double, 6> &lattice_parameters() const {
    return lattice_parameters_;
  }
  void setLatticeParameters(std::array<double, 6>);
  // Lattice Vectors
  const linalg::Matrix3<double> &latticeVectors() const noexcept {
    return lattice_vectors_;
  }
  const linalg::Matrix3<double> &inverseLatticeVectors() const noexcept {
    return inverse_lattice_vectors_;
  }
  // Volume
  const double &volume() const noexcept { return volume_; }
  // Atoms
  const std::vector<Atom> &atoms() const noexcept { return atoms_; }
  // Elements
  const std::vector<Element> &elements() const { return elements_; }

  size_t atomCount() const noexcept { return atoms_.size(); }
  bool isEmpty() const noexcept { return atoms_.empty(); }

  /**
   * @brief Finds the Element properties for a given element symbol.
   * @param symbol The element symbol (e.g., "Si").
   * @return An optional containing the Element struct if found, otherwise
   * std::nullopt.
   */
  std::optional<Element> findElement(const std::string &symbol) const;

  // --- Mutating Commands ---

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

private:
  void updateLattice(const linalg::Matrix3<double> &new_lattice);
  void updateLatticeParametersFromVectors();
  ElementID getOrRegisterElement(const std::string &symbol);

  linalg::Matrix3<double> lattice_vectors_;
  linalg::Matrix3<double> inverse_lattice_vectors_;
  std::array<double, 6> lattice_parameters_;
  double volume_{0.0};

  std::vector<Atom> atoms_;
  std::vector<Element> elements_;
};
