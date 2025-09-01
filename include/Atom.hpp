#ifndef INCLUDE_ATOM_HPP_
#define INCLUDE_ATOM_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include <string>
#include <vector>

#include "../include/LinearAlgebra.hpp"

class Atom {
  /* --------------------------------------------------------------------------
   * This object represents every atom in the cell.
   *
   * The atributes consist of:
   *
   * ID (Unique Identifier),
   * Element (string and id),
   * Postion (array of three doubles),
   * Bonded Atoms (Images of the bonded atoms)
   * Element ID (id of element in Cell._elements)
   *
   * As well as several methods to calculate distance between pairs of atoms,
   * and angle between terns of atoms
   * --------------------------------------------------------------------------
   */

private:
  int id_;
  linalg::Vector3<double> position_;
  std::vector<Atom> bonded_atoms_;
  std::string element_;
  int element_id_;

public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//
  explicit Atom(std::string, linalg::Vector3<double>, int, int);
  explicit Atom(std::string, linalg::Vector3<double>, int);
  explicit Atom(std::string, linalg::Vector3<double>);
  explicit Atom();

  void resetPositionAndElement(std::string, linalg::Vector3<double>);

  //-------------------------------------------------------------------------//
  //------------------------------- Accessors -------------------------------//
  //-------------------------------------------------------------------------//

  int id() const { return id_; }
  void setID(int num) { id_ = num; }

  const linalg::Vector3<double> &position() const { return position_; }
  void setPosition(linalg::Vector3<double> pos) { position_ = pos; }

  const std::vector<Atom> &bonded_atoms() const { return bonded_atoms_; }
  void setBondedAtoms(const std::vector<Atom> &a) { bonded_atoms_ = a; }
  void addBondedAtom(const Atom &);

  const std::string &element() const { return element_; }
  void setElement(const std::string &ele) { element_ = ele; }

  int element_id() const { return element_id_; }
  void setElementID(int num) { element_id_ = num; }

  //-------------------------------------------------------------------------//
  //------------------------------- Methods ---------------------------------//
  //-------------------------------------------------------------------------//

  // Default functions for the object Atom
  [[nodiscard]] double distance(const Atom &) const;

  // get the angle formed between two other Atoms and this Central Atom
  [[nodiscard]] double angle(const Atom &, const Atom &) const;
};
#endif // INCLUDE_ATOM_HPP_
