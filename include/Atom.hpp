#ifndef INCLUDE_ATOM_HPP_
#define INCLUDE_ATOM_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include <array>
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
  Vector3D position_;
  std::vector<Atom> bonded_atoms_;
  std::string element_;
  int element_id_;

public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//
  explicit Atom(std::string, std::array<double, 3>, int, int);
  explicit Atom(std::string, std::array<double, 3>, int);
  explicit Atom(std::string, std::array<double, 3>);
  explicit Atom();

  void setAll(std::string, Vector3D);

  //-------------------------------------------------------------------------//
  //------------------------------- Accessors -------------------------------//
  //-------------------------------------------------------------------------//

  int id() const { return id_; }
  void setID(int num) { id_ = num; }

  const Vector3D position() const { return position_; }
  void setPosition(Vector3D pos) { position_ = pos; }

  const std::vector<Atom> bonded_atoms() const { return bonded_atoms_; }
  void setBondedAtom( const std::vector<Atom> &a) {bonded_atoms_ = a;}
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
  [[nodiscard]] double getAngle(const Atom &, const Atom &) const;
};
#endif // INCLUDE_ATOM_HPP_
