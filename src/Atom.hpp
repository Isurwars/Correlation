#ifndef SRC_ATOM_H_
#define SRC_ATOM_H_
/* ----------------------------------------------------------------------------
 * Correlation: An Analysis Tool for Liquids and for Amorphous Solids
 * Copyright (c) 2013-2025 Isaías Rodríguez <isurwars@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the MIT License version as published in:
 * https://github.com/Isurwars/Correlation/blob/main/LICENSE
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 * ----------------------------------------------------------------------------
 */
#include <array>
#include <string>
#include <vector>

// Minimal structure that represents an atom
struct Atom_Img {
  int element_id;
  int atom_id;
  std::array<double, 3> position;
};

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
   *
   * As well as several methods to calculate distance between pairs of atoms,
   * and angle between terns of atoms
   * --------------------------------------------------------------------------
   */
private:
  int id;
  std::array<double, 3> _position_;
  std::vector<Atom_Img> _bonded_atoms_;
  int _element_id_;
  std::string _element_;
  static int _num_of_atoms_;

public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//

  Atom(std::string, std::array<double, 3>);
  Atom();
  void setAll(std::string, std::array<double, 3>);

  //-------------------------------------------------------------------------//
  //--------------------------- Setters & Getters ---------------------------//
  //-------------------------------------------------------------------------//

  int getID() { return this->id; }
  void setID(int num) { this->id = num; }

  std::array<double, 3> position() { return this->_position_; }
  void setPosition(std::array<double, 3> pos) { this->_position_ = pos; }

  std::vector<Atom_Img> bonded_atoms() { return this->_bonded_atoms_; }
  std::vector<int> getBondedAtomsID();
  void addBondedAtom(Atom_Img);

  int element_id() { return this->_element_id_; }
  void setElementID(int ele_id) { this->_element_id_ = ele_id; }

  std::string element() { return this->_element_; }
  void setElement(const std::string &ele) { this->_element_ = ele; }
  static int getNumberOfAtoms() { return _num_of_atoms_; }

  //-------------------------------------------------------------------------//
  //------------------------------- Methods ---------------------------------//
  //-------------------------------------------------------------------------//

  // Default functions for the object Atom
  double distance(Atom &);
  // Produce a minimal structure to compute the bond angle.
  Atom_Img getImage();
  // get the angle between other atom (Atom_Img) and this object
  double getAngle(Atom_Img, Atom_Img);
};
#endif // SRC_ATOM_H_
