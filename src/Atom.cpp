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
#include "Atom.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <numeric>

int Atom::_num_of_atoms_ = 0;

//---------------------------------------------------------------------------//
//------------------------------ Constructors -------------------------------//
//---------------------------------------------------------------------------//
void Atom::setAll(std::string ele, std::array<double, 3> pos) {
  this->_element_ = ele;
  this->_position_ = pos;
}
// Complete constructor for the Atom object
Atom::Atom(std::string ele, std::array<double, 3> pos) {
  this->_element_ = ele;
  this->_position_ = pos;
  this->id = Atom::_num_of_atoms_;
  Atom::_num_of_atoms_++;
}
// Default constructor for the Atom object
Atom::Atom() {
  this->_element_ = "H";
  this->_position_ = {0.0, 0.0, 0.0};
  this->id = Atom::_num_of_atoms_;
  Atom::_num_of_atoms_++;
}

//----------------------------------------------------------------------------//
//--------------------------------- Methods ----------------------------------//
//----------------------------------------------------------------------------//
// Distance to the other atom
double Atom::distance(const Atom &other_atom) const {
  return std::hypot(this->position_[0] - other_atom.position_[0],
                   this->position_[1] - other_atom.position_[1],
                   this->position_[2] - other_atom.position_[2]);
}

// Get Image
void Atom::getImage(Atom_Img& img) const {
  img.element_id = this->element_id_;
  img.atom_id = this->id;
  img.position = this->position_;
}

// Get Bond Angle
double Atom::getAngle(const Atom_Img& atom_A, const Atom_Img& atom_B) const {
  // Use const references to avoid unnecessary copies
  const std::vector<double>& pos_A = atom_A.position;
  const std::vector<double>& pos_B = atom_B.position;
  const std::vector<double>& pos_C = this->_position_;

  // Calculate the vectors representing the bonds
  std::vector<double> vA = {pos_A[0] - pos_C[0], pos_A[1] - pos_C[1], pos_A[2] - pos_C[2]};
  std::vector<double> vB = {pos_B[0] - pos_C[0], pos_B[1] - pos_C[1], pos_B[2] - pos_C[2]};

  // Calculate the magnitudes of the vectors
  double vA_ = std::sqrt(vA[0] * vA[0] + vA[1] * vA[1] + vA[2] * vA[2]);
  double vB_ = std::sqrt(vB[0] * vB[0] + vB[1] * vB[1] + vB[2] * vB[2]);

  // Handle zero magnitude vectors
  if (vA_ == 0.0 || vB_ == 0.0) {
    return 0.0; // Or throw an exception, depending on desired behavior
  }

  // Calculate the dot product and normalize
  double aux = (vA[0] * vB[0] + vA[1] * vB[1] + vA[2] * vB[2]) / (vA_ * vB_);

  // Clamp the value to the valid range [-1, 1]
  aux = std::clamp(aux, -1.0, 1.0);

  return std::acos(aux);
} // Get Bond Angle

void Atom::addBondedAtom(const Atom_Img& atom_A) {
  try {
    this->bonded_atoms_.push_back(atom_A);
  } catch (const std::bad_alloc& e) {
    throw std::runtime_error("Memory allocation failed in addBondedAtom: " + std::string(e.what()));
  }
} // addBondedAtom

// Get a vector containing the IDs of all bonded atoms
std::vector<int> Atom::getBondedAtomsID() const {
  std::vector<int> atoms_ids(this->bonded_atoms_.size());
  std::transform(this->bonded_atoms_.begin(), this->bonded_atoms_.end(), atoms_ids.begin(),
                 [](const Atom_Img& atom_img) { return atom_img.atom_id; });
  return atoms_ids;
} // getBondedAtomsID
