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
double Atom::distance(Atom &other_atom) {
  std::array<double, 3> b_pos = other_atom.position();
  return sqrt(pow(this->_position_[0] - b_pos[0], 2) +
	      pow(this->_position_[1] - b_pos[1], 2) +
	      pow(this->_position_[2] - b_pos[2], 2));
}

// Get Image
Atom_Img Atom::getImage() {
  Atom_Img temp_img;

  temp_img.element_id = this->_element_id_;
  temp_img.atom_id = this->id;
  temp_img.position = this->_position_;
  return temp_img;
}

// Get Bond Angle
double Atom::getAngle(Atom_Img atom_A, Atom_Img atom_B) {
  std::vector<double> vA, vB;
  double vA_, vB_, aux;

  vA = {atom_A.position[0] - this->_position_[0],
	atom_A.position[1] - this->_position_[1],
	atom_A.position[2] - this->_position_[2]};
  vB = {atom_B.position[0] - this->_position_[0],
	atom_B.position[1] - this->_position_[1],
	atom_B.position[2] - this->_position_[2]};
  vA_ = sqrt(vA[0] * vA[0] + vA[1] * vA[1] + vA[2] * vA[2]);
  vB_ = sqrt(vB[0] * vB[0] + vB[1] * vB[1] + vB[2] * vB[2]);
  aux = vA[0] * vB[0] + vA[1] * vB[1] + vA[2] * vB[2];
  aux /= vA_ * vB_;

  if (aux > 1.0)
    aux = 1.0;

  if (aux < -1.0)
    aux = -1.0;
  return acos(aux);
} // Get Bond Angle

void Atom::addBondedAtom(Atom_Img atom_A) {
  this->_bonded_atoms_.push_back(atom_A);
} // Add Bonded Atom

// Get a vector containing the IDs of all bonded atoms
std::vector<int> Atom::getBondedAtomsID() {
  std::vector<Atom_Img>::iterator atom_it;
  // int n_ = this->_bonded_atoms_.size();
  std::vector<int> atoms_ids;
  for (atom_it = this->_bonded_atoms_.begin();
       atom_it != this->_bonded_atoms_.end(); ++atom_it) {
    atoms_ids.push_back(atom_it->atom_id);
  }
  return atoms_ids;
} // getBondedAtomsID
