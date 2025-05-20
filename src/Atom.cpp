// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/Atom.hpp"

#include <algorithm>
#include <cmath>

//---------------------------------------------------------------------------//
//------------------------------ Constructors -------------------------------//
//---------------------------------------------------------------------------//

Atom::Atom(std::string element, Vector3D pos, int id,
           int element_id)
    : id_{id}, position_{pos}, element_{std::move(element)},
      element_id_{element_id} {}

Atom::Atom(std::string element, Vector3D pos, int id)
    : Atom(std::move(element), pos, id, -1) {}

Atom::Atom(std::string element, Vector3D pos)
    : Atom(std::move(element), pos, -1, -1) {}

Atom::Atom() : Atom("H", {0.0, 0.0, 0.0}, -1, -1) {}

//---------------------------------------------------------------------------//
//-------------------------------- Accessors --------------------------------//
//---------------------------------------------------------------------------//

void Atom::addBondedAtom(const Atom &atom_A) {
    bonded_atoms_.push_back(atom_A);
} // addBondedAtom

void Atom::setAll(std::string ele, Vector3D pos) {
  element_ = ele;
  position_ = pos;
}
//----------------------------------------------------------------------------//
//--------------------------------- Methods ----------------------------------//
//----------------------------------------------------------------------------//
// Distance to the other atom
double Atom::distance(const Atom &other_atom) const {
  return std::hypot(position_[0] - other_atom.position_[0],
                    position_[1] - other_atom.position_[1],
                    position_[2] - other_atom.position_[2]);
}

// Get Bond Angle
double Atom::getAngle(const Atom &atom_A, const Atom &atom_B) const {
  // Use const references to avoid unnecessary copies
  const Vector3D pos_A = atom_A.position();
  const Vector3D pos_B = atom_B.position();
  const Vector3D pos_C = position_;

  // Calculate the vectors representing the bonds
  Vector3D vA = {pos_A[0] - pos_C[0], pos_A[1] - pos_C[1], pos_A[2] - pos_C[2]};
  Vector3D vB = {pos_B[0] - pos_C[0], pos_B[1] - pos_C[1], pos_B[2] - pos_C[2]};

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
