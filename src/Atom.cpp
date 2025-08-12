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

Atom::Atom(std::string element, Vector3D pos, int id, int element_id)
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

void Atom::resetPositionAndElement(std::string ele, Vector3D pos) {
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
double Atom::angle(const Atom &atom_A, const Atom &atom_B) const {
  // Use const references to avoid unnecessary copies
  const Vector3D &pos_A = atom_A.position();
  const Vector3D &pos_B = atom_B.position();
  const Vector3D &pos_C = position_;

  // Calculate the vectors representing the bonds
  Vector3D vA = pos_A - pos_C;
  Vector3D vB = pos_B - pos_C;

  // Calculate the magnitudes of the vectors
  double norm_A = norm(vA);
  double norm_B = norm(vB);

  // Handle zero magnitude vectors
  if (norm_A == 0.0 || norm_B == 0.0) {
    return 0.0; // Or throw an exception, depending on desired behavior
  }

  // Calculate the dot product and normalize
  double cos_theta = (vA * vB) / (norm_A * norm_B);

  // Clamp the value to the valid range [-1, 1]
  cos_theta = std::clamp(cos_theta, -1.0, 1.0);

  return std::acos(cos_theta);
} // Get Bond Angle
