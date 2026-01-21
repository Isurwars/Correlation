// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/Cell.hpp"

#include <algorithm>
#include <cmath>

#include "../include/Atom.hpp"
#include "../include/LinearAlgebra.hpp"
#include "../include/PhysicalData.hpp"

//---------------------------------------------------------------------------//
//------------------------------- Constructors ------------------------------//
//---------------------------------------------------------------------------//
Cell::Cell(const linalg::Vector3<double> &a, const linalg::Vector3<double> &b,
           const linalg::Vector3<double> &c) {
  updateLattice(linalg::Matrix3<double>(a, b, c));
}

Cell::Cell(const std::array<double, 6> &params) {
  setLatticeParameters(params);
}

// Move Constructor
Cell::Cell(Cell &&other) noexcept
    : lattice_vectors_(std::move(other.lattice_vectors_)),
      inverse_lattice_vectors_(std::move(other.inverse_lattice_vectors_)),
      lattice_parameters_(std::move(other.lattice_parameters_)),
      volume_(std::move(other.volume_)), atoms_(std::move(other.atoms_)),
      elements_(std::move(other.elements_)) {}

// Move Assignment Operator
Cell &Cell::operator=(Cell &&other) noexcept {
  if (this != &other) {
    lattice_vectors_ = std::move(other.lattice_vectors_);
    inverse_lattice_vectors_ = std::move(other.inverse_lattice_vectors_);
    lattice_parameters_ = std::move(other.lattice_parameters_);
    volume_ = std::move(other.volume_);
    atoms_ = std::move(other.atoms_);
    elements_ = std::move(other.elements_);
  }
  return *this;
}

//---------------------------------------------------------------------------//
//-------------------------------- Accessors --------------------------------//
//---------------------------------------------------------------------------//

void Cell::setLatticeParameters(std::array<double, 6> params) {
  lattice_parameters_ = params;
  const double a = params[0], b = params[1], c = params[2];
  const double alpha = params[3] * constants::deg2rad;
  const double beta = params[4] * constants::deg2rad;
  const double gamma = params[5] * constants::deg2rad;

  if (a <= 0 || b <= 0 || c <= 0) {
    throw std::invalid_argument("Lattice parameters a, b, c must be positive.");
  }

  const double cos_g = std::cos(gamma);
  const double sin_g = std::sin(gamma);

  // Standard conversion from lattice parameters to vectors
  linalg::Vector3<double> v_a = {a, 0.0, 0.0};
  linalg::Vector3<double> v_b = {b * cos_g, b * sin_g, 0.0};
  linalg::Vector3<double> v_c = {
      c * std::cos(beta),
      c * (std::cos(alpha) - std::cos(beta) * cos_g) / sin_g,
      0.0 // z-component is calculated from volume
  };
  // Re-calculate z-component of v_c from volume to ensure orthogonality
  const double volume =
      a * b * c *
      std::sqrt(1.0 - std::pow(std::cos(alpha), 2) -
                std::pow(std::cos(beta), 2) - std::pow(std::cos(gamma), 2) +
                2.0 * std::cos(alpha) * std::cos(beta) * std::cos(gamma));
  v_c.z() = volume / (a * b * sin_g);
  updateLattice(linalg::Matrix3<double>(v_a, v_b, v_c));
}

void Cell::updateLattice(const linalg::Matrix3<double> &new_lattice) {
  lattice_vectors_ = new_lattice;
  volume_ = linalg::determinant(lattice_vectors_);
  if (volume_ <= 1e-9) {
    throw std::logic_error("Cell volume must be positive.");
  }
  inverse_lattice_vectors_ =
      linalg::transpose(linalg::invert(lattice_vectors_));
  updateLatticeParametersFromVectors();
}

void Cell::updateLatticeParametersFromVectors() {
  const auto &a_vec = lattice_vectors_[0];
  const auto &b_vec = lattice_vectors_[1];
  const auto &c_vec = lattice_vectors_[2];

  const double a = linalg::norm(a_vec);
  const double b = linalg::norm(b_vec);
  const double c = linalg::norm(c_vec);

  if (a < 1e-9 || b < 1e-9 || c < 1e-9) {
    // Handle case of zero-length vectors, though updateLattice would likely
    // throw first.
    lattice_parameters_ = {0, 0, 0, 0, 0, 0};
    return;
  }

  const double alpha_rad = std::acos(linalg::dot(b_vec, c_vec) / (b * c));
  const double beta_rad = std::acos(linalg::dot(a_vec, c_vec) / (a * c));
  const double gamma_rad = std::acos(linalg::dot(a_vec, b_vec) / (a * b));

  lattice_parameters_ = {a,
                         b,
                         c,
                         alpha_rad * constants::rad2deg,
                         beta_rad * constants::rad2deg,
                         gamma_rad * constants::rad2deg};
}

//---------------------------------------------------------------------------//
//--------------------------------- Methods ---------------------------------//
//---------------------------------------------------------------------------//

void Cell::precomputeBondCutoffs() {
  const size_t num_elements = elements_.size();
  bond_cutoffs_sq_.resize(num_elements, std::vector<double>(num_elements));

  for (size_t i = 0; i < num_elements; ++i) {
    const double radius_A = CovalentRadii::get(elements_[i].symbol);
    for (size_t j = i; j < num_elements; ++j) {
      const double radius_B = CovalentRadii::get(elements_[j].symbol);
      const double max_bond_dist = (radius_A + radius_B) * 1.2;
      const double max_bond_dist_sq = max_bond_dist * max_bond_dist;
      bond_cutoffs_sq_[i][j] = max_bond_dist_sq;
      bond_cutoffs_sq_[j][i] = max_bond_dist_sq;
    }
  }
}

std::optional<Element> Cell::findElement(const std::string &symbol) const {
  auto it = std::find_if(elements_.begin(), elements_.end(),
                         [&](const Element &e) { return e.symbol == symbol; });
  if (it != elements_.end()) {
    return *it;
  }
  return std::nullopt;
}

double Cell::getBondCutoff(int Type_A, int Type_B) {
  if (bond_cutoffs_sq_.empty())
    precomputeBondCutoffs();
  return sqrt(bond_cutoffs_sq_[Type_A][Type_B]);
}

ElementID Cell::getOrRegisterElement(const std::string &symbol) {
  if (auto existing_element = findElement(symbol)) {
    return existing_element->id;
  }
  // Register the new element
  ElementID new_id{static_cast<int>(elements_.size())};
  elements_.push_back({symbol, new_id});
  return new_id;
}

Atom &Cell::addAtom(const std::string &symbol,
                    const linalg::Vector3<double> &position) {
  ElementID element_id = getOrRegisterElement(symbol);
  auto element_it =
      std::find_if(elements_.begin(), elements_.end(), [&](const Element &e) {
        return e.id.value == element_id.value;
      });

  AtomID new_atom_id{static_cast<std::uint32_t>(atoms_.size())};
  atoms_.emplace_back(*element_it, position, new_atom_id);
  return atoms_.back();
}

void Cell::wrapPositions() {
  for (Atom &atom : atoms_) {
    linalg::Vector3<double> frac_pos =
        inverse_lattice_vectors_ * atom.position();
    frac_pos.x() -= std::floor(frac_pos.x());
    frac_pos.y() -= std::floor(frac_pos.y());
    frac_pos.z() -= std::floor(frac_pos.z());
    atom.setPosition(lattice_vectors_ * frac_pos);
  }
}
