/**
 * @file Cell.cpp
 * @brief Implementation of the simulation cell and periodic boundary logic.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "core/Cell.hpp"
#include "core/Atom.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"

#include <algorithm>
#include <cmath>

namespace correlation::core {

//---------------------------------------------------------------------------//
//------------------------------- Constructors ------------------------------//
//---------------------------------------------------------------------------//
Cell::Cell(const math::Vector3<double> &a, const math::Vector3<double> &b,
           const math::Vector3<double> &c) {
  updateLattice(math::Matrix3<double>(a, b, c));
}

Cell::Cell(const std::array<double, 6> &params) {
  setLatticeParameters(params);
}

// Move Constructor
Cell::Cell(Cell &&other) noexcept
    : lattice_vectors_(std::move(other.lattice_vectors_)),
      inverse_lattice_vectors_(std::move(other.inverse_lattice_vectors_)),
      lattice_parameters_(std::move(other.lattice_parameters_)),
      volume_(std::move(other.volume_)), energy_(std::move(other.energy_)),
      atoms_(std::move(other.atoms_)), elements_(std::move(other.elements_)) {}

// Move Assignment Operator
Cell &Cell::operator=(Cell &&other) noexcept {
  if (this != &other) {
    lattice_vectors_ = std::move(other.lattice_vectors_);
    inverse_lattice_vectors_ = std::move(other.inverse_lattice_vectors_);
    lattice_parameters_ = std::move(other.lattice_parameters_);
    volume_ = std::move(other.volume_);
    energy_ = std::move(other.energy_);
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
  const double alpha = params[3] * math::deg_to_rad;
  const double beta = params[4] * math::deg_to_rad;
  const double gamma = params[5] * math::deg_to_rad;

  if (a <= 0 || b <= 0 || c <= 0) {
    throw std::invalid_argument("Lattice parameters a, b, c must be positive.");
  }

  const double cos_g = std::cos(gamma);
  const double sin_g = std::sin(gamma);

  // Standard conversion from lattice parameters (lengths and angles) to lattice
  // vectors. We align 'a' with the x-axis, and 'b' in the xy-plane.
  math::Vector3<double> v_a = {a, 0.0, 0.0};
  math::Vector3<double> v_b = {b * cos_g, b * sin_g, 0.0};
  math::Vector3<double> v_c = {
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
  updateLattice(math::Matrix3<double>(v_a, v_b, v_c));
}

void Cell::updateLattice(const math::Matrix3<double> &new_lattice) {
  lattice_vectors_ = new_lattice;
  volume_ = math::determinant(lattice_vectors_);
  if (volume_ <= 1e-9) {
    throw std::logic_error("Cell volume must be positive.");
  }
  inverse_lattice_vectors_ = math::transpose(math::invert(lattice_vectors_));
  updateLatticeParametersFromVectors();
}

void Cell::updateLatticeParametersFromVectors() {
  const auto &a_vec = lattice_vectors_[0];
  const auto &b_vec = lattice_vectors_[1];
  const auto &c_vec = lattice_vectors_[2];

  const double a = math::norm(a_vec);
  const double b = math::norm(b_vec);
  const double c = math::norm(c_vec);

  if (a < 1e-9 || b < 1e-9 || c < 1e-9) {
    // Handle case of zero-length vectors, though updateLattice would likely
    // throw first.
    lattice_parameters_ = {0, 0, 0, 0, 0, 0};
    return;
  }

  const double alpha_rad = std::acos(math::dot(b_vec, c_vec) / (b * c));
  const double beta_rad = std::acos(math::dot(a_vec, c_vec) / (a * c));
  const double gamma_rad = std::acos(math::dot(a_vec, b_vec) / (a * b));

  lattice_parameters_ = {a,
                         b,
                         c,
                         alpha_rad * math::rad_to_deg,
                         beta_rad * math::rad_to_deg,
                         gamma_rad * math::rad_to_deg};
}

//---------------------------------------------------------------------------//
//--------------------------------- Methods ---------------------------------//
//---------------------------------------------------------------------------//

std::optional<Element> Cell::findElement(const std::string &symbol) const {
  auto it = std::find_if(elements_.begin(), elements_.end(),
                         [&](const Element &e) { return e.symbol == symbol; });
  if (it != elements_.end()) {
    return *it;
  }
  return std::nullopt;
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
                    const math::Vector3<double> &position) {
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
    math::Vector3<double> frac_pos = inverse_lattice_vectors_ * atom.position();
    frac_pos.x() -= std::floor(frac_pos.x());
    frac_pos.y() -= std::floor(frac_pos.y());
    frac_pos.z() -= std::floor(frac_pos.z());
    atom.setPosition(lattice_vectors_ * frac_pos);
  }
}

} // namespace correlation::core
