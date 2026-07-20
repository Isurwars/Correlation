/**
 * @file Cell.cpp
 * @brief Implementation of the simulation cell and periodic boundary logic.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/Cell.hpp"
#include "core/Atom.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"

#include <algorithm>
#include <cmath>

namespace correlation::core {

Cell::Cell(const math::Vector3<real_t> &vec_a, const math::Vector3<real_t> &vec_b, const math::Vector3<real_t> &vec_c) {
  updateLattice(math::Matrix3<real_t>(vec_a, vec_b, vec_c));
}

Cell::Cell(const std::array<real_t, 6> &params) { setLatticeParameters(params); }

void Cell::setLatticeParameters(std::array<real_t, 6> params) {
  lattice_parameters_ = params;
  const real_t len_a = params[0];
  const real_t len_b = params[1];
  const real_t len_c = params[2];
  const real_t alpha = params[3] * static_cast<real_t>(math::deg_to_rad);
  const real_t beta = params[4] * static_cast<real_t>(math::deg_to_rad);
  const real_t gamma = params[5] * static_cast<real_t>(math::deg_to_rad);

  if (std::isnan(len_a) || std::isnan(len_b) || std::isnan(len_c) || len_a <= 0 || len_b <= 0 || len_c <= 0) {
    throw std::invalid_argument("Lattice parameters a, b, c must be positive.");
  }

  if (std::isnan(params[3]) || std::isnan(params[4]) || std::isnan(params[5]) || params[3] <= 0 || params[3] >= 180 ||
      params[4] <= 0 || params[4] >= 180 || params[5] <= 0 || params[5] >= 180) {
    throw std::invalid_argument("Lattice angles must be between 0 and 180 degrees.");
  }

  const real_t cos_g = std::cos(gamma);
  const real_t sin_g = std::sin(gamma);

  math::Vector3<real_t> const v_a = {len_a, 0.0, 0.0};
  math::Vector3<real_t> const v_b = {len_b * cos_g, len_b * sin_g, 0.0};
  math::Vector3<real_t> v_c = {len_c * std::cos(beta), len_c * (std::cos(alpha) - std::cos(beta) * cos_g) / sin_g, 0.0};
  const auto volume = static_cast<real_t>(len_a * len_b * len_c *
                                          std::sqrt(1.0 - std::pow(std::cos(alpha), 2) - std::pow(std::cos(beta), 2) -
                                                    std::pow(std::cos(gamma), 2) +
                                                    2.0 * std::cos(alpha) * std::cos(beta) * std::cos(gamma)));
  v_c.z() = volume / (len_a * len_b * sin_g);
  updateLattice(math::Matrix3<real_t>(v_a, v_b, v_c));
}

void Cell::updateLattice(const math::Matrix3<real_t> &new_lattice) {
  lattice_vectors_ = new_lattice;
  volume_ = math::determinant(lattice_vectors_);
  if (std::isnan(volume_) || volume_ <= 1e-9) {
    throw std::logic_error("Cell volume must be positive and finite.");
  }
  inverse_lattice_vectors_ = math::invert(lattice_vectors_);
  updateLatticeParametersFromVectors();
}

void Cell::updateLatticeParametersFromVectors() {
  const auto &a_vec = lattice_vectors_[0];
  const auto &b_vec = lattice_vectors_[1];
  const auto &c_vec = lattice_vectors_[2];

  const real_t len_a = math::norm(a_vec);
  const real_t len_b = math::norm(b_vec);
  const real_t len_c = math::norm(c_vec);

  if (std::isnan(len_a) || std::isnan(len_b) || std::isnan(len_c) || len_a < 1e-9 || len_b < 1e-9 || len_c < 1e-9) {
    lattice_parameters_ = {0, 0, 0, 0, 0, 0};
    return;
  }

  const real_t cos_alpha = math::dot(b_vec, c_vec) / (len_b * len_c);
  const real_t cos_beta = math::dot(a_vec, c_vec) / (len_a * len_c);
  const real_t cos_gamma = math::dot(a_vec, b_vec) / (len_a * len_b);

  if (std::isnan(cos_alpha) || std::isnan(cos_beta) || std::isnan(cos_gamma)) {
    lattice_parameters_ = {0, 0, 0, 0, 0, 0};
    return;
  }

  const real_t alpha_rad = std::acos(std::clamp(cos_alpha, static_cast<real_t>(-1.0), static_cast<real_t>(1.0)));
  const real_t beta_rad = std::acos(std::clamp(cos_beta, static_cast<real_t>(-1.0), static_cast<real_t>(1.0)));
  const real_t gamma_rad = std::acos(std::clamp(cos_gamma, static_cast<real_t>(-1.0), static_cast<real_t>(1.0)));

  lattice_parameters_ = {len_a,
                         len_b,
                         len_c,
                         alpha_rad * static_cast<real_t>(math::rad_to_deg),
                         beta_rad * static_cast<real_t>(math::rad_to_deg),
                         gamma_rad * static_cast<real_t>(math::rad_to_deg)};
}

std::optional<Element> Cell::findElement(const std::string &symbol) const {
  auto iter =
      std::find_if(elements_.begin(), elements_.end(), [&](const Element &elem) { return elem.symbol == symbol; });
  if (iter != elements_.end()) {
    return *iter;
  }
  return std::nullopt;
}

ElementID Cell::getOrRegisterElement(const std::string &symbol) {
  if (auto existing_element = findElement(symbol)) {
    return existing_element->id;
  }
  // Register the new element
  ElementID new_id{static_cast<int>(elements_.size())};
  elements_.push_back({.symbol = symbol, .id = new_id});
  return new_id;
}

Atom &Cell::addAtom(const std::string &symbol, const math::Vector3<real_t> &position) {
  ElementID element_id = getOrRegisterElement(symbol);
  auto element_it = std::find_if(elements_.begin(), elements_.end(),
                                 [&](const Element &elem) { return elem.id.value == element_id.value; });

  AtomID const new_atom_id{static_cast<std::uint32_t>(atoms_.size())};
  atoms_.emplace_back(*element_it, position, new_atom_id);
  return atoms_.back();
}

math::Vector3<real_t> Cell::minimumImage(const math::Vector3<real_t> &distance) const {
  // Convert Cartesian distance to fractional coordinates
  math::Vector3<real_t> frac_dist = inverse_lattice_vectors_ * distance;

  // Apply minimum image convention: shift to [-0.5, 0.5)
  frac_dist.x() -= std::round(frac_dist.x());
  frac_dist.y() -= std::round(frac_dist.y());
  frac_dist.z() -= std::round(frac_dist.z());

  // Convert back to Cartesian coordinates
  return lattice_vectors_ * frac_dist;
}

void Cell::wrapPositions() {
  for (Atom &atom : atoms_) {
    math::Vector3<real_t> frac_pos = inverse_lattice_vectors_ * atom.position();
    frac_pos.x() -= std::floor(frac_pos.x());
    frac_pos.y() -= std::floor(frac_pos.y());
    frac_pos.z() -= std::floor(frac_pos.z());
    atom.setPosition(lattice_vectors_ * frac_pos);
  }
}

} // namespace correlation::core
