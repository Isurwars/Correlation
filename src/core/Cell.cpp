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
#include <limits>
#include <stdexcept>

namespace correlation::core {

Cell::Cell(const math::Vector3<real_t> &vec_a, const math::Vector3<real_t> &vec_b,
           const math::Vector3<real_t> &vec_c) {
  updateLattice(math::Matrix3<real_t>(vec_a, vec_b, vec_c));
}

Cell::Cell(const std::array<real_t, 6> &params) { setLatticeParameters(params); }

void Cell::setLatticeParameters(const std::array<real_t, 6> &params) {
  lattice_parameters_ = params;
  const real_t len_a = params[0];
  const real_t len_b = params[1];
  const real_t len_c = params[2];
  const real_t alpha = params[3] * static_cast<real_t>(math::deg_to_rad);
  const real_t beta = params[4] * static_cast<real_t>(math::deg_to_rad);
  const real_t gamma = params[5] * static_cast<real_t>(math::deg_to_rad);

  if (std::isnan(len_a) || std::isnan(len_b) || std::isnan(len_c) || len_a <= 0 || len_b <= 0 ||
      len_c <= 0) {
    throw std::invalid_argument("Lattice parameters a, b, c must be positive.");
  }

  if (std::isnan(params[3]) || std::isnan(params[4]) || std::isnan(params[5]) || params[3] <= 0 ||
      params[3] >= 180 || params[4] <= 0 || params[4] >= 180 || params[5] <= 0 ||
      params[5] >= 180) {
    throw std::invalid_argument("Lattice angles must be between 0 and 180 degrees.");
  }

  const real_t cos_g = std::cos(gamma);
  const real_t sin_g = std::sin(gamma);

  math::Vector3<real_t> const v_a = {len_a, 0.0, 0.0};
  math::Vector3<real_t> const v_b = {len_b * cos_g, len_b * sin_g, 0.0};
  math::Vector3<real_t> v_c = {len_c * std::cos(beta),
                               len_c * (std::cos(alpha) - std::cos(beta) * cos_g) / sin_g, 0.0};
  const auto volume =
      static_cast<real_t>(len_a * len_b * len_c *
                          std::sqrt(1.0 - std::pow(std::cos(alpha), 2) -
                                    std::pow(std::cos(beta), 2) - std::pow(std::cos(gamma), 2) +
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

  if (std::isnan(len_a) || std::isnan(len_b) || std::isnan(len_c) || len_a < 1e-9 || len_b < 1e-9 ||
      len_c < 1e-9) {
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

  const real_t alpha_rad =
      std::acos(std::clamp(cos_alpha, static_cast<real_t>(-1.0), static_cast<real_t>(1.0)));
  const real_t beta_rad =
      std::acos(std::clamp(cos_beta, static_cast<real_t>(-1.0), static_cast<real_t>(1.0)));
  const real_t gamma_rad =
      std::acos(std::clamp(cos_gamma, static_cast<real_t>(-1.0), static_cast<real_t>(1.0)));

  lattice_parameters_ = {len_a,
                         len_b,
                         len_c,
                         alpha_rad * static_cast<real_t>(math::rad_to_deg),
                         beta_rad * static_cast<real_t>(math::rad_to_deg),
                         gamma_rad * static_cast<real_t>(math::rad_to_deg)};
}

std::optional<Element> Cell::findElement(std::string_view symbol) const {
  auto iter =
      std::ranges::find_if(elements_, [&](const Element &elem) { return elem.symbol == symbol; });
  if (iter != elements_.end()) {
    return *iter;
  }
  return std::nullopt;
}

const Element &Cell::getOrRegisterElement(std::string_view symbol) {
  auto iter =
      std::ranges::find_if(elements_, [&](const Element &elem) { return elem.symbol == symbol; });
  if (iter != elements_.end()) {
    return *iter;
  }
  // Register the new element
  ElementID const new_id{static_cast<int>(elements_.size())};
  elements_.push_back({
      .symbol = std::string(symbol),
      .id = new_id,
  });
  return elements_.back();
}

Atom &Cell::addAtom(std::string_view symbol, const math::Vector3<real_t> &position) {
  const Element &element = getOrRegisterElement(symbol);
  AtomID const new_atom_id{static_cast<std::uint32_t>(atoms_.size())};
  atoms_.emplace_back(element, position, new_atom_id);
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

std::array<real_t, 3> Cell::perpendicularWidths() const noexcept {
  const real_t eps = std::numeric_limits<real_t>::epsilon();
  if (volume_ <= eps) {
    return {0.0, 0.0, 0.0};
  }
  const auto &lat_v0 = lattice_vectors_[0];
  const auto &lat_v1 = lattice_vectors_[1];
  const auto &lat_v2 = lattice_vectors_[2];

  const real_t norm_bc = math::norm(math::cross(lat_v1, lat_v2));
  const real_t norm_ac = math::norm(math::cross(lat_v0, lat_v2));
  const real_t norm_ab = math::norm(math::cross(lat_v0, lat_v1));

  const real_t w_a = (norm_bc > eps) ? (volume_ / norm_bc) : static_cast<real_t>(0.0);
  const real_t w_b = (norm_ac > eps) ? (volume_ / norm_ac) : static_cast<real_t>(0.0);
  const real_t w_c = (norm_ab > eps) ? (volume_ / norm_ab) : static_cast<real_t>(0.0);

  return {w_a, w_b, w_c};
}

Cell Cell::replicate(int rep_a, int rep_b, int rep_c) const {
  if (rep_a < 1 || rep_b < 1 || rep_c < 1) {
    throw std::invalid_argument(
        "Replication factors must each be >= 1, got: " + std::to_string(rep_a) + "x" +
        std::to_string(rep_b) + "x" + std::to_string(rep_c));
  }
  if (rep_a == 1 && rep_b == 1 && rep_c == 1) {
    return *this;
  }

  const size_t total_rep =
      static_cast<size_t>(rep_a) * static_cast<size_t>(rep_b) * static_cast<size_t>(rep_c);
  const size_t total_atoms = total_rep * atoms_.size();
  constexpr size_t MAX_SUPERCELL_ATOMS = 50000;
  if (total_atoms > MAX_SUPERCELL_ATOMS) {
    throw std::invalid_argument("Supercell atom count (" + std::to_string(total_atoms) +
                                ") exceeds safety ceiling of " +
                                std::to_string(MAX_SUPERCELL_ATOMS) + " atoms");
  }

  const auto &vec_a = lattice_vectors_[0];
  const auto &vec_b = lattice_vectors_[1];
  const auto &vec_c = lattice_vectors_[2];

  Cell new_cell(vec_a * static_cast<real_t>(rep_a), vec_b * static_cast<real_t>(rep_b),
                vec_c * static_cast<real_t>(rep_c));
  new_cell.reserveAtoms(total_atoms);

  for (int idx_a = 0; idx_a < rep_a; ++idx_a) {
    for (int idx_b = 0; idx_b < rep_b; ++idx_b) {
      for (int idx_c = 0; idx_c < rep_c; ++idx_c) {
        const math::Vector3<real_t> displacement = vec_a * static_cast<real_t>(idx_a) +
                                                   vec_b * static_cast<real_t>(idx_b) +
                                                   vec_c * static_cast<real_t>(idx_c);

        for (const auto &atom : atoms_) {
          new_cell.addAtom(atom.element().symbol, atom.position() + displacement);
        }
      }
    }
  }

  return new_cell;
}

Cell Cell::autoSupercell(real_t r_cut, int max_replication, real_t max_radius) const {
  if (r_cut <= static_cast<real_t>(0.0)) {
    return *this;
  }
  if (r_cut > max_radius) {
    throw std::invalid_argument("Cutoff radius r_cut (" + std::to_string(r_cut) +
                                " Å) exceeds safe maximum radius of " + std::to_string(max_radius) +
                                " Å");
  }

  const auto widths = perpendicularWidths();
  const real_t eps = std::numeric_limits<real_t>::epsilon();
  if (widths[0] <= eps || widths[1] <= eps || widths[2] <= eps) {
    return *this;
  }

  const real_t target_width = static_cast<real_t>(2.0) * r_cut;
  const auto req_a = static_cast<int>(std::ceil(target_width / widths[0]));
  const auto req_b = static_cast<int>(std::ceil(target_width / widths[1]));
  const auto req_c = static_cast<int>(std::ceil(target_width / widths[2]));

  if (req_a > max_replication || req_b > max_replication || req_c > max_replication) {
    throw std::invalid_argument("Required auto-supercell replication (" + std::to_string(req_a) +
                                "x" + std::to_string(req_b) + "x" + std::to_string(req_c) +
                                ") exceeds safe ceiling of " + std::to_string(max_replication) +
                                " per axis");
  }

  const int rep_a = std::max(1, req_a);
  const int rep_b = std::max(1, req_b);
  const int rep_c = std::max(1, req_c);

  return replicate(rep_a, rep_b, rep_c);
}

} // namespace correlation::core
