/**
 * @file Scene3D.cpp
 * @brief Implementation of Scene3D molecular geometry and wireframe builder.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/Scene3D.hpp"
#include <algorithm>
#include <cmath>

namespace correlation::render3d {

namespace {

constexpr RgbaColor K_UNIT_CELL_COLOR{.r = 190, .g = 195, .b = 205, .a = 220};
constexpr RgbaColor K_DEFAULT_BOND_COLOR{.r = 200, .g = 200, .b = 200, .a = 255};
constexpr float K_MIN_BOND_DISTANCE = 0.4F;
constexpr std::size_t K_MAX_ATOMS_FOR_FULL_BONDS = 5000;

void generateBoxEdges(const std::array<math::Vector3<float>, 8> &vertices,
                      std::vector<LineSegment> &lines) {
  constexpr std::array<std::pair<std::size_t, std::size_t>, 12> K_EDGE_INDICES = {
      {
          {0, 1},
          {1, 3},
          {3, 2},
          {2, 0}, // Bottom
          {4, 5},
          {5, 7},
          {7, 6},
          {6, 4}, // Top
          {0, 4},
          {1, 5},
          {2, 6},
          {3, 7}, // Vertical pillars
      },
  };

  lines.reserve(lines.size() + K_EDGE_INDICES.size());
  for (const auto &[idx_start, idx_end] : K_EDGE_INDICES) {
    lines.push_back(LineSegment{
        .start = vertices.at(idx_start),
        .end = vertices.at(idx_end),
        .color = K_UNIT_CELL_COLOR,
    });
  }
}

void extractUnitCellLines(const math::Matrix3<real_t> &lattice,
                          std::vector<LineSegment> &out_lines) {
  const math::Vector3<float> a_vec(lattice[0]);
  const math::Vector3<float> b_vec(lattice[1]);
  const math::Vector3<float> c_vec(lattice[2]);

  // Only generate unit cell wireframe if vectors are non-zero
  if (math::norm(a_vec) < 1e-4F && math::norm(b_vec) < 1e-4F && math::norm(c_vec) < 1e-4F) {
    return;
  }

  const std::array<math::Vector3<float>, 8> vertices = {
      {
          math::Vector3<float>(0.0F, 0.0F, 0.0F),
          a_vec,
          b_vec,
          a_vec + b_vec,
          c_vec,
          a_vec + c_vec,
          b_vec + c_vec,
          a_vec + b_vec + c_vec,
      },
  };

  generateBoxEdges(vertices, out_lines);
}

void extractAtoms(const std::vector<core::Atom> &cell_atoms, const math::Matrix3<real_t> &lattice,
                  const math::Matrix3<real_t> &inv_lattice, bool has_lattice, bool in_cell,
                  float atom_scale, std::vector<AtomInstance> &out_atoms) {
  out_atoms.reserve(cell_atoms.size());
  for (const auto &atom : cell_atoms) {
    const std::string_view symbol = atom.element().symbol;
    auto pos = atom.position();

    if (in_cell && has_lattice) {
      math::Vector3<real_t> frac_pos = inv_lattice * pos;
      frac_pos.x() -= std::floor(frac_pos.x());
      frac_pos.y() -= std::floor(frac_pos.y());
      frac_pos.z() -= std::floor(frac_pos.z());
      pos = lattice * frac_pos;
    }

    const float cov_radius = getElementDisplayRadius(symbol, 1.0F);
    out_atoms.push_back(AtomInstance{
        .position = math::Vector3<float>(pos),
        .radius = cov_radius * atom_scale,
        .covalent_radius = cov_radius,
        .color = getElementColor(symbol),
        .id = static_cast<uint32_t>(atom.id()),
        .atomic_number = 0,
    });
  }
}

void extractBonds(const std::vector<AtomInstance> &atoms, double cutoff_scale,
                  std::vector<BondInstance> &out_bonds) {
  if (atoms.size() > K_MAX_ATOMS_FOR_FULL_BONDS) {
    return;
  }

  const auto num_atoms = atoms.size();
  for (std::size_t i = 0; i < num_atoms; ++i) {
    for (std::size_t j = i + 1; j < num_atoms; ++j) {
      const math::Vector3<float> diff = atoms[i].position - atoms[j].position;
      const float dist = math::norm(diff);
      const auto max_cutoff =
          static_cast<float>((atoms[i].covalent_radius + atoms[j].covalent_radius) * cutoff_scale);

      if (dist >= K_MIN_BOND_DISTANCE && dist <= max_cutoff) {
        out_bonds.push_back(BondInstance{
            .start = atoms[i].position,
            .end = atoms[j].position,
            .radius = 0.08F,
            .color = K_DEFAULT_BOND_COLOR,
        });
      }
    }
  }
}

bool clipTest(float param_p, float param_q, float &t_min, float &t_max) noexcept {
  if (param_p < 0.0F) {
    const float ratio = param_q / param_p;
    if (ratio > t_max) {
      return false;
    }
    t_min = std::max(ratio, t_min);
  } else if (param_p > 0.0F) {
    const float ratio = param_q / param_p;
    if (ratio < t_min) {
      return false;
    }
    t_max = std::min(ratio, t_max);
  } else if (param_q < 0.0F) {
    return false;
  }
  return true;
}

bool clipSegmentToUnitBox(const math::Vector3<float> &seg_start,
                          const math::Vector3<float> &seg_end, math::Vector3<float> &out_start,
                          math::Vector3<float> &out_end) noexcept {
  float t_min = 0.0F;
  float t_max = 1.0F;
  const float delta_x = seg_end.x() - seg_start.x();
  const float delta_y = seg_end.y() - seg_start.y();
  const float delta_z = seg_end.z() - seg_start.z();

  const bool inside = clipTest(-delta_x, seg_start.x(), t_min, t_max) &&
                      clipTest(delta_x, 1.0F - seg_start.x(), t_min, t_max) &&
                      clipTest(-delta_y, seg_start.y(), t_min, t_max) &&
                      clipTest(delta_y, 1.0F - seg_start.y(), t_min, t_max) &&
                      clipTest(-delta_z, seg_start.z(), t_min, t_max) &&
                      clipTest(delta_z, 1.0F - seg_start.z(), t_min, t_max) && (t_min <= t_max);

  if (!inside) {
    return false;
  }

  out_start = math::Vector3<float>(seg_start.x() + t_min * delta_x, seg_start.y() + t_min * delta_y,
                                   seg_start.z() + t_min * delta_z);
  out_end = math::Vector3<float>(seg_start.x() + t_max * delta_x, seg_start.y() + t_max * delta_y,
                                 seg_start.z() + t_max * delta_z);
  return true;
}

math::Vector3<float> toCartesian(const math::Vector3<float> &frac,
                                 const math::Matrix3<real_t> &lattice) noexcept {
  return math::Vector3<float>(lattice[0] * frac.x() + lattice[1] * frac.y() +
                              lattice[2] * frac.z());
}

void appendInCellBonds(const AtomInstance &atom_i, const AtomInstance &atom_j,
                       const math::Vector3<float> &frac_i, const math::Vector3<float> &frac_j,
                       const math::Matrix3<real_t> &lattice, double cutoff_scale,
                       std::vector<BondInstance> &out_bonds) {
  const auto max_cutoff =
      static_cast<float>((atom_i.covalent_radius + atom_j.covalent_radius) * cutoff_scale);

  const math::Vector3<float> delta = frac_j - frac_i;
  const float shift_x = -std::round(delta.x());
  const float shift_y = -std::round(delta.y());
  const float shift_z = -std::round(delta.z());
  const math::Vector3<float> shift_vec(shift_x, shift_y, shift_z);

  const math::Vector3<float> mic_delta = delta + shift_vec;
  const math::Vector3<float> cart_diff = toCartesian(mic_delta, lattice);

  const float dist = math::norm(cart_diff);
  if (dist < K_MIN_BOND_DISTANCE || dist > max_cutoff) {
    return;
  }

  // Intra-cell bond (no periodic boundary crossing)
  if (std::abs(shift_x) < 0.5F && std::abs(shift_y) < 0.5F && std::abs(shift_z) < 0.5F) {
    out_bonds.push_back(BondInstance{
        .start = atom_i.position,
        .end = atom_j.position,
        .radius = 0.08F,
        .color = K_DEFAULT_BOND_COLOR,
    });
    return;
  }

  // Crossing bond: create boundary-clipped stubs
  math::Vector3<float> clipped_start;
  math::Vector3<float> clipped_end;
  if (clipSegmentToUnitBox(frac_i, frac_j + shift_vec, clipped_start, clipped_end)) {
    out_bonds.push_back(BondInstance{
        .start = atom_i.position,
        .end = toCartesian(clipped_end, lattice),
        .radius = 0.08F,
        .color = K_DEFAULT_BOND_COLOR,
    });
  }

  if (clipSegmentToUnitBox(frac_j, frac_i - shift_vec, clipped_start, clipped_end)) {
    out_bonds.push_back(BondInstance{
        .start = atom_j.position,
        .end = toCartesian(clipped_end, lattice),
        .radius = 0.08F,
        .color = K_DEFAULT_BOND_COLOR,
    });
  }
}

void extractInCellBonds(const std::vector<AtomInstance> &atoms,
                        const math::Matrix3<real_t> &lattice,
                        const math::Matrix3<real_t> &inv_lattice, double cutoff_scale,
                        std::vector<BondInstance> &out_bonds) {
  if (atoms.size() > K_MAX_ATOMS_FOR_FULL_BONDS) {
    return;
  }

  const auto num_atoms = atoms.size();
  std::vector<math::Vector3<float>> frac_coords;
  frac_coords.reserve(num_atoms);

  for (const auto &atom : atoms) {
    const math::Vector3<real_t> cart(atom.position.x(), atom.position.y(), atom.position.z());
    math::Vector3<real_t> frac_raw = inv_lattice * cart;
    frac_raw.x() -= std::floor(frac_raw.x());
    frac_raw.y() -= std::floor(frac_raw.y());
    frac_raw.z() -= std::floor(frac_raw.z());
    const auto frac_float = math::Vector3<float>(frac_raw);
    const float clamped_x = std::clamp(frac_float.x(), 0.0F, 0.999999F);
    const float clamped_y = std::clamp(frac_float.y(), 0.0F, 0.999999F);
    const float clamped_z = std::clamp(frac_float.z(), 0.0F, 0.999999F);
    frac_coords.emplace_back(clamped_x, clamped_y, clamped_z);
  }

  for (std::size_t i = 0; i < num_atoms; ++i) {
    for (std::size_t j = i + 1; j < num_atoms; ++j) {
      appendInCellBonds(atoms[i], atoms[j], frac_coords[i], frac_coords[j], lattice, cutoff_scale,
                        out_bonds);
    }
  }
}

} // namespace

void Scene3D::clear() noexcept {
  atoms_.clear();
  bonds_.clear();
  unit_cell_lines_.clear();
  bounds_min_ = math::Vector3<float>(0.0F, 0.0F, 0.0F);
  bounds_max_ = math::Vector3<float>(0.0F, 0.0F, 0.0F);
  has_lattice_ = false;
  in_cell_ = false;
}

void Scene3D::addAtom(const AtomInstance &atom) { atoms_.push_back(atom); }

void Scene3D::addBond(const BondInstance &bond) { bonds_.push_back(bond); }

void Scene3D::addUnitCellLine(const LineSegment &line) { unit_cell_lines_.push_back(line); }

void Scene3D::setAtomScale(float atom_scale) noexcept {
  atom_scale_ = std::clamp(atom_scale, 0.05F, 2.0F);
  for (auto &atom : atoms_) {
    atom.radius = atom.covalent_radius * atom_scale_;
  }
}

void Scene3D::updateBounds() noexcept {
  if (atoms_.empty() && unit_cell_lines_.empty()) {
    bounds_min_ = math::Vector3<float>(-5.0F, -5.0F, -5.0F);
    bounds_max_ = math::Vector3<float>(5.0F, 5.0F, 5.0F);
    return;
  }

  float min_x = 1e9F;
  float min_y = 1e9F;
  float min_z = 1e9F;
  float max_x = -1e9F;
  float max_y = -1e9F;
  float max_z = -1e9F;

  for (const auto &atom : atoms_) {
    min_x = std::min(min_x, atom.position.x() - atom.radius);
    min_y = std::min(min_y, atom.position.y() - atom.radius);
    min_z = std::min(min_z, atom.position.z() - atom.radius);
    max_x = std::max(max_x, atom.position.x() + atom.radius);
    max_y = std::max(max_y, atom.position.y() + atom.radius);
    max_z = std::max(max_z, atom.position.z() + atom.radius);
  }

  for (const auto &line : unit_cell_lines_) {
    min_x = std::min({min_x, line.start.x(), line.end.x()});
    min_y = std::min({min_y, line.start.y(), line.end.y()});
    min_z = std::min({min_z, line.start.z(), line.end.z()});
    max_x = std::max({max_x, line.start.x(), line.end.x()});
    max_y = std::max({max_y, line.start.y(), line.end.y()});
    max_z = std::max({max_z, line.start.z(), line.end.z()});
  }

  bounds_min_ = math::Vector3<float>(min_x, min_y, min_z);
  bounds_max_ = math::Vector3<float>(max_x, max_y, max_z);
}

void Scene3D::updateAtomPositions(const std::vector<core::Atom> &cell_atoms, bool in_cell) {
  in_cell_ = in_cell;
  const auto update_count = std::min(atoms_.size(), cell_atoms.size());
  for (std::size_t i = 0; i < update_count; ++i) {
    auto pos = cell_atoms[i].position();
    if (in_cell_ && has_lattice_) {
      math::Vector3<real_t> frac_pos = inv_lattice_vectors_ * pos;
      frac_pos.x() -= std::floor(frac_pos.x());
      frac_pos.y() -= std::floor(frac_pos.y());
      frac_pos.z() -= std::floor(frac_pos.z());
      pos = lattice_vectors_ * frac_pos;
    }
    atoms_[i].position = math::Vector3<float>(pos);
  }
}

void Scene3D::recomputeBonds(double bond_cutoff_scale, bool in_cell) {
  in_cell_ = in_cell;
  bonds_.clear();
  if (in_cell_ && has_lattice_) {
    extractInCellBonds(atoms_, lattice_vectors_, inv_lattice_vectors_, bond_cutoff_scale, bonds_);
  } else {
    extractBonds(atoms_, bond_cutoff_scale, bonds_);
  }
}

bool Scene3D::updateFromCell(const core::Cell &cell, bool recompute_bonds, double bond_cutoff_scale,
                             bool in_cell) {
  lattice_vectors_ = cell.latticeVectors();
  inv_lattice_vectors_ = cell.inverseLatticeVectors();
  has_lattice_ = (cell.volume() > 1e-4);
  in_cell_ = in_cell;

  const auto &cell_atoms = cell.atoms();
  if (atoms_.size() != cell_atoms.size()) {
    unit_cell_lines_.clear();
    atoms_.clear();
    bonds_.clear();
    extractUnitCellLines(lattice_vectors_, unit_cell_lines_);
    extractAtoms(cell_atoms, lattice_vectors_, inv_lattice_vectors_, has_lattice_, in_cell_,
                 atom_scale_, atoms_);
    if (recompute_bonds) {
      recomputeBonds(bond_cutoff_scale, in_cell_);
    }
    updateBounds();
    return false;
  }

  updateAtomPositions(cell_atoms, in_cell_);

  unit_cell_lines_.clear();
  extractUnitCellLines(lattice_vectors_, unit_cell_lines_);

  if (recompute_bonds) {
    recomputeBonds(bond_cutoff_scale, in_cell_);
  }

  updateBounds();
  return true;
}

std::shared_ptr<Scene3D> Scene3D::fromCell(const core::Cell &cell, double bond_cutoff_scale,
                                           bool generate_bonds, bool in_cell) {
  auto scene = std::make_shared<Scene3D>();
  scene->lattice_vectors_ = cell.latticeVectors();
  scene->inv_lattice_vectors_ = cell.inverseLatticeVectors();
  scene->has_lattice_ = (cell.volume() > 1e-4);
  scene->in_cell_ = in_cell;

  extractUnitCellLines(scene->lattice_vectors_, scene->unit_cell_lines_);
  extractAtoms(cell.atoms(), scene->lattice_vectors_, scene->inv_lattice_vectors_,
               scene->has_lattice_, scene->in_cell_, scene->atom_scale_, scene->atoms_);

  if (generate_bonds) {
    scene->recomputeBonds(bond_cutoff_scale, scene->in_cell_);
  }

  scene->updateBounds();
  return scene;
}

} // namespace correlation::render3d
