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
  const math::Vector3<float> a_vec(static_cast<float>(lattice[0].x()),
                                   static_cast<float>(lattice[0].y()),
                                   static_cast<float>(lattice[0].z()));
  const math::Vector3<float> b_vec(static_cast<float>(lattice[1].x()),
                                   static_cast<float>(lattice[1].y()),
                                   static_cast<float>(lattice[1].z()));
  const math::Vector3<float> c_vec(static_cast<float>(lattice[2].x()),
                                   static_cast<float>(lattice[2].y()),
                                   static_cast<float>(lattice[2].z()));

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

void extractAtoms(const std::vector<core::Atom> &cell_atoms, std::vector<AtomInstance> &out_atoms) {
  out_atoms.reserve(cell_atoms.size());
  for (const auto &atom : cell_atoms) {
    const std::string_view symbol = atom.element().symbol;
    const auto pos = atom.position();

    out_atoms.push_back(AtomInstance{
        .position = math::Vector3<float>(static_cast<float>(pos.x()), static_cast<float>(pos.y()),
                                         static_cast<float>(pos.z())),
        .radius = getElementDisplayRadius(symbol, 0.35F),
        .covalent_radius = getElementDisplayRadius(symbol, 1.0F),
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

} // namespace

void Scene3D::clear() noexcept {
  atoms_.clear();
  bonds_.clear();
  unit_cell_lines_.clear();
  bounds_min_ = math::Vector3<float>(0.0F, 0.0F, 0.0F);
  bounds_max_ = math::Vector3<float>(0.0F, 0.0F, 0.0F);
}

void Scene3D::addAtom(const AtomInstance &atom) { atoms_.push_back(atom); }

void Scene3D::addBond(const BondInstance &bond) { bonds_.push_back(bond); }

void Scene3D::addUnitCellLine(const LineSegment &line) { unit_cell_lines_.push_back(line); }

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

void Scene3D::updateAtomPositions(const std::vector<core::Atom> &cell_atoms) {
  const auto update_count = std::min(atoms_.size(), cell_atoms.size());
  for (std::size_t i = 0; i < update_count; ++i) {
    const auto pos = cell_atoms[i].position();
    atoms_[i].position = math::Vector3<float>(
        static_cast<float>(pos.x()), static_cast<float>(pos.y()), static_cast<float>(pos.z()));
  }
}

void Scene3D::recomputeBonds(double bond_cutoff_scale) {
  bonds_.clear();
  extractBonds(atoms_, bond_cutoff_scale, bonds_);
}

bool Scene3D::updateFromCell(const core::Cell &cell, bool recompute_bonds,
                             double bond_cutoff_scale) {
  const auto &cell_atoms = cell.atoms();
  if (atoms_.size() != cell_atoms.size()) {
    unit_cell_lines_.clear();
    atoms_.clear();
    bonds_.clear();
    extractUnitCellLines(cell.latticeVectors(), unit_cell_lines_);
    extractAtoms(cell_atoms, atoms_);
    if (recompute_bonds) {
      recomputeBonds(bond_cutoff_scale);
    }
    updateBounds();
    return false;
  }

  updateAtomPositions(cell_atoms);

  unit_cell_lines_.clear();
  extractUnitCellLines(cell.latticeVectors(), unit_cell_lines_);

  if (recompute_bonds) {
    recomputeBonds(bond_cutoff_scale);
  }

  updateBounds();
  return true;
}

std::shared_ptr<Scene3D> Scene3D::fromCell(const core::Cell &cell, double bond_cutoff_scale,
                                           bool generate_bonds) {
  auto scene = std::make_shared<Scene3D>();

  extractUnitCellLines(cell.latticeVectors(), scene->unit_cell_lines_);
  extractAtoms(cell.atoms(), scene->atoms_);

  if (generate_bonds) {
    extractBonds(scene->atoms_, bond_cutoff_scale, scene->bonds_);
  }

  scene->updateBounds();
  return scene;
}

} // namespace correlation::render3d
