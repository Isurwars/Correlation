/**
 * @file Scene3D.hpp
 * @brief Scene graph and geometric instances (atoms, bonds, lattice lines) for 3D visualization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Cell.hpp"
#include "render3d/ElementColors.hpp"
#include "render3d/Math3D.hpp"
#include <memory>
#include <string>
#include <vector>

namespace correlation::render3d {

/**
 * @struct AtomInstance
 * @brief Visual instance of an atom for 3D rendering.
 */
struct AtomInstance {
  math::Vector3<float> position; ///< World coordinates (Angstroms).
  float radius{0.5f};            ///< Visualization radius (Angstroms).
  float covalent_radius{0.75f};  ///< Physical covalent radius for bonding (Angstroms).
  RgbaColor color;               ///< CPK element color.
  uint32_t id{0};                ///< Atom index or identifier.
  int atomic_number{0};          ///< Atomic number (Z).
};

/**
 * @struct BondInstance
 * @brief Visual instance of a covalent or ionic bond connecting two atoms.
 */
struct BondInstance {
  math::Vector3<float> start; ///< Starting atom position.
  math::Vector3<float> end;   ///< Ending atom position.
  float radius{0.12f};        ///< Cylinder visual radius.
  RgbaColor color;            ///< Bond color (interpolated or neutral).
};

/**
 * @struct LineSegment
 * @brief 3D line segment for wireframe rendering (e.g. simulation box boundaries).
 */
struct LineSegment {
  math::Vector3<float> start; ///< Starting vertex.
  math::Vector3<float> end;   ///< Ending vertex.
  RgbaColor color;            ///< Line color.
};

/**
 * @class Scene3D
 * @brief Aggregates atomic coordinates, bonds, and unit cell wireframes for rendering.
 */
class Scene3D {
public:
  Scene3D() = default;

  /**
   * @brief Constructs a Scene3D from a simulation Cell.
   * @param cell Atomic cell containing lattice and atoms.
   * @param bond_cutoff_scale Scaling multiplier for covalent radius cutoff.
   * @param generate_bonds Flag enabling covalent bond detection.
   * @return Shared pointer to initialized Scene3D.
   */
  [[nodiscard]] static std::shared_ptr<Scene3D>
  fromCell(const core::Cell &cell, double bond_cutoff_scale = 1.15, bool generate_bonds = true);

  /**
   * @brief Clears all geometry in the scene.
   */
  void clear() noexcept;

  [[nodiscard]] const std::vector<AtomInstance> &atoms() const noexcept { return atoms_; }
  [[nodiscard]] const std::vector<BondInstance> &bonds() const noexcept { return bonds_; }
  [[nodiscard]] const std::vector<LineSegment> &unitCellLines() const noexcept {
    return unit_cell_lines_;
  }

  [[nodiscard]] const math::Vector3<float> &boundsMin() const noexcept { return bounds_min_; }
  [[nodiscard]] const math::Vector3<float> &boundsMax() const noexcept { return bounds_max_; }

  void addAtom(const AtomInstance &atom);
  void addBond(const BondInstance &bond);
  void addUnitCellLine(const LineSegment &line);

  /**
   * @brief Updates atomic positions and unit cell in-place from a new Cell snapshot.
   * @param cell New simulation snapshot.
   * @param recompute_bonds If true, recalculates bond pairs.
   * @param bond_cutoff_scale Scaling multiplier for cutoff.
   * @return True if update succeeded without rebuilding the scene graph.
   */
  bool updateFromCell(const core::Cell &cell, bool recompute_bonds = false,
                      double bond_cutoff_scale = 1.15);

  /**
   * @brief Fast in-place update of atom coordinates from Cell atoms vector.
   * @param cell_atoms Source atom list.
   */
  void updateAtomPositions(const std::vector<core::Atom> &cell_atoms);

  /**
   * @brief Recomputes axis-aligned bounding box from atoms and unit cell vertices.
   */
  void updateBounds() noexcept;

private:
  std::vector<AtomInstance> atoms_;
  std::vector<BondInstance> bonds_;
  std::vector<LineSegment> unit_cell_lines_;
  math::Vector3<float> bounds_min_{0.0f, 0.0f, 0.0f};
  math::Vector3<float> bounds_max_{0.0f, 0.0f, 0.0f};
};

} // namespace correlation::render3d
