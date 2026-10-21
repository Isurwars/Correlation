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
#include <memory>
#include <vector>

namespace correlation::render3d {

/**
 * @struct AtomInstance
 * @brief Visual instance of an atom for 3D rendering.
 */
struct AtomInstance {
  math::Vector3<float> position; ///< World coordinates (Angstroms).
  float radius{0.5F};            ///< Visualization radius (Angstroms).
  float covalent_radius{0.75F};  ///< Physical covalent radius for bonding (Angstroms).
  RgbaColor color;               ///< CPK element color.
  uint32_t id{0};                ///< Atom index or identifier.
  int atomic_number{0};          ///< Atomic number (Z).

  // Phase Search Halo and styling attributes
  RgbaColor halo_color{.r = 0, .g = 0, .b = 0, .a = 0}; ///< a == 0 denotes no halo
  float halo_radius_ratio{1.35F};                       ///< Expansion over visual radius
  float dim_factor{1.0F};                               ///< 1.0 = normal, 0.3 = ghosted
};

/**
 * @struct BondInstance
 * @brief Visual instance of a covalent or ionic bond connecting two atoms.
 */
struct BondInstance {
  math::Vector3<float> start; ///< Starting atom position.
  math::Vector3<float> end;   ///< Ending atom position.
  float radius{0.12F};        ///< Cylinder visual radius.
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
   * @param[in] cell Atomic cell containing lattice and atoms.
   * @param[in] bond_cutoff_scale Scaling multiplier for covalent radius cutoff.
   * @param[in] generate_bonds Flag enabling covalent bond detection.
   * @param[in] in_cell Flag enabling in-cell coordinate folding and periodic boundary clipped
   * bonds.
   * @return Shared pointer to initialized Scene3D.
   */
  [[nodiscard]] static std::shared_ptr<Scene3D> fromCell(const core::Cell &cell,
                                                         double bond_cutoff_scale = 1.15,
                                                         bool generate_bonds = true,
                                                         bool in_cell = false);

  /**
   * @brief Clears all geometry in the scene.
   */
  void clear() noexcept;

  /**
   * @brief Gets the collection of atom visualization instances in the scene.
   * @return Immutable reference to atom instances.
   */
  [[nodiscard]] const std::vector<AtomInstance> &atoms() const noexcept { return atoms_; }

  /**
   * @brief Gets the collection of bond cylinders in the scene.
   * @return Immutable reference to bond instances.
   */
  [[nodiscard]] const std::vector<BondInstance> &bonds() const noexcept { return bonds_; }

  /**
   * @brief Gets the unit cell boundary line segments for wireframe rendering.
   * @return Immutable reference to unit cell lines.
   */
  [[nodiscard]] const std::vector<LineSegment> &unitCellLines() const noexcept {
    return unit_cell_lines_;
  }

  /**
   * @brief Gets the minimum 3D coordinate bounding box limits.
   * @return Vector of minimum x, y, z coordinates.
   */
  [[nodiscard]] const math::Vector3<float> &boundsMin() const noexcept { return bounds_min_; }

  /**
   * @brief Gets the maximum 3D coordinate bounding box limits.
   * @return Vector of maximum x, y, z coordinates.
   */
  [[nodiscard]] const math::Vector3<float> &boundsMax() const noexcept { return bounds_max_; }

  /**
   * @brief Checks whether in-cell PBC coordinate wrapping is currently active.
   * @return True if in-cell mode is enabled.
   */
  [[nodiscard]] bool inCell() const noexcept { return in_cell_; }

  /**
   * @brief Sets whether in-cell PBC coordinate wrapping is active.
   * @param[in] in_cell True to enable in-cell folding and boundary clipping.
   */
  void setInCell(bool in_cell) noexcept { in_cell_ = in_cell; }

  /**
   * @brief Checks if the scene has valid unit cell lattice vectors.
   * @return True if lattice vectors define a positive simulation volume.
   */
  [[nodiscard]] bool hasLattice() const noexcept { return has_lattice_; }

  /**
   * @brief Gets the 3x3 direct lattice vectors matrix.
   * @return Constant reference to lattice matrix.
   */
  [[nodiscard]] const math::Matrix3<real_t> &latticeVectors() const noexcept {
    return lattice_vectors_;
  }

  /**
   * @brief Gets the 3x3 inverse lattice vectors matrix.
   * @return Constant reference to inverse lattice matrix.
   */
  [[nodiscard]] const math::Matrix3<real_t> &inverseLatticeVectors() const noexcept {
    return inv_lattice_vectors_;
  }

  /**
   * @brief Appends an atom visualization instance to the scene.
   * @param[in] atom Visual atom definition.
   */
  void addAtom(const AtomInstance &atom);

  /**
   * @brief Appends a bond cylinder to the scene.
   * @param[in] bond Visual bond definition.
   */
  void addBond(const BondInstance &bond);

  /**
   * @brief Appends a wireframe line segment to the scene.
   * @param[in] line Line segment definition.
   */
  void addUnitCellLine(const LineSegment &line);

  /**
   * @brief Updates atomic positions and unit cell in-place from a new Cell snapshot.
   * @param[in] cell New simulation snapshot.
   * @param[in] recompute_bonds If true, recalculates bond pairs.
   * @param[in] bond_cutoff_scale Scaling multiplier for cutoff.
   * @param[in] in_cell Flag enabling in-cell coordinate folding and periodic boundary clipped
   * bonds.
   * @return True if update succeeded without rebuilding the scene graph.
   */
  bool updateFromCell(const core::Cell &cell, bool recompute_bonds = false,
                      double bond_cutoff_scale = 1.15, bool in_cell = false);

  /**
   * @brief Recomputes covalent bonds for current atomic coordinates.
   * @param[in] bond_cutoff_scale Scaling multiplier for covalent radius cutoff.
   * @param[in] in_cell Flag enabling periodic minimum image convention and boundary clipping.
   */
  void recomputeBonds(double bond_cutoff_scale = 1.15, bool in_cell = false);

  /**
   * @brief Fast in-place update of atom coordinates from Cell atoms vector.
   * @param[in] cell_atoms Source atom list.
   * @param[in] in_cell Whether to wrap coordinates into [0, 1) fractional box.
   */
  void updateAtomPositions(const std::vector<core::Atom> &cell_atoms, bool in_cell = false);

  /**
   * @brief Recomputes axis-aligned bounding box from atoms and unit cell vertices.
   */
  void updateBounds() noexcept;

private:
  std::vector<AtomInstance> atoms_;
  std::vector<BondInstance> bonds_;
  std::vector<LineSegment> unit_cell_lines_;
  math::Vector3<float> bounds_min_{0.0F, 0.0F, 0.0F};
  math::Vector3<float> bounds_max_{0.0F, 0.0F, 0.0F};

  math::Matrix3<real_t> lattice_vectors_;
  math::Matrix3<real_t> inv_lattice_vectors_;
  bool has_lattice_{false};
  bool in_cell_{false};
};

} // namespace correlation::render3d
