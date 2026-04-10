/**
 * @file StructureAnalyzer.hpp
 * @brief Structural analysis utilities for neighbor detection and bond
 * topology.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

namespace correlation::analysis {

/**
 * @class StructureAnalyzer
 * @brief Computes and stores pairwise distances and bond angles in tensors
 * indexed by element type.
 *
 * This class efficiently finds all atom pairs within a cutoff radius and all
 * unique bond angles, storing the results in multi-dimensional vectors
 * (tensors) suitable for calculating partial distribution functions.
 * The distance and angle calculation loops are parallelized with OpenMP.
 */
class StructureAnalyzer {
public:
  /** @brief Tensor for storing pair distances [element1][element2][pair_index]. */
  using DistanceTensor = std::vector<std::vector<std::vector<double>>>;
  /** @brief Tensor for storing bond angles [center][outer1][outer2][angle_index]. */
  using AngleTensor =
      std::vector<std::vector<std::vector<std::vector<double>>>>;
  /** @brief Tensor for storing dihedrals [e1][e2][e3][e4][dihedral_index]. */
  using DihedralTensor =
      std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>;

  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//
  /**
   * @brief Constructs a StructureAnalyzer for a single frame (Cell).
   *
   * @param cell The periodic cell containing atomic positions.
   * @param cutoff The neighbor search cutoff radius (Angstrom).
   * @param bond_cutoffs_sq Squared bond cutoffs per element pair.
   * @param ignore_periodic_self_interactions Flag to ignore a-a image pairs.
   */
  explicit StructureAnalyzer(
      correlation::core::Cell &cell, double cutoff,
      const std::vector<std::vector<double>> &bond_cutoffs_sq,
      bool ignore_periodic_self_interactions = true);

  //-------------------------------------------------------------------------//
  //------------------------------- Accessors -------------------------------//
  //-------------------------------------------------------------------------//
  /**
   * @brief Gets a multi-dimensional tensor containing pair distances.
   * @return The distance tensor `[element1][element2][pair_index]`.
   */
  const DistanceTensor &distances() const { return distance_tensor_; }

  /**
   * @brief Gets a multi-dimensional tensor containing bond angles.
   * @return The angle tensor
   * `[center_element][outer_element1][outer_element2][angle_index]`.
   */
  const AngleTensor &angles() const { return angle_tensor_; }

  /**
   * @brief Gets a multi-dimensional tensor containing dihedral angles.
   * @return The dihedral tensor
   * `[element1][element2][element3][element4][dihedral_index]`.
   */
  const DihedralTensor &dihedrals() const {
    return dihedral_tensor_;
  }

  /**
   * @brief Gets the corresponding neighbor graph capturing topological
   * connections.
   * @return Constant reference to the correlation::core::NeighborGraph object.
   */
  const correlation::core::NeighborGraph &neighborGraph() const {
    return neighbor_graph_;
  }

private:
  correlation::core::Cell &cell_;              ///< Reference to the current periodic cell.
  double cutoff_sq_;                         ///< Squared cutoff for efficiency.
  std::vector<std::vector<double>> bond_cutoffs_sq_; ///< Internal squared bond cutoffs.

  bool ignore_periodic_self_interactions_;    ///< Interaction guard.

  correlation::core::NeighborGraph neighbor_graph_; ///< Graph of topological bonds.
  DistanceTensor distance_tensor_;           ///< Pairwise distance storage.
  AngleTensor angle_tensor_;                 ///< Bond angle storage.
  DihedralTensor dihedral_tensor_;           ///< Dihedral angle storage.
};

} // namespace correlation::analysis
