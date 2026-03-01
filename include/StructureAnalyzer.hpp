// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "Cell.hpp"
#include "NeighborGraph.hpp"
#include "calculators/DihedralCalculator.hpp"

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
  using DistanceTensor = std::vector<std::vector<std::vector<double>>>;
  using AngleTensor =
      std::vector<std::vector<std::vector<std::vector<double>>>>;

  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//
  explicit StructureAnalyzer(
      Cell &cell, double cutoff,
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
  const calculators::DihedralTensor &dihedrals() const {
    return dihedral_tensor_;
  }

  /**
   * @brief Gets the corresponding neighbor graph capturing topological
   * connections.
   * @return Constant reference to the NeighborGraph object.
   */
  const NeighborGraph &neighborGraph() const { return neighbor_graph_; }

private:
  Cell &cell_;
  double cutoff_sq_;
  std::vector<std::vector<double>> bond_cutoffs_sq_;

  bool ignore_periodic_self_interactions_;

  NeighborGraph neighbor_graph_;
  DistanceTensor distance_tensor_;
  AngleTensor angle_tensor_;
  calculators::DihedralTensor dihedral_tensor_;
};
