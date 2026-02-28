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
  const DistanceTensor &distances() const { return distance_tensor_; }
  const AngleTensor &angles() const { return angle_tensor_; }
  const calculators::DihedralTensor &dihedrals() const {
    return dihedral_tensor_;
  }
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
