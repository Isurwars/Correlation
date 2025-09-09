#ifndef INCLUDE_STRUCTUREANALIZER_HPP_
#define INCLUDE_STRUCTUREANALIZER_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <vector>

#include "Cell.hpp"

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
  using NeighborTensor = std::vector<std::vector<Neighbor>>;
  using DistanceTensor = std::vector<std::vector<std::vector<double>>>;
  using AngleTensor =
      std::vector<std::vector<std::vector<std::vector<double>>>>;

  /**
   * @brief Computes all distances and angles within the cutoff radius.
   * @param cell The cell containing the atoms and lattice information.
   * @param cutoff The maximum distance to search for neighbors.
   * @param bond_factor A factor to multiply with the sum of covalent radii to
   * determine if two atoms are bonded.
   */
  explicit StructureAnalyzer(const Cell &cell, double cutoff = 20.0,
                             double bond_factor = 1.2);

  // --- Accessors for the computed data tensors ---
  const DistanceTensor &distances() const { return distance_tensor_; }
  const AngleTensor &angles() const { return angle_tensor_; }
  const NeighborTensor &neighbors() const { return neighbor_tensor_; }

private:
  const Cell &cell_;
  double cutoff_sq_;
  double bond_factor_;

  NeighborTensor neighbor_tensor_;
  DistanceTensor distance_tensor_;
  AngleTensor angle_tensor_;
  std::vector<std::vector<double>> bond_cutoffs_sq_;

  /**
   * @brief Pre-calculates the squared bond cutoff distances for every pair of
   * element types.
   */
  void precomputeBondCutoffs();

  /**
   * @brief Computes pairwise distances and identifies bonded neighbors.
   * This method is parallelized using OpenMP.
   * @return A list of all neighbors within the cutoff for each atom, used for
   * angle calculations.
   */
  void computeDistances();

  /**
   * @brief Computes all unique bond angles.
   * This method is parallelized using OpenMP.
   */
  void computeAngles();
};

#endif // INCLUDE_STRUCTUREANALYZER_HPP_
