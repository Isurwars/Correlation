// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "../Cell.hpp"
#include "../NeighborGraph.hpp"
#include <vector>

namespace calculators {

// [Element A] [Element B] [Element C] [Element D] -> List of Angles
using DihedralTensor =
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>;

/**
 * @class DihedralCalculator
 * @brief Computes all unique 4-body dihedral (torsion) angles in the simulation
 * cell.
 *
 * This calculator identifies all bonded quadruplets (A-B-C-D) and calculates
 * the dihedral angle defined by the two intersecting planes: plane(A,B,C) and
 * plane(B,C,D). The computation is accelerated using Intel TBB
 * (`tbb::parallel_for`) with thread-local storage.
 */
class DihedralCalculator {
public:
  /**
   * @brief Computes and populates the dihedral angles tensor.
   *
   * @param cell The simulation cell containing atoms and positions.
   * @param graph The pre-computed neighbor graph containing local connections.
   * @param out_dihedrals A 5D tensor `[Element A][Element B][Element C][Element
   * D][angle_idx]` where elements are the indices, populated with angles in
   * radians [-pi, pi].
   */
  static void compute(const Cell &cell, const NeighborGraph &graph,
                      DihedralTensor &out_dihedrals);
};

} // namespace calculators
