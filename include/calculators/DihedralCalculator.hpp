/**
 * @file DihedralCalculator.hpp
 * @brief Dihedral angle geometry utilities.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include <vector>

// Forward declarations
class DistributionFunctions;
struct AnalysisSettings;

namespace calculators {

// [correlation::core::Element A] [correlation::core::Element B]
// [correlation::core::Element C] [correlation::core::Element D] -> List of
// Angles
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
class DihedralCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "Dihedral"; }
  std::string getShortName() const override { return "Dihedral"; }
  std::string getGroup() const override { return "Structural"; }
  std::string getDescription() const override {
    return "Computes all unique 4-body dihedral (torsion) angles.";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(DistributionFunctions &df,
                      const AnalysisSettings &settings) const override;

  /**
   * @brief Computes and populates the dihedral angles tensor.
   *
   * @param cell The simulation cell containing atoms and positions.
   * @param graph The pre-computed neighbor graph containing local connections.
   * @param out_dihedrals A 5D tensor `[correlation::core::Element
   * A][correlation::core::Element B][correlation::core::Element
   * C][correlation::core::Element D][angle_idx]` where elements are the
   * indices, populated with angles in radians [-pi, pi].
   */
  static void compute(const correlation::core::Cell &cell,
                      const correlation::core::NeighborGraph &graph,
                      DihedralTensor &out_dihedrals);
};

} // namespace calculators
