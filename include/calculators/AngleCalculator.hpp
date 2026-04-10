/**
 * @file AngleCalculator.hpp
 * @brief Calculator for angular distribution functions.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"
#include <vector>

namespace calculators {

using AngleTensor = std::vector<std::vector<std::vector<std::vector<double>>>>;

/**
 * @class AngleCalculator
 * @brief Computes all unique 3-body bond angles in the simulation cell.
 *
 * This calculator iterates over every central atom and finds all unique
 * combinations of its neighbors. It calculates the bond angle formed by the
 * triad (Neighbor1 - Central_Atom - Neighbor2) taking into account
 * periodic boundary conditions.
 *
 * The computation is accelerated using Intel TBB (`tbb::parallel_for`) with
 * thread-local storage (`tbb::enumerable_thread_specific`) to avoid locks when
 * writing to the output tensor.
 */
class AngleCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "Angle"; }
  std::string getShortName() const override { return "BAD"; }
  std::string getGroup() const override { return "Structural"; }
  std::string getDescription() const override {
    return "Computes the Bond Angle Distribution (BAD).";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(DistributionFunctions &df,
                      const AnalysisSettings &settings) const override;

  /**
   * @brief Computes and populates the bond angles tensor.
   *
   * @param cell The simulation cell containing atoms and positions.
   * @param graph The pre-computed neighbor graph containing local connections.
   * @param out_angles A 4D tensor `[outer1][central][outer2][angle_idx]`
   * populated with angles in radians.
   */
  static void compute(const correlation::core::Cell &cell, const correlation::core::NeighborGraph &graph,
                      AngleTensor &out_angles);
};

} // namespace calculators
