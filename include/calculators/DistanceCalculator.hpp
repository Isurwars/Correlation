/**
 * @file DistanceCalculator.hpp
 * @brief Pairwise distance calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include <vector>

namespace correlation::calculators {

/** @brief Type definition for a distance tensor [e1][e2][pair_index]. */
using DistanceTensor = std::vector<std::vector<std::vector<double>>>;

/**
 * @class DistanceCalculator
 * @brief High-performance calculator for pairwise atomic distances.
 */
class DistanceCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "Distance"; }
  std::string getShortName() const override { return "Distance"; }
  std::string getGroup() const override { return "Structural"; }
  std::string getDescription() const override {
    return "Computes all unique 2-body distances.";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(
      correlation::analysis::DistributionFunctions &df,
      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief High-performance computation of the distance tensor and neighbor graph.
   * 
   * @param cell The periodic cell container.
   * @param cutoff_sq Squared cutoff radius for neighbor search.
   * @param bond_cutoffs_sq Squared bond cutoffs per element pair.
   * @param ignore_periodic_self_interactions Flag to prevent self-images.
   * @param out_distances Tensor to be populated with pairwise distances.
   * @param out_graph Graph object to be populated with adjacency data.
   */
  static void compute(const correlation::core::Cell &cell, double cutoff_sq,
                      const std::vector<std::vector<double>> &bond_cutoffs_sq,
                      bool ignore_periodic_self_interactions,
                      DistanceTensor &out_distances,
                      correlation::core::NeighborGraph &out_graph);
};

} // namespace correlation::calculators
