// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "BaseCalculator.hpp"
#include "../Cell.hpp"
#include "../NeighborGraph.hpp"
#include <vector>

namespace calculators {

using DistanceTensor = std::vector<std::vector<std::vector<double>>>;

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

  void calculateFrame(DistributionFunctions &df,
                      const AnalysisSettings &settings) const override;

  static void compute(const Cell &cell, double cutoff_sq,
                      const std::vector<std::vector<double>> &bond_cutoffs_sq,
                      bool ignore_periodic_self_interactions,
                      DistanceTensor &out_distances, NeighborGraph &out_graph);
};

} // namespace calculators
