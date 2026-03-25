// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "BaseCalculator.hpp"
#include "DistributionFunctions.hpp"

/**
 * @class StructureFactorCalculator
 * @brief Computes the static structure factor S(Q) using plane wave summation
 *        over reciprocal lattice vectors for periodic systems.
 *
 * For each reciprocal lattice vector q, computes:
 *   rho(q) = sum_j exp(i q . r_j)
 *   S(q)   = |rho(q)|^2 / N
 *
 * Results are averaged over q-vectors with similar magnitudes (|q| shells).
 */
class StructureFactorCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "S(K)"; }
  std::string getShortName() const override { return "S_K"; }
  std::string getGroup() const override { return "Radial"; }
  std::string getDescription() const override {
    return "Computes the static structure factor S(K) via plane wave summation "
           "over reciprocal lattice vectors (for periodic systems).";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(DistributionFunctions &df,
                      const AnalysisSettings &settings) const override;
};
