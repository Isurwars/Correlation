/**
 * @file RDCalculator.hpp
 * @brief Reduced distribution function calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "DistributionFunctions.hpp"

/**
 * @class RDCalculator
 * @brief Computes the Ring Distribution (RD).
 */
class RDCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "RD"; }
  std::string getShortName() const override { return "RD"; }
  std::string getGroup() const override { return "Rings"; }
  std::string getDescription() const override {
    return "Computes the Ring Distribution (RD).";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(DistributionFunctions &df,
                      const AnalysisSettings &settings) const override;

  static Histogram calculate(const NeighborGraph &graph, size_t max_ring_size);
};
