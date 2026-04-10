/**
 * @file RDCalculator.hpp
 * @brief Reduced distribution function calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

namespace correlation::calculators {

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

  void calculateFrame(
      correlation::analysis::DistributionFunctions &df,
      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief High-performance computation of the Ring Distribution (RD).
   * 
   * @param graph The pre-computed neighbor graph.
   * @param max_ring_size The maximum number of atoms in a single ring.
   * @return A histogram representing the ring size distribution.
   */
  static correlation::analysis::Histogram
  calculate(const correlation::core::NeighborGraph &graph,
            size_t max_ring_size);
};

} // namespace correlation::calculators
