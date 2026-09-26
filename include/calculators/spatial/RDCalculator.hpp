/**
 * @file RDCalculator.hpp
 * @brief Ring distribution calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "calculators/BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

namespace correlation::calculators {

/**
 * @class RDCalculator
 * @brief Computes the Ring Distribution (RD).
 */
class RDCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string_view getName() const override { return "RD"; }
  [[nodiscard]] std::string_view getShortName() const override { return "RD"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Rings"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Computes the Ring Distribution (RD).";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  /**
   * @brief Dispatches the calculation for a single configuration frame.
   *
   * @param[in,out] dists Distribution functions container to append results to.
   * @param[in] settings Current analysis configuration settings.
   */
  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief High-performance computation of the Ring Distribution (RD).
   *
   * @param[in] graph The pre-computed neighbor graph.
   * @param[in] max_ring_size The maximum number of atoms in a single ring.
   * @return A histogram representing the ring size distribution.
   * @throws std::invalid_argument If @p max_ring_size is less than 3.
   */
  static correlation::analysis::Histogram calculate(const correlation::core::NeighborGraph &graph,
                                                    size_t max_ring_size);
};

} // namespace correlation::calculators
