/**
 * @file CNCalculator.hpp
 * @brief Coordination number (CN) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

namespace correlation::calculators {

/**
 * @class CNCalculator
 * @brief Computes the Coordination Number (CN) distribution.
 */
class CNCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string_view getName() const override { return "CN"; }
  [[nodiscard]] std::string_view getShortName() const override { return "CN"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Structural"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Computes the Coordination Number (CN).";
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
   * @brief High-performance computation of the Coordination Number (CN).
   *
   * @param[in] cell The periodic cell.
   * @param[in] neighbors Structural analyzer containing the neighbor graph.
   * @return A histogram representing the CN distribution (count vs coordination).
   * @throws std::logic_error If @p neighbors is nullptr.
   */
  static correlation::analysis::Histogram
  calculate(const correlation::core::Cell &cell,
            const correlation::analysis::StructureAnalyzer *neighbors);
};

} // namespace correlation::calculators
