/**
 * @file ISpatialCalculator.hpp
 * @brief Contract for single-frame spatial calculators.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "calculators/BaseCalculator.hpp"

namespace correlation::calculators {

/**
 * @brief Abstract contract for calculators operating on per-frame spatial coordinates and cells.
 * Enforces isFrameCalculator() = true and provides a pure virtual calculateFrame contract.
 */
class ISpatialCalculator : public BaseCalculator {
public:
  ~ISpatialCalculator() override = default;

  [[nodiscard]] bool isFrameCalculator() const noexcept override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const noexcept override { return false; }

  /**
   * @brief Calculate per-frame properties.
   *        Called concurrently on distinct DistributionFunctions objects.
   *        Must be thread-safe (const).
   * @param[in,out] dists The DistributionFunctions object to store the results.
   * @param[in] settings The analysis settings and parameters.
   */
  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override = 0;
};

} // namespace correlation::calculators
