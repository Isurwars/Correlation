/**
 * @file IFrameCalculator.hpp
 * @brief Abstract interface for per-frame spatial correlation calculators.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "calculators/contracts/ICalculator.hpp"

namespace correlation::analysis {
class DistributionFunctions;
struct AnalysisSettings;
} // namespace correlation::analysis

namespace correlation::calculators {

/**
 * @brief Segregated contract for calculators computing per-frame structural or spatial properties.
 */
class IFrameCalculator : public virtual ICalculator {
public:
  ~IFrameCalculator() override = default;

  /**
   * @brief Calculates properties for a single structure frame.
   * @param[in,out] dists DistributionFunctions container to accumulate partial and total distributions.
   * @param[in] settings Configuration parameters for analysis execution.
   */
  virtual void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                              const correlation::analysis::AnalysisSettings &settings) const = 0;
};

} // namespace correlation::calculators
