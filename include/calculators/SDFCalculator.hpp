/**
 * @file SDFCalculator.hpp
 * @brief Spatial Distribution Function (SDF) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

#include <string>

namespace correlation::calculators {

/**
 * @class SDFCalculator
 * @brief Computes the Spatial Distribution Function (SDF) in 3D.
 *
 * It calculates the 3D density of particles. The results are stored
 * as flattened 3D grids in the DistributionFunctions::Histogram to
 * enable standard accumulation and averaging.
 */
class SDFCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "Spatial Distribution Function (3D)"; }
  std::string getShortName() const override { return "SDF"; }
  std::string getGroup() const override { return "Spatial"; }
  std::string getDescription() const override {
    return "Computes the 3D Spatial Distribution Function (SDF) density grid.";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &df,
                      const correlation::analysis::AnalysisSettings &settings) const override;
};

} // namespace correlation::calculators
