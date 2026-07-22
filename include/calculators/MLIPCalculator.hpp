/**
 * @file MLIPCalculator.hpp
 * @brief Calculator for machine learning interatomic potential energy, forces, and stress.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "mlip/MLIPInterface.hpp"

#include <string>

namespace correlation::calculators {

/**
 * @class MLIPCalculator
 * @brief Computes atomic energy, forces, and stress using ORB-v3 / ORB-mol model family.
 */
class MLIPCalculator : public BaseCalculator {
public:
  MLIPCalculator() = default;

  [[nodiscard]] std::string getName() const override { return "ML Interatomic Potential (ORB-v3)"; }
  [[nodiscard]] std::string getShortName() const override { return "MLIP"; }
  [[nodiscard]] std::string getGroup() const override { return "Machine Learning"; }
  [[nodiscard]] std::string getDescription() const override {
    return "Evaluates atomic energy, forces, and stress using ORB-v3 / ORB-mol model family.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Evaluates potential properties for a single unit cell using an MLIP engine.
   * @param[in] cell The atomic unit cell.
   * @param[in] model Optional custom MLIPInterface engine pointer.
   * @return Generated MLIPOutput struct.
   */
  static correlation::mlip::MLIPOutput calculate(const correlation::core::Cell &cell,
                                                 const correlation::mlip::MLIPInterface *model = nullptr);
};

} // namespace correlation::calculators
