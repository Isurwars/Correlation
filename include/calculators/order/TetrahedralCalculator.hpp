/**
 * @file TetrahedralCalculator.hpp
 * @brief Chau-Hardwick / Errington-Debenedetti tetrahedral order parameter calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/BaseCalculator.hpp"

#include <vector>

namespace correlation::calculators {

/**
 * @struct TetrahedralParams
 * @brief Configuration parameters for tetrahedral order parameter calculation.
 */
struct TetrahedralParams {
  size_t num_bins = 50;                         ///< Number of bins for 1D histogram.
  real_t min_domain = static_cast<real_t>(0.0); ///< Lower domain boundary.
  real_t max_domain = static_cast<real_t>(1.0); ///< Upper domain boundary.
};

/**
 * @class TetrahedralCalculator
 * @brief Computes Chau-Hardwick tetrahedral order parameter (q_tet) quantifying 4-fold angular
 * coordination.
 */
class TetrahedralCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string_view getName() const override { return "Tetrahedral"; }
  [[nodiscard]] std::string_view getShortName() const override { return "q_tet"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Structural"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Computes Chau-Hardwick tetrahedral order parameter (q_tet) for 4-fold networks.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Computes per-atom tetrahedral order parameter (q_tet).
   *
   * @param[in] cell The periodic cell.
   * @param[in] neighbors Optional structural analyzer containing the neighbor graph.
   * @param[in] params Calculation parameters.
   * @return Vector of q_tet per atom.
   */
  [[nodiscard]] static std::vector<real_t>
  calculatePerAtom(const correlation::core::Cell &cell,
                   const correlation::analysis::StructureAnalyzer *neighbors = nullptr,
                   TetrahedralParams params = {});

  /**
   * @brief Computes tetrahedral order parameter distribution histogram.
   *
   * @param[in] cell The periodic cell.
   * @param[in] neighbors Optional structural analyzer containing the neighbor graph.
   * @param[in] params Calculation parameters.
   * @return Histogram of q_tet distribution.
   */
  [[nodiscard]] static correlation::analysis::Histogram
  calculate(const correlation::core::Cell &cell,
            const correlation::analysis::StructureAnalyzer *neighbors = nullptr,
            TetrahedralParams params = {});
};

} // namespace correlation::calculators
