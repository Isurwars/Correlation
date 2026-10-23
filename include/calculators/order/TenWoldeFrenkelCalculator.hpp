/**
 * @file TenWoldeFrenkelCalculator.hpp
 * @brief Ten Wolde-Frenkel (TWF) bond coherence and solid-like clustering calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "calculators/BaseCalculator.hpp"
#include "calculators/order/SteinhardtCalculator.hpp"

#include <complex>
#include <span>
#include <vector>

namespace correlation::calculators {

/**
 * @struct TenWoldeFrenkelParams
 * @brief Configuration parameters for Ten Wolde-Frenkel bond coherence order parameter.
 */
struct TenWoldeFrenkelParams {
  real_t chi_cut = static_cast<real_t>(0.7); ///< Minimum bond correlation threshold.
  size_t nconn_min = 7; ///< Minimum coherent bonds for solid-like classification.
};

/**
 * @class TenWoldeFrenkelCalculator
 * @brief Computes Ten Wolde-Frenkel mutual bond coherence (N_conn) and crystalline nucleus
 * fractions.
 */
class TenWoldeFrenkelCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string_view getName() const override { return "Ten Wolde-Frenkel"; }
  [[nodiscard]] std::string_view getShortName() const override { return "TWF"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Structural"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Computes Ten Wolde-Frenkel mutual bond coherence (N_conn) and solid-like fractions.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Computes per-atom coherent neighbor count (N_conn).
   *
   * @param[in] cell The periodic cell.
   * @param[in] neighbors Structural analyzer containing the neighbor graph.
   * @param[in] params Coherence threshold parameters.
   * @param[in] steinhardt_params Optional pre-computed Steinhardt params with normalized q6m.
   * @return Vector of N_conn per atom.
   */
  [[nodiscard]] static std::vector<real_t>
  calculatePerAtom(const correlation::core::Cell &cell,
                   const correlation::analysis::StructureAnalyzer *neighbors,
                   TenWoldeFrenkelParams params = {},
                   const SteinhardtCalculator::SteinhardtParams *steinhardt_params = nullptr);

  /**
   * @brief Computes TWF N_conn distribution histogram.
   *
   * @param[in] cell The periodic cell.
   * @param[in] neighbors Structural analyzer containing the neighbor graph.
   * @param[in] params Coherence threshold parameters.
   * @param[in] steinhardt_params Optional pre-computed Steinhardt params.
   * @return Histogram of N_conn distribution.
   */
  [[nodiscard]] static correlation::analysis::Histogram
  calculate(const correlation::core::Cell &cell,
            const correlation::analysis::StructureAnalyzer *neighbors,
            TenWoldeFrenkelParams params = {},
            const SteinhardtCalculator::SteinhardtParams *steinhardt_params = nullptr);
};

} // namespace correlation::calculators
