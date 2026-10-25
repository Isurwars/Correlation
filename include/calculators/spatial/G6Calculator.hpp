/**
 * @file G6Calculator.hpp
 * @brief Radial bond-orientational correlation function g6(r) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/BaseCalculator.hpp"
#include "calculators/order/SteinhardtCalculator.hpp"
#include "core/Cell.hpp"

#include <cmath>
#include <span>
#include <string_view>
#include <vector>

namespace correlation::calculators {

/**
 * @struct CorrelationLengthFit
 * @brief Results of exponential decay fitting for bond-orientational correlation g6(r).
 */
struct CorrelationLengthFit {
  real_t xi_6{0.0};       ///< Correlation length in Angstroms (> 0, or inf if plateau).
  real_t r_squared{0.0};  ///< Coefficient of determination R^2 of log-linear regression.
  bool is_plateau{false}; ///< True if crystalline plateau without exponential decay.
};

/**
 * @class G6Calculator
 * @brief Computes distance-resolved bond-orientational correlation function g6(r) and extracts
 * xi_6.
 */
class G6Calculator : public BaseCalculator {
public:
  [[nodiscard]] std::string_view getName() const override {
    return "Bond-Orientational Correlation (G6)";
  }
  [[nodiscard]] std::string_view getShortName() const override { return "G6"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Spatial"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Computes radial bond-orientational correlation g6(r) and extracts correlation length xi_6.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  /**
   * @brief Dispatches calculation for a single frame, storing histogram as "G6_r".
   *
   * @param[in,out] dists Distribution functions container.
   * @param[in] settings Active analysis settings containing r_max and r_bin_width.
   */
  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Fits the bond-orientational correlation length xi_6 from radial profile g6(r).
   *
   * Uses log-linear regression over envelope peaks: ln(g6) = a - r / xi_6.
   *
   * @param[in] r_bins Radial distance bins in Angstroms.
   * @param[in] g6_values Bond correlation function values g6(r).
   * @return Fitted correlation length parameters.
   */
  [[nodiscard]] static CorrelationLengthFit
  fitCorrelationLength(std::span<const real_t> r_bins, std::span<const real_t> g6_values) noexcept;

  /**
   * @brief High-performance evaluation of distance-resolved g6(r).
   *
   * @param[in] cell The periodic cell.
   * @param[in] neighbors Structural analyzer containing neighbor graph.
   * @param[in] r_max Maximum radial distance.
   * @param[in] r_bin_width Radial bin width in Angstroms.
   * @param[in] steinhardt_params Optional pre-computed Steinhardt params with normalized q6m.
   * @return 1D Histogram containing radial g6(r) profile and fit metadata.
   * @throws std::logic_error If neighbors is nullptr.
   * @throws std::invalid_argument If r_max <= 0 or r_bin_width <= 0.
   */
  [[nodiscard]] static correlation::analysis::Histogram
  calculate(const correlation::core::Cell &cell,
            const correlation::analysis::StructureAnalyzer *neighbors, real_t r_max,
            real_t r_bin_width,
            const SteinhardtCalculator::SteinhardtParams *steinhardt_params = nullptr);
};

} // namespace correlation::calculators
