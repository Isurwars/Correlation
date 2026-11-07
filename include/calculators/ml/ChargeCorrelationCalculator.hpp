/**
 * @file ChargeCorrelationCalculator.hpp
 * @brief Charge-weighted pair distribution function g_qq(r) and screening integral calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "calculators/BaseCalculator.hpp"
#include "core/Cell.hpp"
#include "mlip/MLIPInterface.hpp"

#include <span>
#include <vector>

namespace correlation::calculators {

/**
 * @struct ChargeCorrelationResult
 * @brief Output container for charge-weighted correlation and screening distributions.
 */
struct alignas(64) ChargeCorrelationResult {
  std::vector<real_t> r_bins;          /**< Radial distance centers in Angstroms. */
  std::vector<real_t> g_qq;            /**< Charge-weighted pair correlation function g_qq(r). */
  std::vector<real_t> charge_integral; /**< Running charge neutrality screening integral Q(r). */
  std::vector<uint64_t> pair_counts;   /**< Total pair observations in each radial bin. */
};

/**
 * @class ChargeCorrelationCalculator
 * @brief Computes charge-weighted pair distribution function g_qq(r) from atomic partial charges.
 */
class ChargeCorrelationCalculator : public BaseCalculator {
public:
  ChargeCorrelationCalculator() = default;

  [[nodiscard]] std::string_view getName() const override { return "Charge Pair Distribution"; }
  [[nodiscard]] std::string_view getShortName() const override { return "ChargePDF"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Machine Learning"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Evaluates charge-weighted pair distribution g_qq(r) and cumulative screening integral Q(r).";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Computes charge correlation from a Cell and atomic charges.
   * @param[in] cell Atomic configuration cell.
   * @param[in] charges Per-atom partial / Born charges (elementary charges e).
   * @param[in] r_max Radial cutoff distance in Angstroms (default: 10.0).
   * @param[in] r_bin_width Radial histogram bin width in Angstroms (default: 0.05).
   * @return Calculated ChargeCorrelationResult.
   */
  [[nodiscard]] static ChargeCorrelationResult
  calculate(const correlation::core::Cell &cell, std::span<const real_t> charges,
            real_t r_max = static_cast<real_t>(10.0),
            real_t r_bin_width = static_cast<real_t>(0.05));

  /**
   * @brief Computes charge correlation directly from an MLIP evaluation output.
   * @param[in] cell Atomic configuration cell.
   * @param[in] mlip_out Evaluated MLIPOutput containing charges.
   * @param[in] r_max Radial cutoff distance in Angstroms.
   * @param[in] r_bin_width Radial histogram bin width in Angstroms.
   * @return Calculated ChargeCorrelationResult.
   */
  [[nodiscard]] static ChargeCorrelationResult calculateFromMLIP(
      const correlation::core::Cell &cell, const correlation::mlip::MLIPOutput &mlip_out,
      real_t r_max = static_cast<real_t>(10.0), real_t r_bin_width = static_cast<real_t>(0.05));
};

} // namespace correlation::calculators
