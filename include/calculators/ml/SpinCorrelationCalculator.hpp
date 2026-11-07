/**
 * @file SpinCorrelationCalculator.hpp
 * @brief Spatial spin-spin and magnetic moment correlation calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "calculators/BaseCalculator.hpp"
#include "core/Cell.hpp"
#include "math/LinearAlgebra.hpp"
#include "mlip/MLIPInterface.hpp"

#include <span>
#include <vector>

namespace correlation::calculators {

/**
 * @struct SpinCorrelationResult
 * @brief Output container for spatial spin and magnetic correlation functions.
 */
struct alignas(64) SpinCorrelationResult {
  std::vector<real_t> r_bins;           /**< Radial distance centers in Angstroms. */
  std::vector<real_t> spin_correlation; /**< Normalized spin alignment <S_i · S_j>(r) in [-1, 1]. */
  std::vector<real_t>
      magnetic_moment_corr; /**< Unnormalized moment correlation <m_i · m_j>(r) (mu_B^2). */
  std::vector<uint64_t> pair_counts; /**< Total pair observations in each radial bin. */
};

/**
 * @class SpinCorrelationCalculator
 * @brief Computes spatial spin-spin correlation <S_i · S_j>(r) from site magnetic moments.
 */
class SpinCorrelationCalculator : public BaseCalculator {
public:
  SpinCorrelationCalculator() = default;

  [[nodiscard]] std::string_view getName() const override { return "Spin Correlation"; }
  [[nodiscard]] std::string_view getShortName() const override { return "SpinCorr"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Machine Learning"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Evaluates spatial spin-spin correlation <S_i · S_j>(r) and magnetic moment RDF.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Computes spin-spin correlation from a Cell and magnetic moment vectors.
   * @param[in] cell Atomic configuration cell.
   * @param[in] magmoms Per-atom magnetic moment vectors in Bohr magnetons (mu_B).
   * @param[in] r_max Radial cutoff distance in Angstroms (default: 10.0).
   * @param[in] r_bin_width Radial histogram bin width in Angstroms (default: 0.05).
   * @return Calculated SpinCorrelationResult.
   */
  [[nodiscard]] static SpinCorrelationResult
  calculate(const correlation::core::Cell &cell,
            std::span<const correlation::math::Vector3<real_t>> magmoms,
            real_t r_max = static_cast<real_t>(10.0),
            real_t r_bin_width = static_cast<real_t>(0.05));

  /**
   * @brief Computes spin-spin correlation directly from an MLIP evaluation output.
   * @param[in] cell Atomic configuration cell.
   * @param[in] mlip_out Evaluated MLIPOutput containing magmoms.
   * @param[in] r_max Radial cutoff distance in Angstroms.
   * @param[in] r_bin_width Radial histogram bin width in Angstroms.
   * @return Calculated SpinCorrelationResult.
   */
  [[nodiscard]] static SpinCorrelationResult calculateFromMLIP(
      const correlation::core::Cell &cell, const correlation::mlip::MLIPOutput &mlip_out,
      real_t r_max = static_cast<real_t>(10.0), real_t r_bin_width = static_cast<real_t>(0.05));
};

} // namespace correlation::calculators
