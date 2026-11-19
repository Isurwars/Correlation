/**
 * @file HyperuniformityCalculator.hpp
 * @brief Calculator for local number variance and hyperuniformity index.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "calculators/BaseCalculator.hpp"

#include <cstddef>
#include <map>
#include <span>
#include <string>
#include <string_view>

namespace correlation::calculators {

/**
 * @brief Classification result for a hyperuniformity check.
 */
struct HyperuniformityVerdict {
  std::string_view label; ///< Human-readable regime classification.
  int severity_level{0}; ///< 0: Non-hyperuniform, 1: Effectively hyperuniform, 2: Near-HU, -1: N/A.
};

/**
 * @brief Parameters for the hyperuniformity calculation.
 *
 * Bundles the sampling and binning parameters into a single struct to
 * prevent accidental argument swapping (bugprone-easily-swappable-parameters).
 */
struct HyperuniformityParams {
  size_t num_samples; ///< Number of random sample points.
  real_t r_bin_width; ///< Bin width for R values (Å).
};

/**
 * @class HyperuniformityCalculator
 * @brief Computes the local number variance σ²_N(R) and hyperuniformity
 * index χ_H(R) by sampling random observation windows in the periodic box.
 *
 * For each window radius R (from 2.0 Å to L/2, where L is the minimum box
 * length), random points are placed in the periodic box and the number of
 * atoms N(R) within each spherical window is counted. The statistics are:
 *
 *   σ²_N(R) = ⟨N²(R)⟩ − ⟨N(R)⟩²
 *   χ_H(R)  = σ²_N(R) / ⟨N(R)⟩
 *
 * For a hyperuniform system, σ²_N(R) grows as R² (surface-dominated),
 * while for an uncorrelated (Poisson) system it grows as R³ (volume-dominated).
 */
class HyperuniformityCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string_view getName() const override { return "Hyperuniformity Diagnostics"; }
  [[nodiscard]] std::string_view getShortName() const override { return "Hyperuniformity"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Advanced"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Computes local number variance sigma2_N(R) and hyperuniformity index chi_H(R).";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Computes σ²_N(R) and χ_H(R) histograms.
   *
   * @param[in] cell The simulation cell.
   * @param[in] params Sampling and binning parameters.
   * @return A map containing "sigma2_N" and "chi_H" histograms.
   * @throws std::invalid_argument if bin width or num_samples is non-positive.
   */
  [[nodiscard]] static std::map<std::string, correlation::analysis::Histogram>
  calculate(const correlation::core::Cell &cell, const HyperuniformityParams &params);

  /**
   * @brief Extracts real-space hyperuniformity metrics (β and large-R χ_H) from σ²_N(R) and χ_H(R).
   *
   * @param[in] radii Window radii R (Å).
   * @param[in] sigma2_n Local number variance σ²_N(R).
   * @param[in] chi_h Hyperuniformity index χ_H(R) = σ²_N(R) / ⟨N(R)⟩.
   * @return Populated real-space fields of HyperuniformityDiagnostics.
   */
  [[nodiscard]] static correlation::analysis::HyperuniformityDiagnostics
  computeRealSpaceDiagnostics(std::span<const real_t> radii, std::span<const real_t> sigma2_n,
                              std::span<const real_t> chi_h) noexcept;

  /**
   * @brief Extracts reciprocal-space hyperuniformity metrics (H = S(0⁺)/S(q_peak) and α) from S(Q).
   *
   * @param[in] q_bins Reciprocal wavevector magnitudes q (Å⁻¹).
   * @param[in] sq_total Total static structure factor S(q).
   * @return Populated reciprocal-space fields of HyperuniformityDiagnostics.
   */
  [[nodiscard]] static correlation::analysis::HyperuniformityDiagnostics
  computeReciprocalDiagnostics(std::span<const real_t> q_bins,
                               std::span<const real_t> sq_total) noexcept;

  /**
   * @brief Classifies the hyperuniformity regime from real- and/or reciprocal-space diagnostics.
   *
   * @param[in] diag Computed scalar hyperuniformity diagnostics.
   * @return Regime label and severity level (1 = hyperuniform, 2 = near-HU, 0 = non-HU, -1 = N/A).
   */
  [[nodiscard]] static HyperuniformityVerdict
  classifyHyperuniformity(const correlation::analysis::HyperuniformityDiagnostics &diag) noexcept;
};

} // namespace correlation::calculators
