/**
 * @file HyperuniformityCalculator.hpp
 * @brief Calculator for local number variance and hyperuniformity index.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

#include <cstddef>
#include <map>
#include <string>

namespace correlation::calculators {

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
  std::string getName() const override { return "σ²_N(R), χ_H(R)"; }
  std::string getShortName() const override { return "Hyperuniformity"; }
  std::string getGroup() const override { return "Advanced"; }
  std::string getDescription() const override { return "Computes local number variance and hyperuniformity index."; }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Computes σ²_N(R) and χ_H(R) histograms.
   *
   * @param cell The simulation cell.
   * @param params Sampling and binning parameters.
   * @return A map containing "sigma2_N" and "chi_H" histograms.
   */
  static std::map<std::string, correlation::analysis::Histogram> calculate(const correlation::core::Cell &cell,
                                                                           const HyperuniformityParams &params);
};

} // namespace correlation::calculators
