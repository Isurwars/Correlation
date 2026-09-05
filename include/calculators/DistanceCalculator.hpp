/**
 * @file DistanceCalculator.hpp
 * @brief Pairwise distance calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/AnalysisTypes.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include <vector>

namespace correlation::calculators {

/** @brief Type definition for a 3D tensor of raw distance histogram bins [e1][e2][bin_index]. */
using RawHistogramTensor = std::vector<std::vector<std::vector<real_t>>>;

/**
 * @struct DistanceCalculationConfig
 * @brief Parameters for optional inline distance histogram accumulation.
 */
struct DistanceCalculationConfig {
  real_t r_max = 0.0;
  real_t r_bin_width = 0.0;
  size_t num_bins = 0;
};

/**
 * @class DistanceCalculator
 * @brief High-performance calculator for pairwise atomic distances.
 */
class DistanceCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string_view getName() const override { return "Distance"; }
  [[nodiscard]] std::string_view getShortName() const override { return "Distance"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Structural"; }
  [[nodiscard]] std::string_view getDescription() const override { return "Computes all unique 2-body distances."; }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief High-performance computation of the neighbor graph and optional inline raw distance histograms.
   *
   * @param cell The periodic cell container.
   * @param cutoff_sq Squared cutoff radius for neighbor search.
   * @param bond_cutoffs Matrix of minimum and maximum squared bond cutoffs per element pair.
   * @param ignore_periodic_self_interactions Flag to prevent self-images.
   * @param out_graph Graph object to be populated with adjacency data.
   * @param out_histograms Optional 3D tensor populated with raw pairwise distance histogram bins [e1][e2][bin].
   * @param hist_config Optional histogram grid configuration (r_max, r_bin_width, num_bins).
   */
  static void compute(const correlation::core::Cell &cell, real_t cutoff_sq,
                      const correlation::analysis::BondCutoffMatrix &bond_cutoffs,
                      bool ignore_periodic_self_interactions, correlation::core::NeighborGraph &out_graph,
                      RawHistogramTensor *out_histograms = nullptr, DistanceCalculationConfig hist_config = {});
};

} // namespace correlation::calculators
