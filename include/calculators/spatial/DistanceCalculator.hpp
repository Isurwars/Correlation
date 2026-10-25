/**
 * @file DistanceCalculator.hpp
 * @brief Pairwise distance calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "calculators/BaseCalculator.hpp"
#include "calculators/spatial/CellListTraversal.hpp"
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
  [[nodiscard]] std::string_view getName() const override { return "Pair Distances"; }
  [[nodiscard]] std::string_view getShortName() const override { return "Distance"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Structural"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Computes all unique 2-body pairwise distances.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  /**
   * @brief Dispatches the calculation for a single configuration frame.
   *
   * @note Currently a no-op as DistanceCalculator is invoked directly during StructureAnalyzer
   * construction.
   * @param[in,out] dists Distribution functions container.
   * @param[in] settings Current analysis configuration settings.
   */
  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief High-performance computation of the neighbor graph and optional inline raw distance
   * histograms.
   *
   * @param[in] cell The periodic cell container.
   * @param[in] cutoff_sq Squared cutoff radius for neighbor search.
   * @param[in] bond_cutoffs Matrix of minimum and maximum squared bond cutoffs per element pair.
   * @param[in] ignore_periodic_self_interactions Flag to prevent self-images.
   * @param[out] out_graph Graph object to be populated with adjacency data.
   * @param[in,out] out_histograms Optional 3D tensor populated with raw pairwise distance histogram
   * bins [e1][e2][bin].
   * @param[in] hist_config Optional histogram grid configuration (r_max, r_bin_width, num_bins).
   * @note Pairwise distance histogram normalization contract:
   * Same-species pairs (A == B) are accumulated once in out_histograms[A][A].
   * Cross-species pairs (A != B) are accumulated symmetrically in both out_histograms[A][B] and
   * out_histograms[B][A]. Downstream distribution calculators (e.g. RDFCalculator) scale [A][A]
   * bins by 2.0 to account for A_1 -> A_2 and A_2 -> A_1 interactions.
   * @throws std::invalid_argument If @p cutoff_sq is negative or if bond cutoffs are invalid.
   */
  static void compute(const correlation::core::Cell &cell, real_t cutoff_sq,
                      const correlation::analysis::BondCutoffMatrix &bond_cutoffs,
                      bool ignore_periodic_self_interactions,
                      correlation::core::NeighborGraph &out_graph,
                      RawHistogramTensor *out_histograms = nullptr,
                      DistanceCalculationConfig hist_config = {});

  /**
   * @brief High-performance parallel traversal of all unique atomic pairs within @p r_max.
   *
   * Traverses pairs (i, j) within cutoff distance @p r_max using bin-sorted cell lists
   * and SIMD distance filtering.
   *
   * @tparam Visitor Invocable with signature:
   *   void(size_t atom_i, size_t atom_j, real_t d_sq, const correlation::math::Vector3<real_t>
   * &r_ij)
   * @param[in] cell The periodic cell container.
   * @param[in] r_max Maximum interaction distance.
   * @param[in] ignore_periodic_self_interactions If true, ignores self-interactions in periodic
   * images.
   * @param[in] visitor Callable executed for each pair within cutoff.
   */
  template <typename Visitor>
  static void forEachPairWithin(const correlation::core::Cell &cell, real_t r_max,
                                bool ignore_periodic_self_interactions, Visitor &&visitor) {
    detail::forEachPairWithin(cell, r_max, ignore_periodic_self_interactions,
                              std::forward<Visitor>(visitor));
  }
};

} // namespace correlation::calculators
