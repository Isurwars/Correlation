/**
 * @file RDCalculator.hpp
 * @brief Ring distribution calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "calculators/BaseCalculator.hpp"

namespace correlation::calculators {

/**
 * @class RDCalculator
 * @brief Computes the Ring Distribution (RD).
 */
class RDCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string_view getName() const override { return "RD"; }
  [[nodiscard]] std::string_view getShortName() const override { return "RD"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Rings"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Computes the Ring Distribution (RD).";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  /**
   * @brief Dispatches the calculation for a single configuration frame.
   *
   * @param[in,out] dists Distribution functions container to append results to.
   * @param[in] settings Current analysis configuration settings.
   */
  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @struct RDParams
   * @brief Configuration parameters for Ring Distribution calculation.
   */
  struct RDParams {
    size_t max_ring_size = 8; ///< Maximum ring size to search for.
    correlation::analysis::RingType ring_type =
        correlation::analysis::RingType::King; ///< Ring criterion (King or Franzblau).
    correlation::analysis::RingProjectionMode projection_mode =
        correlation::analysis::RingProjectionMode::Direct; ///< Projection or filtering strategy.
    std::string network_former;   ///< Network-former element symbol (e.g. "Si").
    std::string bridging_element; ///< Bridging element symbol (e.g. "O").
    bool report_polyhedra_size =
        true; ///< If true in bridged/alternating mode, reports size in former polyhedra.
  };

  /**
   * @brief High-performance computation of the Ring Distribution (RD).
   *
   * @param[in] graph The pre-computed neighbor graph.
   * @param[in] max_ring_size The maximum number of atoms in a single ring.
   * @param[in] ring_type Ring criterion (King or Franzblau).
   * @return A histogram representing the ring size distribution.
   * @throws std::invalid_argument If @p max_ring_size is less than 3.
   */
  static correlation::analysis::Histogram
  calculate(const correlation::core::NeighborGraph &graph, size_t max_ring_size,
            correlation::analysis::RingType ring_type = correlation::analysis::RingType::King);

  /**
   * @brief Comprehensive computation of Ring Distribution with network projections and
   * normalization.
   *
   * @param[in] graph The pre-computed neighbor graph.
   * @param[in] cell The atomic cell containing elements and geometry.
   * @param[in] params Ring calculation and network filtering parameters.
   * @return A histogram containing normalized fractions, raw counts, and per-former statistics.
   * @throws std::invalid_argument If max_ring_size is less than 3.
   */
  static correlation::analysis::Histogram calculate(const correlation::core::NeighborGraph &graph,
                                                    const correlation::core::Cell &cell,
                                                    const RDParams &params);
};

} // namespace correlation::calculators
