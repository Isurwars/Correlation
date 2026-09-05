/**
 * @file ClusterCalculator.hpp
 * @brief Calculator for cluster analysis using Union-Find on the neighbor graph.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "calculators/BaseCalculator.hpp"

namespace correlation::calculators {

/**
 * @class ClusterCalculator
 * @brief Identifies connected components in the atomic structure.
 *
 * Uses a disjoint-set (Union-Find) algorithm to traverse the neighbor graph
 * and group atoms into clusters based on the connectivity radius. It then
 * outputs the cluster size distribution.
 */
class ClusterCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string_view getName() const override { return "Cluster Analysis"; }
  [[nodiscard]] std::string_view getShortName() const override { return "Clusters"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Topology"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Identifies connected components (clusters) using the neighbor graph.";
  }
  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;
};

} // namespace correlation::calculators
