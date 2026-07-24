/**
 * @file VoronoiCalculator.hpp
 * @brief Voronoi Tessellation order parameter and metric calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "core/Cell.hpp"
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace correlation::calculators {

/**
 * @class VoronoiCalculator
 * @brief Computes Voronoi cell volumes, coordination numbers, sphericity, and polyhedral signatures.
 *
 * Space is partitioned around atoms using Voronoi cells, from which physical properties
 * and coordination topologies are extracted without the need for predefined cutoffs.
 */
class VoronoiCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string getName() const override { return "Voronoi Tessellation"; }
  [[nodiscard]] std::string getShortName() const override { return "Voronoi"; }
  [[nodiscard]] std::string getGroup() const override { return "Structural"; }
  [[nodiscard]] std::string getDescription() const override {
    return "Computes Voronoi cell volumes, coordination numbers, sphericity, and polyhedral signatures.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Computes Voronoi tessellation and extracts metric histograms.
   *
   * @param cell The periodic simulation cell.
   * @param neighbors Optional structure analyzer for API consistency.
   * @return A map of histogram names to their respective Histogram data.
   */
  static std::map<std::string, correlation::analysis::Histogram>
  calculate(const correlation::core::Cell &cell, const correlation::analysis::StructureAnalyzer *neighbors = nullptr);

private:
  /**
   * @struct BinRange
   * @brief Describes the binning range for a Voronoi metric histogram.
   */
  struct BinRange {
    real_t bin_width{0.0}; ///< Width of individual histogram bins.
    size_t num_bins{0};    ///< Total number of histogram bins.
    real_t range_min{0.0}; ///< Lower bound of histogram range.
    real_t range_max{0.0}; ///< Upper bound of histogram range.
  };

  /**
   * @struct CellData
   * @brief Aggregated metric data extracted from all Voronoi cells in a configuration.
   */
  struct CellData {
    std::vector<real_t> volumes;           ///< Per-atom Voronoi cell volumes (Angstrom^3).
    std::vector<real_t> sphericities;      ///< Per-atom cell sphericities (dimensionless).
    std::vector<int> coordination_numbers; ///< Per-atom Voronoi face coordination counts.
    std::vector<std::string> signatures;   ///< Per-atom polyhedral index signatures <n3,n4,n5,n6>.
  };

  /// Build the Voronoi tessellation and extract per-atom cell data.
  static CellData computeVoronoiCells(const correlation::core::Cell &cell);

  /// Build the signature frequency map and return sorted signatures + description.
  static std::pair<std::vector<std::pair<std::string, int>>, std::string>
  buildSignatureMap(const std::vector<std::string> &signatures);

  /// Create an initialized histogram with bins and empty partials.
  static correlation::analysis::Histogram makeHistogram(const std::string &title, const std::string &x_label,
                                                        const std::string &y_label, const std::string &x_unit,
                                                        const std::string &y_unit, const std::string &description,
                                                        const std::string &file_suffix, const std::vector<real_t> &bins,
                                                        const std::vector<std::string> &element_symbols);

  /// Populate histogram bins from per-atom numeric data.
  static void populateHistogram(correlation::analysis::Histogram &hist, const BinRange &range,
                                const std::vector<real_t> &values, const std::vector<correlation::core::Atom> &atoms);
};

} // namespace correlation::calculators
