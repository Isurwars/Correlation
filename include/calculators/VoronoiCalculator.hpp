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
  calculate(const correlation::core::Cell &cell,
            const correlation::analysis::StructureAnalyzer *neighbors = nullptr);

private:
  /// Describes the binning range for a histogram.
  struct BinRange {
    double bin_width;
    size_t num_bins;
    double range_min;
    double range_max;
  };

  /// Data extracted from a single Voronoi cell for each atom.
  struct CellData {
    std::vector<double> volumes;
    std::vector<double> sphericities;
    std::vector<int> coordination_numbers;
    std::vector<std::string> signatures;
  };

  /// Build the Voronoi tessellation and extract per-atom cell data.
  static CellData computeVoronoiCells(const correlation::core::Cell &cell);

  /// Build the signature frequency map and return sorted signatures + description.
  static std::pair<std::vector<std::pair<std::string, int>>, std::string>
  buildSignatureMap(const std::vector<std::string> &signatures);

  /// Create an initialized histogram with bins and empty partials.
  static correlation::analysis::Histogram
  makeHistogram(const std::string &title, const std::string &x_label,
                const std::string &y_label, const std::string &x_unit,
                const std::string &y_unit, const std::string &description,
                const std::string &file_suffix, const std::vector<double> &bins,
                const std::vector<std::string> &element_symbols);

  /// Populate histogram bins from per-atom numeric data.
  static void populateHistogram(correlation::analysis::Histogram &hist,
                                const BinRange &range,
                                const std::vector<double> &values,
                                const std::vector<correlation::core::Atom> &atoms);
};

} // namespace correlation::calculators
