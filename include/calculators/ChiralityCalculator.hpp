/**
 * @file ChiralityCalculator.hpp
 * @brief Chiral order parameter (COP) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

namespace correlation::calculators {

/**
 * @class ChiralityCalculator
 * @brief Computes the Chiral Order Parameter (COP) for all atoms.
 *
 * For each atom, COP quantifies the local mirror asymmetry (chirality)
 * of its local motif based on its 3 nearest neighbors.
 */
class ChiralityCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string getName() const override { return "Chiral Order Parameter"; }
  [[nodiscard]] std::string getShortName() const override { return "COP"; }
  [[nodiscard]] std::string getGroup() const override { return "Structural"; }
  [[nodiscard]] std::string getDescription() const override {
    return "Computes the Chiral Order Parameter (COP) distribution.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Computes the Chiral Order Parameter for a single atom.
   *
   * @param atom_idx Index of the central atom.
   * @param cell The periodic cell.
   * @param neighbors Structural analyzer containing the neighbor graph.
   * @return The normalized scalar triple product chirality value.
   */
  static double computeSingleAtomChirality(size_t atom_idx,
                                           const correlation::core::Cell &cell,
                                           const correlation::analysis::StructureAnalyzer *neighbors);

  /**
   * @brief Computes the Chiral Order Parameter distribution for all atoms in the cell.
   *
   * @param cell The periodic cell.
   * @param neighbors Structural analyzer containing the neighbor graph.
   * @return A histogram of chirality values in the range [-1.0, 1.0].
   */
  static correlation::analysis::Histogram
  calculate(const correlation::core::Cell &cell,
            const correlation::analysis::StructureAnalyzer *neighbors);
};

} // namespace correlation::calculators
