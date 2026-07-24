/**
 * @file LocalEntropyCalculator.hpp
 * @brief Local Entropic Fingerprint order parameter calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

namespace correlation::calculators {

/**
 * @brief Parameters for the Local Entropic Fingerprint (LEF) calculation.
 */
struct LocalEntropyParams {
  real_t cutoff = 5.0; ///< Integration cutoff radius Rc in Angstroms.
  real_t sigma = 0.2;  ///< Gaussian kernel standard deviation in Angstroms.
};

/**
 * @class LocalEntropyCalculator
 * @brief Computes the Local Entropic Fingerprint (LEF) for all atoms.
 *
 * The local entropy measures the degree of local order around each atom
 * based on its local pair correlation function.
 */
class LocalEntropyCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string getName() const override { return "Local Entropy"; }
  [[nodiscard]] std::string getShortName() const override { return "LEF"; }
  [[nodiscard]] std::string getGroup() const override { return "Structural"; }
  [[nodiscard]] std::string getDescription() const override {
    return "Computes the Local Entropic Fingerprint (LEF) distribution.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Computes the Local Entropic Fingerprint for all atoms in the cell.
   *
   * @param cell The periodic cell.
   * @param neighbors Structural analyzer containing the neighbor graph.
   * @param params Integration cutoff and Gaussian smoothing parameters.
   * @return A histogram of local entropy values.
   */
  static correlation::analysis::Histogram calculate(const correlation::core::Cell &cell,
                                                    const correlation::analysis::StructureAnalyzer *neighbors,
                                                    LocalEntropyParams params);
};

} // namespace correlation::calculators
