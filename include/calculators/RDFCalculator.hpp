/**
 * @file RDFCalculator.hpp
 * @brief Radial distribution function (RDF) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "DistributionFunctions.hpp"
#include <map>
#include <string>

/**
 * @class RDFCalculator
 * @brief Computes the Radial Distribution Function (RDF), also known as the
 * pair correlation function `g(r)`.
 *
 * It calculates the local particle density around a reference atom as a
 * function of distance `r`. From `g(r)`, it also calculates:L
 * - `G(r)`: The reduced pair distribution function, `G(r) = 4 * \pi * r *
 * \rho_0 * (g(r) - 1)`.
 * - `J(r)`: The radial distribution function related to the coordination
 * number, `J(r) = 4 * \pi * r^2 * \rho_0 * g(r)`.
 */
class RDFCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "g(r), J(r), G(r)"; }
  std::string getShortName() const override { return "RDF"; }
  std::string getGroup() const override { return "Radial"; }
  std::string getDescription() const override {
    return "Computes the Radial Distribution Function g_r, J_r, and G_r.";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(DistributionFunctions &df,
                      const AnalysisSettings &settings) const override;

  /**
   * @brief Calculates the J(r), g(r), and G(r) histograms for the given cell
   * and pair distances.
   *
   * @param cell The simulation Cell containing atoms and lattice information.
   * @param neighbors The StructureAnalyzer containing pre-computed pair
   * distances.
   * @param ashcroft_weights The Ashcroft-Langreth scattering weights for
   * assembling the total `g(r)` from partials.
   * @param r_max The maximum radius to calculate the RDF up to.
   * @param r_bin_width The histogram bin width in Angstroms.
   * @return A map containing the `J(r)`, `g(r)`, and `G(r)` histograms, both
   * partial and total.
   */
  static std::map<std::string, Histogram>
  calculate(const Cell &cell, const StructureAnalyzer *neighbors,
            const std::map<std::string, double> &ashcroft_weights, double r_max,
            double r_bin_width);
};
