// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "DistributionFunctions.hpp"
#include <map>
#include <string>

/**
 * @class RDFCalculator
 * @brief Computes the Radial Distribution Function (RDF), also known as the
 * pair correlation function `g(r)`.
 *
 * It calculates the local particle density around a reference atom as a
 * function of distance `r`. From `g(r)`, it also calculates:
 * - `G(r)`: The reduced pair distribution function, `G(r) = 4 * \pi * r *
 * \rho_0 * (g(r) - 1)`.
 * - `J(r)`: The radial distribution function related to the coordination
 * number, `J(r) = 4 * \pi * r^2 * \rho_0 * g(r)`.
 */
class RDFCalculator {
public:
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
