/**
 * @file SteinhardtCalculator.hpp
 * @brief Steinhardt bond-order parameter calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

#include <complex>

namespace correlation::calculators {

/**
 * @class SteinhardtCalculator
 * @brief Computes Steinhardt Bond-Orientational Parameters (Q4, Q6, W6).
 */
class SteinhardtCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "Steinhardt Parameter"; }
  std::string getShortName() const override { return "Steinhardt"; }
  std::string getGroup() const override { return "Structural"; }
  std::string getDescription() const override {
    return "Computes Steinhardt Bond-Orientational Parameters (Q4, Q6, W6).";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(
      correlation::analysis::DistributionFunctions &df,
      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief High-performance computation of Steinhardt parameters (Q4, Q6, W6).
   * 
   * @param cell The periodic cell.
   * @param neighbors Structural analyzer containing the neighbor graph.
   * @return A map of histograms for each requested parameter.
   */
  static std::map<std::string, correlation::analysis::Histogram>
  calculate(const correlation::core::Cell &cell,
            const correlation::analysis::StructureAnalyzer *neighbors);

  /**
   * @brief Computes the complex spherical harmonic @f$ Y_l^m(\theta, \phi) @f$.
   * 
   * @param l Degree of the harmonic.
   * @param m Order of the harmonic.
   * @param theta Inclination angle (radians).
   * @param phi Azimuthal angle (radians).
   * @return The complex value of the spherical harmonic.
   */
  static std::complex<double> sphericalHarmonic(int l, int m, double theta,
                                                double phi);

  /**
   * @brief Computes the Wigner 3-j symbol using the Racah formula.
   * 
   * @param j1, j2, j3 Angular momenta.
   * @param m1, m2, m3 Magnetic projections.
   * @return The scalar value of the 3-j symbol.
   */
  static double wigner3j(int j1, int j2, int j3, int m1, int m2, int m3);
};

} // namespace correlation::calculators
