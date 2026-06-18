/**
 * @file SteinhardtCalculator.hpp
 * @brief Steinhardt bond-order parameter calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
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
  [[nodiscard]] std::string getName() const override { return "Steinhardt Parameter"; }
  [[nodiscard]] std::string getShortName() const override { return "Steinhardt"; }
  [[nodiscard]] std::string getGroup() const override { return "Structural"; }
  [[nodiscard]] std::string getDescription() const override {
    return "Computes Steinhardt Bond-Orientational Parameters (Q4, Q6, W6).";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief High-performance computation of Steinhardt parameters (Q4, Q6, W6).
   *
   * @param cell The periodic cell.
   * @param neighbors Structural analyzer containing the neighbor graph.
   * @return A map of histograms for each requested parameter.
   */
  static std::map<std::string, correlation::analysis::Histogram>
  calculate(const correlation::core::Cell &cell, const correlation::analysis::StructureAnalyzer *neighbors);

  /**
   * @brief Computes the complex spherical harmonic @f$ Y_l^m(\theta, \phi) @f$.
   *
   * @param degree Degree of the harmonic.
   * @param order Order of the harmonic.
   * @param theta Inclination angle (radians).
   * @param phi Azimuthal angle (radians).
   * @return The complex value of the spherical harmonic.
   */
  static std::complex<double> sphericalHarmonic(int degree, int order, double theta, double phi);

  /**
   * @brief Computes the Wigner 3-j symbol using the Racah formula.
   *
   * @param j_one, j_two, j_three Angular momenta.
   * @param m_one, m_two, m_three Magnetic projections.
   * @return The scalar value of the 3-j symbol.
   */
  static double wigner3j(int j_one, int j_two, int j_three, int m_one, int m_two, int m_three);
};

} // namespace correlation::calculators
