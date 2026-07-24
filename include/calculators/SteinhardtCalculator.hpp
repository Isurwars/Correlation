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
   * @struct SphericalAngles
   * @brief Spherical polar coordinates (theta, phi) for angular calculations.
   */
  struct SphericalAngles {
    real_t theta{0.0};
    real_t phi{0.0};
  };

  /**
   * @struct SingleAtomSteinhardt
   * @brief Holds calculated Steinhardt order parameters for an individual atom.
   */
  struct SingleAtomSteinhardt {
    real_t Q4{0.0};     ///< Q4 bond-orientational order parameter.
    real_t Q6{0.0};     ///< Q6 bond-orientational order parameter.
    real_t W6_hat{0.0}; ///< Normalized W6 bond-orientational order parameter.
  };

  /**
   * @struct SteinhardtParams
   * @brief Per-atom arrays of Steinhardt order parameters across an atomic configuration frame.
   */
  struct SteinhardtParams {
    std::vector<real_t> Q4;     ///< Per-atom Q4 values.
    std::vector<real_t> Q6;     ///< Per-atom Q6 values.
    std::vector<real_t> W6_hat; ///< Per-atom W6_hat values.
  };

  /**
   * @brief Computes the complex spherical harmonic @f$ Y_l^m(\theta, \phi) @f$.
   *
   * @param degree Degree of the harmonic.
   * @param order Order of the harmonic.
   * @param angles Structure containing the theta and phi angles.
   * @return The complex value of the spherical harmonic.
   */
  static std::complex<real_t> sphericalHarmonic(int degree, int order, SphericalAngles angles);

  /**
   * @struct Wigner3jParams
   * @brief Bundles angular momentum indices for the Wigner 3-j symbol calculation.
   */
  struct Wigner3jParams {
    int j_one{0};   ///< First angular momentum quantum number.
    int j_two{0};   ///< Second angular momentum quantum number.
    int j_three{0}; ///< Third angular momentum quantum number.
    int m_one{0};   ///< First magnetic projection quantum number.
    int m_two{0};   ///< Second magnetic projection quantum number.
    int m_three{0}; ///< Third magnetic projection quantum number.
  };

  /**
   * @brief Computes the Wigner 3-j symbol using the Racah formula.
   *
   * @param params Wigner 3-j symbol quantum numbers.
   * @return The scalar value of the 3-j symbol.
   */
  static real_t wigner3j(Wigner3jParams params);
};

} // namespace correlation::calculators
