/**
 * @file SteinhardtCalculator.hpp
 * @brief Steinhardt bond-order parameter calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "DistributionFunctions.hpp"
#include <complex>

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

  void calculateFrame(DistributionFunctions &df,
                      const AnalysisSettings &settings) const override;

  // Static function to calculate Steinhardt histograms
  static std::map<std::string, Histogram>
  calculate(const Cell &cell, const StructureAnalyzer *neighbors);

  // Helper math functions
  static std::complex<double> sphericalHarmonic(int l, int m, double theta,
                                                double phi);
  static double wigner3j(int j1, int j2, int j3, int m1, int m2, int m3);
};
