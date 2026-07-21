/**
 * @file VDOSCalculator.hpp
 * @brief Vibrational density of states (VDOS) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

namespace correlation::calculators {

/**
 * @brief Parameters for the Vibrational Density of States (VDOS) calculation.
 *
 * Bundles frequency bound parameters into a single struct to
 * prevent accidental argument swapping (bugprone-easily-swappable-parameters).
 */
struct VDOSParams {
  real_t max_imag_freq = static_cast<real_t>(5.0);  ///< Maximum negative (imaginary) frequency bound in THz.
  real_t max_real_freq = static_cast<real_t>(10.0); ///< Maximum positive (real) frequency bound in THz.
};

/**
 * @class VDOSCalculator
 * @brief Computes the Vibrational Density of States (VDOS) from the VACF.
 */
class VDOSCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "vDoS"; }
  std::string getShortName() const override { return "vDoS"; }
  std::string getGroup() const override { return "Dynamic"; }
  std::string getDescription() const override { return "Computes the Vibrational Density of States (vDoS)."; }

  bool isFrameCalculator() const override { return false; }
  bool isTrajectoryCalculator() const override { return true; }

  void calculateTrajectory(correlation::analysis::DistributionFunctions &dists,
                           const correlation::core::Trajectory &traj,
                           const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Computes the Vibrational Density of States (VDOS).
   * @param vacf_hist Input Velocity Autocorrelation Function (VACF) histogram.
   * @param params Frequency bounds parameters for VDOS calculation.
   * @return A histogram representing intensity vs frequency (THz/cm^-1).
   */
  static correlation::analysis::Histogram calculate(const correlation::analysis::Histogram &vacf_hist,
                                                    const VDOSParams &params = {});

  /**
   * @brief Computes the Vibrational Density of States (VDOS) with explicit frequency bounds.
   * @param vacf_hist Input Velocity Autocorrelation Function (VACF) histogram.
   * @param max_imag_freq Maximum negative (imaginary) frequency bound in THz.
   * @param max_real_freq Maximum positive (real) frequency bound in THz.
   * @return A histogram representing intensity vs frequency (THz/cm^-1).
   */
  static correlation::analysis::Histogram calculate(const correlation::analysis::Histogram &vacf_hist,
                                                    real_t max_imag_freq, real_t max_real_freq);
};

} // namespace correlation::calculators
