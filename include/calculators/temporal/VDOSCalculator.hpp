/**
 * @file VDOSCalculator.hpp
 * @brief Vibrational density of states (VDOS) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "calculators/BaseCalculator.hpp"

namespace correlation::core {
class Cell;
}

namespace correlation::calculators {

/**
 * @brief Parameters for the Vibrational Density of States (VDOS) calculation.
 *
 * Bundles frequency bound parameters into a single struct to
 * prevent accidental argument swapping (bugprone-easily-swappable-parameters).
 */
struct VDOSParams {
  real_t max_imag_freq =
      static_cast<real_t>(5.0); ///< Maximum negative (imaginary) frequency bound in THz.
  real_t max_real_freq =
      static_cast<real_t>(10.0); ///< Maximum positive (real) frequency bound in THz.
};

/**
 * @class VDOSCalculator
 * @brief Computes the Vibrational Density of States (VDOS) from the VACF.
 */
class VDOSCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string_view getName() const override {
    return "Vibrational Density of States";
  }
  [[nodiscard]] std::string_view getShortName() const override { return "VDOS"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Dynamical"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Computes the Vibrational Density of States (VDOS).";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return false; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return true; }

  void calculateTrajectory(correlation::analysis::DistributionFunctions &dists,
                           const correlation::core::Trajectory &traj,
                           const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Computes the Vibrational Density of States (VDOS).
   * @param[in] vacf_hist Input Velocity Autocorrelation Function (VACF) histogram.
   * @param[in] params Frequency bounds parameters for VDOS calculation.
   * @return A histogram representing intensity vs frequency (THz/cm^-1).
   * @throws std::logic_error if VACF data has fewer than 2 data points.
   */
  [[nodiscard]] static correlation::analysis::Histogram
  calculate(const correlation::analysis::Histogram &vacf_hist, const VDOSParams &params = {});

  /**
   * @brief Computes the Vibrational Density of States (VDOS) with explicit frequency bounds.
   * @param[in] vacf_hist Input Velocity Autocorrelation Function (VACF) histogram.
   * @param[in] max_imag_freq Maximum negative (imaginary) frequency bound in THz.
   * @param[in] max_real_freq Maximum positive (real) frequency bound in THz.
   * @return A histogram representing intensity vs frequency (THz/cm^-1).
   * @throws std::logic_error if VACF data has fewer than 2 data points.
   */
  [[nodiscard]] static correlation::analysis::Histogram
  calculate(const correlation::analysis::Histogram &vacf_hist, real_t max_imag_freq,
            real_t max_real_freq);

  /**
   * @brief Computes the recommended maximum positive vDoS frequency based on the lightest element.
   *
   * Continuous mass scaling: ν_max = max(10.0, 50.0 / sqrt(m_min)), bounded by Nyquist frequency.
   *
   * @param[in] cell The atomic cell containing chemical elements.
   * @param[in] time_step The trajectory time step in femtoseconds.
   * @return Recommended maximum real frequency in THz.
   */
  [[nodiscard]] static real_t computeRecommendedMaxFrequency(const correlation::core::Cell &cell,
                                                             real_t time_step) noexcept;
};

} // namespace correlation::calculators
