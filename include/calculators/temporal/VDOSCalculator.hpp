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
 */
struct VDOSParams {
  real_t max_real_freq = static_cast<real_t>(10.0); ///< Maximum positive frequency bound in THz.
  correlation::analysis::VDOSWindowType window_type =
      correlation::analysis::VDOSWindowType::Hann; ///< Time-domain tapering half-window.
  size_t zero_pad_factor = 4; ///< Zero-padding multiplier for FFT interpolation (>= 1).
  bool normalize = true;      ///< Area-normalize VDOS over [0, E_Nyquist] to 1 (units: meV⁻¹).
};

/**
 * @brief Result of 2PT spectral decomposition into vibrational and translational modes.
 */
struct VDOSDecomposition {
  correlation::analysis::Histogram
      vibrational; ///< Vibrational-only VDOS g_vib(ω), g_vib(0) = 0, ∫g_vib d(ħω) = 1.
  correlation::analysis::Histogram
      translational;      ///< Translational/diffusive DoS g_trans(ω), ∫g_trans d(ħω) = 1.
  real_t fluidicity{0.0}; ///< 2PT fluidicity parameter f ∈ [0, 1] (0 = solid, 1 = gas).
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
   * @brief Computes 2PT-decomposed vibrational and translational DoS histograms on ω >= 0.
   * @param[in] vacf_hist Input Velocity Autocorrelation Function (VACF) histogram.
   * @param[in] params Frequency bound, windowing, zero-padding, and normalization parameters.
   * @param[in] cell Optional reference cell for mass-weighted partials and 2PT fluidicity.
   * @return VDOSDecomposition containing independently normalized vibrational and translational
   *         histograms plus the 2PT fluidicity f ∈ [0, 1].
   * @throws std::logic_error if VACF data has fewer than 2 data points.
   */
  [[nodiscard]] static VDOSDecomposition
  calculateDecomposed(const correlation::analysis::Histogram &vacf_hist,
                      const VDOSParams &params = {}, const correlation::core::Cell *cell = nullptr);

  /**
   * @brief Computes the vibrational-only VDOS g_vib(ω) on real energies (meV, >= 0).
   * @param[in] vacf_hist Input Velocity Autocorrelation Function (VACF) histogram.
   * @param[in] params Frequency bound, windowing, zero-padding, and normalization parameters.
   * @param[in] cell Optional reference cell for mass-weighted per-element partial VDOS.
   * @return A histogram representing g_vib(ω) vs phonon energy (meV).
   * @throws std::logic_error if VACF data has fewer than 2 data points.
   */
  [[nodiscard]] static correlation::analysis::Histogram
  calculate(const correlation::analysis::Histogram &vacf_hist, const VDOSParams &params = {},
            const correlation::core::Cell *cell = nullptr);

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
