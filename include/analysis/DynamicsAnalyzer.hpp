/**
 * @file DynamicsAnalyzer.hpp
 * @brief Static methods for dynamics properties (VACF, VDOS, MSD).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
#include "core/Trajectory.hpp"

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace correlation::analysis {

/**
 * @struct PartialVACFResult
 * @brief Container for total (unweighted and mass-weighted) and per-element VACF curves.
 */
struct PartialVACFResult {
  std::vector<real_t> total_unweighted;    ///< Unweighted total VACF (Å² fs⁻²), for Green-Kubo D.
  std::vector<real_t> total_mass_weighted; ///< Mass-weighted total VACF (Å² fs⁻²), for phonon VDOS.
  std::map<std::string, std::vector<real_t>> element_vacfs; ///< Per-element unweighted VACFs.
  std::map<std::string, real_t> element_mass_weights;       ///< Normalized mass weights w_α.
};

/**
 * @brief Static class providing dynamics analysis methods.
 *
 * This class contains stateless static methods for calculating time-dependent
 * properties like VACF and VDOS from correlation::core::Trajectory data.
 */
class DynamicsAnalyzer {
public:
  DynamicsAnalyzer() = delete;

  /** @name Calculation Methods */
  ///@{
  /**
   * @brief Calculates the Velocity Autocorrelation Function (VACF).
   *
   * @f$ C(t) = \langle \mathbf{v}(t_0) \cdot \mathbf{v}(t_0 + t) \rangle @f$
   *
   * @param traj The trajectory containing pre-calculated velocities.
   * @param max_correlation_frames The maximum time lag (in frames) to calculate
   * correlation for.
   * @param start_frame First frame to include (default: 0).
   * @param end_frame One-past-last frame to include (default: all frames).
   * @return A vector containing the VACF values for lag times 0 to
   * max_correlation_frames.
   */
  [[nodiscard]] static std::vector<real_t>
  calculateVACF(const correlation::core::Trajectory &traj, MaxFrames max_correlation_frames,
                StartFrame start_frame = {0}, EndFrame end_frame = {static_cast<size_t>(-1)});

  /**
   * @brief Calculates the Normalized Velocity Autocorrelation Function.
   * Divides the VACF by its initial value (at t=0).
   *
   * @param traj The trajectory containing pre-calculated velocities.
   * @param max_correlation_frames The maximum time lag.
   * @param start_frame First frame to include (default: 0).
   * @param end_frame One-past-last frame to include (default: all frames).
   * @return A vector containing the normalized VACF values.
   */
  [[nodiscard]] static std::vector<real_t>
  calculateNormalizedVACF(const correlation::core::Trajectory &traj,
                          MaxFrames max_correlation_frames, StartFrame start_frame = {0},
                          EndFrame end_frame = {static_cast<size_t>(-1)});

  /**
   * @brief Calculates the Mean Squared Displacement (MSD).
   *
   * MSD(t) = < |r_i(t0 + t) - r_i(t0)|^2 >
   *
   * Uses unwrapped (continuous) trajectories to correctly handle atoms
   * crossing periodic boundary conditions. The displacement is accumulated
   * frame-by-frame using the minimum image convention, so no stored unwrapped
   * positions are needed. Time-averaging over all valid origins is used for
   * maximum statistical efficiency.
   *
   * @param traj The trajectory containing atomic positions.
   * @param max_correlation_frames The maximum lag (in frames) to compute.
   *        Pass -1 to use half the trajectory length (recommended).
   * @param start_frame First frame to include (default: 0).
   * @param end_frame One-past-last frame to include (default: all frames).
   * @return A vector of MSD values indexed by lag (in Å²).
   */
  [[nodiscard]] static std::vector<real_t>
  calculateMSD(const correlation::core::Trajectory &traj, MaxFrames max_correlation_frames,
               StartFrame start_frame = {0}, EndFrame end_frame = {static_cast<size_t>(-1)});

  /**
   * @brief Calculates total (unweighted and mass-weighted) and per-element VACFs.
   *
   * @param traj The trajectory containing pre-calculated velocities.
   * @param max_correlation_frames Maximum time lag in frames (-1 = half trajectory).
   * @param start_frame First frame to include (default: 0).
   * @param end_frame One-past-last frame to include (default: all frames).
   * @return PartialVACFResult containing total and per-element VACFs and mass weights.
   */
  [[nodiscard]] static PartialVACFResult
  calculatePartialVACF(const correlation::core::Trajectory &traj, MaxFrames max_correlation_frames,
                       StartFrame start_frame = {0},
                       EndFrame end_frame = {static_cast<size_t>(-1)});

  /**
   * @brief Calculates the Vibrational Density of States (VDOS) from the VACF.
   *
   * Evaluates the real one-sided Fourier cosine transform via a zero-padded
   * even-symmetric FFT (DCT-I equivalent) after applying a smooth time-domain
   * half-window:
   * @f[ g(\nu) = 2 \int_0^{t_{\max}} C(t)\,w(t)\,\cos(2\pi \nu t)\,\mathrm{d}t @f]
   *
   * @param vacf The Velocity Autocorrelation Function.
   * @param time_step The time step between frames (in femtoseconds).
   * @param window_type Time-domain tapering half-window (default: Half-Hann).
   * @param zero_pad_factor Zero-padding multiplier for frequency interpolation (>= 1, default: 4).
   * @return A pair: {frequencies (THz, >= 0), intensities (>= 0)}.
   */
  [[nodiscard]] static std::pair<std::vector<real_t>, std::vector<real_t>>
  calculateVDOS(const std::vector<real_t> &vacf, real_t time_step,
                VDOSWindowType window_type = VDOSWindowType::Hann, size_t zero_pad_factor = 4);

  /**
   * @brief Computes the self-diffusion coefficient D from MSD using linear regression on the second
   * half of the time series.
   * @param time The time values (bins).
   * @param msd The MSD values.
   * @return The diffusion coefficient in Å²/fs.
   */
  [[nodiscard]] static real_t computeDiffusionCoefficientMSD(const std::vector<real_t> &time,
                                                             const std::vector<real_t> &msd);

  /**
   * @brief Computes the self-diffusion coefficient D from VACF using Green-Kubo integration
   * (trapezoidal rule).
   * @param time The time values (bins).
   * @param vacf The VACF values.
   * @return The diffusion coefficient in Å²/fs.
   */
  [[nodiscard]] static real_t computeDiffusionCoefficientVACF(const std::vector<real_t> &time,
                                                              const std::vector<real_t> &vacf);

  /**
   * @brief Computes the relaxation time tau from the normalized VACF using trapezoidal integration.
   * @param time The time values (bins).
   * @param normalized_vacf The normalized VACF values.
   * @return The relaxation time in fs.
   */
  [[nodiscard]] static real_t computeRelaxationTime(const std::vector<real_t> &time,
                                                    const std::vector<real_t> &normalized_vacf);

  ///@}
};

} // namespace correlation::analysis
