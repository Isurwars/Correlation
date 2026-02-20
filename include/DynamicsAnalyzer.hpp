// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "Trajectory.hpp"
#include <tuple>
#include <vector>

/**
 * @brief Static class providing dynamics analysis methods.
 *
 * This class contains stateless static methods for calculating time-dependent
 * properties like VACF and VDOS from Trajectory data.
 */
class DynamicsAnalyzer {
public:
  //-------------------------------------------------------------------------//
  //--------------------------- Calculation Methods -------------------------//
  //-------------------------------------------------------------------------//
  /**
   * @brief Calculates the Velocity Autocorrelation Function (VACF).
   *
   * C(t) = < v(t0) * v(t0 + t) >
   *
   * @param traj The trajectory containing pre-calculated velocities.
   * @param max_correlation_frames The maximum time lag (in frames) to calculate
   * correlation for.
   * @return A vector containing the VACF values for lag times 0 to
   * max_correlation_frames.
   */
  static std::vector<double>
  calculateVACF(const Trajectory &traj, int max_correlation_frames,
                size_t start_frame = 0,
                size_t end_frame = static_cast<size_t>(-1));

  /**
   * @brief Calculates the Normalized Velocity Autocorrelation Function.
   * Divides the VACF by its initial value (at t=0).
   *
   * @param traj The trajectory containing pre-calculated velocities.
   * @param max_correlation_frames The maximum time lag.
   * @return A vector containing the normalized VACF values.
   */
  static std::vector<double>
  calculateNormalizedVACF(const Trajectory &traj, int max_correlation_frames,
                          size_t start_frame = 0,
                          size_t end_frame = static_cast<size_t>(-1));

  /**
   * @brief Calculates the Vibrational Density of States (VDOS) from the VACF.
   *
   * Real part (Cosine Transform): \int_0^\infty VACF(t) * cos(omega * t) dt
   * Imaginary part (Sine Transform): \int_0^\infty VACF(t) * sin(omega * t) dt
   *
   * @param vacf The Velocity Autocorrelation Function.
   * @param dt The time step between frames (in femtoseconds).
   * @param params Optional parameters for windowing/padding.
   * @return A tuple: {frequencies (THz), real_intensities, imag_intensities}.
   */
  static std::tuple<std::vector<double>, std::vector<double>,
                    std::vector<double>>
  calculateVDOS(const std::vector<double> &vacf, double dt);
};
