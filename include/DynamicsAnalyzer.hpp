// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "Trajectory.hpp"
#include <vector>

class DynamicsAnalyzer {
public:
  /**
   * @brief Calculates the Velocity Autocorrelation Function (VACF).
   * 
   * C(t) = < v(t0) * v(t0 + t) >
   * 
   * @param traj The trajectory containing pre-calculated velocities.
   * @param max_correlation_frames The maximum time lag (in frames) to calculate correlation for.
   * @return A vector containing the VACF values for lag times 0 to max_correlation_frames.
   */
  static std::vector<double> calculateVACF(const Trajectory &traj, int max_correlation_frames);

  /**
   * @brief Calculates the Normalized Velocity Autocorrelation Function.
   * Divides the VACF by its initial value (at t=0).
   * 
   * @param traj The trajectory containing pre-calculated velocities.
   * @param max_correlation_frames The maximum time lag.
   * @return A vector containing the normalized VACF values.
   */
  static std::vector<double> calculateNormalizedVACF(const Trajectory &traj, int max_correlation_frames);
};
