// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "DynamicsAnalyzer.hpp"
#include <iostream>

std::vector<double> DynamicsAnalyzer::calculateVACF(const Trajectory &traj, int max_correlation_frames) {
  const auto &velocities = traj.getVelocities();
  
  if (velocities.empty()) {
      return {};
  }

  size_t num_frames = velocities.size();
  size_t num_atoms = velocities[0].size();
  
  // Clamp max correlation frames
  if (max_correlation_frames < 0 || static_cast<size_t>(max_correlation_frames) >= num_frames) {
      max_correlation_frames = static_cast<int>(num_frames) - 1;
  }

  std::vector<double> vacf(max_correlation_frames + 1, 0.0);
  std::vector<int> counts(max_correlation_frames + 1, 0);

  // Iterate over time origins (t0)
  // To get good statistics, we can average over all possible start times for each lag
  for (size_t t0 = 0; t0 < num_frames; ++t0) {
      for (int lag = 0; lag <= max_correlation_frames; ++lag) {
          size_t t_plus_lag = t0 + lag;
          
          if (t_plus_lag < num_frames) {
              double dot_product_sum = 0.0;
              for (size_t i = 0; i < num_atoms; ++i) {
                  dot_product_sum += linalg::dot(velocities[t0][i], velocities[t_plus_lag][i]);
              }
              vacf[lag] += dot_product_sum;
              counts[lag]++;
          }
      }
  }

  // Normalize by number of time origins and number of atoms
  for (int lag = 0; lag <= max_correlation_frames; ++lag) {
      if (counts[lag] > 0) {
          vacf[lag] /= (counts[lag] * static_cast<double>(num_atoms));
      }
  }

  return vacf;
}

std::vector<double> DynamicsAnalyzer::calculateNormalizedVACF(const Trajectory &traj, int max_correlation_frames) {
    std::vector<double> vacf = calculateVACF(traj, max_correlation_frames);
    if (!vacf.empty() && vacf[0] != 0.0) {
        double normalization_factor = 1.0 / vacf[0];
        for (double &val : vacf) {
            val *= normalization_factor;
        }
    }
    return vacf;
}
