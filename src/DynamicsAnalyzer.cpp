// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "DynamicsAnalyzer.hpp"
#include "PhysicalData.hpp"

//---------------------------------------------------------------------------//
//--------------------------- Calculation Methods ---------------------------//
//---------------------------------------------------------------------------//

std::vector<double>
DynamicsAnalyzer::calculateVACF(const Trajectory &traj,
                                int max_correlation_frames) {
  const auto &velocities = traj.getVelocities();

  if (velocities.empty()) {
    return {};
  }

  size_t num_frames = velocities.size();
  size_t num_atoms = velocities[0].size();

  // Clamp max correlation frames
  if (max_correlation_frames < 0 ||
      static_cast<size_t>(max_correlation_frames) >= num_frames) {
    max_correlation_frames = static_cast<int>(num_frames) - 1;
  }

  std::vector<double> vacf(max_correlation_frames + 1, 0.0);
  std::vector<int> counts(max_correlation_frames + 1, 0);

  // Iterate over time origins (t0)
  // To get good statistics, we can average over all possible start times for
  // each lag
  for (size_t t0 = 0; t0 < num_frames; ++t0) {
    for (int lag = 0; lag <= max_correlation_frames; ++lag) {
      size_t t_plus_lag = t0 + lag;

      if (t_plus_lag < num_frames) {
        double dot_product_sum = 0.0;
        for (size_t i = 0; i < num_atoms; ++i) {
          dot_product_sum +=
              linalg::dot(velocities[t0][i], velocities[t_plus_lag][i]);
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

std::vector<double>
DynamicsAnalyzer::calculateNormalizedVACF(const Trajectory &traj,
                                          int max_correlation_frames) {
  std::vector<double> vacf = calculateVACF(traj, max_correlation_frames);
  if (!vacf.empty() && vacf[0] != 0.0) {
    double normalization_factor = 1.0 / vacf[0];
    for (double &val : vacf) {
      val *= normalization_factor;
    }
  }
  return vacf;
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
DynamicsAnalyzer::calculateVDOS(const std::vector<double> &vacf, double dt) {
  if (vacf.empty()) {
    return {};
  }

  size_t num_frames = vacf.size();
  double t_max = (num_frames - 1) * dt; // Total time of correlation

  // Define frequency range
  // Max frequency (Nyquist) is 1 / (2 * dt)
  // But we are interested in physical vibrations, usually up to ~100 THz or
  // ~3000 cm^-1 dt is in fs. 1 fs = 10^-15 s. 1 THz = 10^12 Hz. Nyquist freq in
  // THz = 1 / (2 * dt * 10^-15) * 10^-12 = 1000 / (2 * dt)

  double nyquist_thz = 1000.0 / (2.0 * dt);

  // Number of frequency points. Let's make it high resolution.
  size_t num_freq_points = 2000;
  double d_nu = nyquist_thz / num_freq_points; // Frequency step in THz

  std::vector<double> frequencies(num_freq_points);
  std::vector<double> intensities_real(num_freq_points);
  std::vector<double> intensities_imag(num_freq_points);

  for (size_t k = 0; k < num_freq_points; ++k) {
    double nu = k * d_nu; // Frequency in THz
    frequencies[k] = nu;

    double integral_real = 0.0;
    double integral_imag = 0.0;

    // Trapezoidal integration
    for (size_t i = 0; i < num_frames; ++i) {
      double t = i * dt; // Time in fs
      double val = vacf[i];

      // Apply Hann Window to reduce spectral leakage
      // W(t) = 0.5 * (1 + cos(pi * t / t_max))

      double window = 0.5 * (1.0 + std::cos(constants::pi * t / t_max));
      val *= window;

      // Cosine transform term: cos(2 * pi * nu * t)
      // nu is in THz (10^12/s), t in fs (10^-15 s).
      // nu * t = (10^12) * (10^-15) = 10^-3
      // So arg = 2 * pi * nu * t * 0.001
      double arg = 2.0 * constants::pi * nu * t * 0.001;

      double term_real = val * std::cos(arg);
      double term_imag = val * std::sin(arg);

      if (i == 0 || i == num_frames - 1) {
        integral_real += 0.5 * term_real;
        integral_imag += 0.5 * term_imag;
      } else {
        integral_real += term_real;
        integral_imag += term_imag;
      }
    }
    // Multiply by frequency to sharpen high-frequency features (and suppress
    // DC)
    intensities_real[k] = integral_real * dt;
    intensities_imag[k] = integral_imag * dt;
  }

  return {frequencies, intensities_real, intensities_imag};
}
