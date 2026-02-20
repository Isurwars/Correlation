// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "DynamicsAnalyzer.hpp"
#include "PhysicalData.hpp"
#include <algorithm>
#include <execution>
#include <mutex>
#include <numeric>

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

  if (max_correlation_frames < 0 ||
      static_cast<size_t>(max_correlation_frames) >= num_frames) {
    max_correlation_frames = static_cast<int>(num_frames) - 1;
  }

  // Get masses for COM removal
  const auto &frames = traj.getFrames();
  if (frames.empty())
    return {};
  const auto &atoms = frames[0].atoms();
  std::vector<double> masses(num_atoms);
  double total_mass = 0.0;
  for (size_t i = 0; i < num_atoms; ++i) {
    try {
      masses[i] = AtomicMasses::get(atoms[i].element().symbol);
    } catch (const std::out_of_range &) {
      // Fallback if mass not found (shouldn't happen for standard elements)
      masses[i] = 1.0;
    }
    total_mass += masses[i];
  }

  // Create corrected velocities (COM removed), transposed to [atom][frame]
  std::vector<std::vector<linalg::Vector3<double>>> atom_velocities(
      num_atoms, std::vector<linalg::Vector3<double>>(num_frames));

  for (size_t t = 0; t < num_frames; ++t) {
    linalg::Vector3<double> momentum_sum = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < num_atoms; ++i) {
      momentum_sum += velocities[t][i] * masses[i];
    }
    linalg::Vector3<double> v_com = momentum_sum / total_mass;

    for (size_t i = 0; i < num_atoms; ++i) {
      atom_velocities[i][t] = velocities[t][i] - v_com;
    }
  }

  std::vector<double> vacf(max_correlation_frames + 1, 0.0);
  std::vector<int> counts(max_correlation_frames + 1, 0);

  // Count time origins per lag analytically
  for (int lag = 0; lag <= max_correlation_frames; ++lag) {
    counts[lag] = num_frames - lag;
  }

  // Parallel computation over atoms
  std::vector<size_t> atom_indices(num_atoms);
  std::iota(atom_indices.begin(), atom_indices.end(), 0);

  std::mutex vacf_mutex;

  std::for_each(
      std::execution::par, atom_indices.begin(), atom_indices.end(),
      [&](size_t i) {
        std::vector<double> local_vacf(max_correlation_frames + 1, 0.0);
        const auto &v_i = atom_velocities[i];

        for (size_t t0 = 0; t0 < num_frames; ++t0) {
          int max_lag = std::min(max_correlation_frames,
                                 static_cast<int>(num_frames - 1 - t0));
          for (int lag = 0; lag <= max_lag; ++lag) {
            local_vacf[lag] += linalg::dot(v_i[t0], v_i[t0 + lag]);
          }
        }

        std::lock_guard<std::mutex> lock(vacf_mutex);
        for (int lag = 0; lag <= max_correlation_frames; ++lag) {
          vacf[lag] += local_vacf[lag];
        }
      });

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

  // Pre-calculate windowed VACF
  std::vector<double> windowed_vacf(num_frames);
  for (size_t i = 0; i < num_frames; ++i) {
    double t = i * dt;
    double window = 0.5 * (1.0 + std::cos(constants::pi * t / t_max));
    windowed_vacf[i] = vacf[i] * window;
  }

  std::vector<size_t> freq_indices(num_freq_points);
  std::iota(freq_indices.begin(), freq_indices.end(), 0);

  std::for_each(
      std::execution::par, freq_indices.begin(), freq_indices.end(),
      [&](size_t k) {
        double nu = k * d_nu; // Frequency in THz
        frequencies[k] = nu;

        double integral_real = 0.0;
        double integral_imag = 0.0;

        // Trapezoidal integration
        for (size_t i = 0; i < num_frames; ++i) {
          double t = i * dt; // Time in fs
          double val = windowed_vacf[i];

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
        intensities_real[k] = integral_real * dt;
        intensities_imag[k] = integral_imag * dt;
      });

  return {frequencies, intensities_real, intensities_imag};
}
