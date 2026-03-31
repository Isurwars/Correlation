// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "DynamicsAnalyzer.hpp"
#include "math/Constants.hpp"
#include "math/FFTUtils.hpp"
#include "math/LinearAlgebra.hpp"
#include "physics/PhysicalData.hpp"
#include <algorithm>
#include <numeric>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

//---------------------------------------------------------------------------//
//--------------------------- Calculation Methods ---------------------------//
//---------------------------------------------------------------------------//

std::vector<double> DynamicsAnalyzer::calculateVACF(const Trajectory &traj,
                                                    int max_correlation_frames,
                                                    size_t start_frame,
                                                    size_t end_frame) {
  const auto &velocities = traj.getVelocities();

  if (velocities.empty()) {
    return {};
  }

  size_t total_frames = velocities.size();
  size_t num_atoms = velocities[0].size();

  start_frame = std::min(start_frame, total_frames > 0 ? total_frames - 1 : 0);
  end_frame = std::min(end_frame, total_frames);
  if (start_frame >= end_frame) {
    return {};
  }

  size_t num_frames = end_frame - start_frame;

  if (max_correlation_frames < 0 ||
      static_cast<size_t>(max_correlation_frames) >= num_frames) {
    max_correlation_frames = static_cast<int>(num_frames) - 1;
  }

  // Get masses for COM removal
  const auto &frames = traj.getFrames();
  if (frames.empty())
    return {};
  // Use the structure at start_frame for consistent masses
  const auto &atoms = frames[start_frame].atoms();
  std::vector<double> masses(num_atoms);
  double total_mass = 0.0;
  for (size_t i = 0; i < num_atoms; ++i) {
    try {
      masses[i] = correlation::physics::getAtomicMass(
          atoms[i].element().symbol);
    } catch (const std::out_of_range &) {
      masses[i] = 1.0;
    }
    total_mass += masses[i];
  }

  // Create corrected velocities (COM removed) mapped relative to start_frame
  std::vector<std::vector<correlation::math::Vector3<double>>>
      atom_velocities(
          num_atoms,
          std::vector<correlation::math::Vector3<double>>(num_frames));

  for (size_t t = 0; t < num_frames; ++t) {
    const size_t traj_t = start_frame + t;
    correlation::math::Vector3<double> momentum_sum = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < num_atoms; ++i) {
      momentum_sum += velocities[traj_t][i] * masses[i];
    }
    correlation::math::Vector3<double> v_com =
        momentum_sum / total_mass;

    for (size_t i = 0; i < num_atoms; ++i) {
      atom_velocities[i][t] = velocities[traj_t][i] - v_com;
    }
  }

  std::vector<double> vacf(max_correlation_frames + 1, 0.0);
  std::vector<int> counts(max_correlation_frames + 1, 0);

  // Count time origins per lag analytically
  for (int lag = 0; lag <= max_correlation_frames; ++lag) {
    counts[lag] = num_frames - lag;
  }

  // Thread-local VACF accumulators: eliminates all mutex acquisitions.
  tbb::enumerable_thread_specific<std::vector<double>> ets_vacf(
      [&] { return std::vector<double>(max_correlation_frames + 1, 0.0); });

  // Thread-local FFT workspace: one heap allocation per thread, reused
  // per-atom across all three autocorrelate calls.
  tbb::enumerable_thread_specific<std::vector<std::complex<double>>> ets_ws;

  tbb::parallel_for(size_t(0), num_atoms, [&](size_t i) {
    auto &local_vacf = ets_vacf.local();
    auto &ws = ets_ws.local();
    const auto &v_i = atom_velocities[i];

    std::vector<double> vx(num_frames), vy(num_frames), vz(num_frames);
    for (size_t t = 0; t < num_frames; ++t) {
      vx[t] = v_i[t].x();
      vy[t] = v_i[t].y();
      vz[t] = v_i[t].z();
    }

    auto S2_x = correlation::math::autocorrelate(vx, ws);
    auto S2_y = correlation::math::autocorrelate(vy, ws);
    auto S2_z = correlation::math::autocorrelate(vz, ws);

    for (int lag = 0; lag <= max_correlation_frames; ++lag) {
      local_vacf[lag] += S2_x[lag] + S2_y[lag] + S2_z[lag];
    }
  });

  // Serial merge of per-thread accumulators into the shared vacf[]
  for (const auto &local_vacf : ets_vacf) {
    for (int lag = 0; lag <= max_correlation_frames; ++lag) {
      vacf[lag] += local_vacf[lag];
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

std::vector<double> DynamicsAnalyzer::calculateMSD(const Trajectory &traj,
                                                   int max_correlation_frames,
                                                   size_t start_frame,
                                                   size_t end_frame) {
  const auto &frames = traj.getFrames();
  if (frames.empty()) {
    return {};
  }

  size_t total_frames = frames.size();
  size_t num_atoms = frames[0].atoms().size();

  start_frame = std::min(start_frame, total_frames > 0 ? total_frames - 1 : 0);
  end_frame = std::min(end_frame, total_frames);
  if (start_frame >= end_frame) {
    return {};
  }

  size_t num_frames = end_frame - start_frame;

  // Default: use half the trajectory length to ensure good statistics
  if (max_correlation_frames < 0 ||
      static_cast<size_t>(max_correlation_frames) >= num_frames) {
    max_correlation_frames = static_cast<int>(num_frames / 2);
  }
  if (max_correlation_frames == 0) {
    return {};
  }

  // --- Build unwrapped trajectories using minimum image convention ---
  // For each atom i and frame t (relative to start_frame):
  //   unwrapped[i][t] = sum_{s=0}^{t-1} min_image( r(s+1) - r(s) )
  // This correctly handles PBC crossings without needing explicit unwrapping.

  std::vector<std::vector<correlation::math::Vector3<double>>>
      unwrapped(num_atoms,
                std::vector<correlation::math::Vector3<double>>(
                    num_frames, {0.0, 0.0, 0.0}));

  for (size_t t = 1; t < num_frames; ++t) {
    const size_t tf = start_frame + t;
    const size_t tf_prev = start_frame + t - 1;

    // Get box vectors for minimum image (use current frame)
    const auto &lattice = frames[tf].latticeVectors();
    correlation::math::Vector3<double> box = {
        lattice[0][0], lattice[1][1], lattice[2][2]};
    bool use_pbc = (box[0] > 0.0 && box[1] > 0.0 && box[2] > 0.0);

    const auto &curr_atoms = frames[tf].atoms();
    const auto &prev_atoms = frames[tf_prev].atoms();

    for (size_t i = 0; i < num_atoms; ++i) {
      correlation::math::Vector3<double> dr =
          curr_atoms[i].position() - prev_atoms[i].position();

      // Apply minimum image convention
      if (use_pbc) {
        for (int d = 0; d < 3; ++d) {
          if (dr[d] > box[d] * 0.5)
            dr[d] -= box[d];
          if (dr[d] < -box[d] * 0.5)
            dr[d] += box[d];
        }
      }

      // Accumulate unwrapped position: unwrapped[i][t] = unwrapped[i][t-1] + dr
      unwrapped[i][t] = unwrapped[i][t - 1] + dr;
    }
  }

  // --- Compute MSD using Wiener-Khinchin Fast Correlation Algorithm ---
  // MSD(k) = (S1(k) - 2 * S2(k)) / (N - k)
  // S1(k) = sum_{i=0}^{N-1-k} (r_{i+k}^2 + r_i^2)
  // S2(k) = sum_{i=0}^{N-1-k} r_{i+k} * r_i

  std::vector<double> counts(max_correlation_frames + 1, 0.0);
  for (int lag = 0; lag <= max_correlation_frames; ++lag) {
    counts[lag] = static_cast<double>(num_frames - lag);
  }

  // Thread-local MSD accumulators to eliminate mutex contention
  tbb::enumerable_thread_specific<std::vector<double>> ets_msd(
      [&] { return std::vector<double>(max_correlation_frames + 1, 0.0); });

  // Thread-local FFT workspace reused across per-atom autocorrelate calls.
  tbb::enumerable_thread_specific<std::vector<std::complex<double>>> ets_ws_msd;

  tbb::parallel_for(size_t(0), num_atoms, [&](size_t i) {
    auto &local_msd = ets_msd.local();
    auto &ws = ets_ws_msd.local();
    const auto &u_i = unwrapped[i];

    std::vector<double> x(num_frames), y(num_frames), z(num_frames);
    std::vector<double> r_sq(num_frames);
    for (size_t t = 0; t < num_frames; ++t) {
      x[t] = u_i[t].x();
      y[t] = u_i[t].y();
      z[t] = u_i[t].z();
      r_sq[t] = x[t] * x[t] + y[t] * y[t] + z[t] * z[t];
    }

    auto S2_x = correlation::math::autocorrelate(x, ws);
    auto S2_y = correlation::math::autocorrelate(y, ws);
    auto S2_z = correlation::math::autocorrelate(z, ws);

    std::vector<double> S1(max_correlation_frames + 1, 0.0);
    S1[0] = 2.0 * std::accumulate(r_sq.begin(), r_sq.end(), 0.0);
    for (int k = 1; k <= max_correlation_frames; ++k) {
      S1[k] = S1[k - 1] - r_sq[k - 1] - r_sq[num_frames - k];
    }

    for (int lag = 1; lag <= max_correlation_frames; ++lag) {
      double S2 = S2_x[lag] + S2_y[lag] + S2_z[lag];
      local_msd[lag] += (S1[lag] - 2.0 * S2);
    }
  });

  std::vector<double> msd(max_correlation_frames + 1, 0.0);
  // Serial merge of per-thread accumulators
  for (const auto &local_msd : ets_msd) {
    for (int lag = 1; lag <= max_correlation_frames; ++lag) {
      msd[lag] += local_msd[lag];
    }
  }

  // Normalize by number of time origins and number of atoms
  msd[0] = 0.0;
  for (int lag = 1; lag <= max_correlation_frames; ++lag) {
    if (counts[lag] > 0) {
      msd[lag] /= (counts[lag] * static_cast<double>(num_atoms));
    }
  }

  return msd;
}

std::vector<double> DynamicsAnalyzer::calculateNormalizedVACF(
    const Trajectory &traj, int max_correlation_frames, size_t start_frame,
    size_t end_frame) {
  std::vector<double> vacf =
      calculateVACF(traj, max_correlation_frames, start_frame, end_frame);
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
  // Define standard deviation for the Gaussian window (e.g., stopping at 3
  // sigma at t_max)
  double sigma = t_max >= 1e-6 ? t_max / 3.0 : 1.0;
  for (size_t i = 0; i < num_frames; ++i) {
    double t = i * dt;
    double window = std::exp(-0.5 * std::pow(t / sigma, 2.0));
    windowed_vacf[i] = vacf[i] * window;
  }

  std::vector<size_t> freq_indices(num_freq_points);
  std::iota(freq_indices.begin(), freq_indices.end(), 0);

  tbb::parallel_for(size_t(0), num_freq_points, [&](size_t k) {
    double nu = k * d_nu; // Frequency in THz
    frequencies[k] = nu;

    double integral_real = 0.0;
    double integral_imag = 0.0;

    // $\omega = 2 * \pi * \nu * 0.001$ (to handle THz to fs
    // scale)
    double theta = correlation::math::two_pi * nu * dt * 0.001;

    if (std::abs(theta) < 1e-6) {
      for (size_t i = 0; i < num_frames; ++i) {
        double val = windowed_vacf[i];
        double arg = theta * i;

        double term_real = val * std::cos(arg);

        if (i == 0 || i == num_frames - 1) {
          integral_real += 0.5 * term_real;
        } else {
          integral_real += term_real;
        }
      }
      integral_real *= dt;
    } else {
      double theta2 = theta * theta;
      double theta3 = theta2 * theta;
      double sin_theta = std::sin(theta);
      double cos_theta = std::cos(theta);

      double alpha = 1.0 / theta + sin_theta * cos_theta / theta2 -
                     2.0 * sin_theta * sin_theta / theta3;
      double beta = 2.0 * ((1.0 + cos_theta * cos_theta) / theta2 -
                           2.0 * sin_theta * cos_theta / theta3);
      double gamma = 4.0 * (sin_theta / theta3 - cos_theta / theta2);

      size_t two_N = (num_frames % 2 == 0) ? num_frames - 2 : num_frames - 1;
      if (two_N < 0)
        two_N = 0;

      double f_0 = windowed_vacf[0];
      double f_2N = windowed_vacf[two_N];

      double arg_2N = theta * two_N;

      double sum_odd_cos = 0.0, sum_odd_sin = 0.0;
      for (size_t i = 1; i < two_N; i += 2) {
        double arg = theta * i;
        sum_odd_cos += windowed_vacf[i] * std::cos(arg);
        sum_odd_sin += windowed_vacf[i] * std::sin(arg);
      }

      double sum_even_cos = 0.0, sum_even_sin = 0.0;
      for (size_t i = 2; i < two_N; i += 2) {
        double arg = theta * i;
        sum_even_cos += windowed_vacf[i] * std::cos(arg);
        sum_even_sin += windowed_vacf[i] * std::sin(arg);
      }

      double even_cos_trapz =
          sum_even_cos + 0.5 * (f_0 * std::cos(0.0) + f_2N * std::cos(arg_2N));
      double even_sin_trapz =
          sum_even_sin + 0.5 * (f_0 * std::sin(0.0) + f_2N * std::sin(arg_2N));

      double bound_cos = f_2N * std::sin(arg_2N);
      double bound_sin = f_0 * std::cos(0.0) - f_2N * std::cos(arg_2N);

      integral_real = dt * (alpha * bound_cos + beta * even_cos_trapz +
                            gamma * sum_odd_cos);
      integral_imag = dt * (alpha * bound_sin + beta * even_sin_trapz +
                            gamma * sum_odd_sin);

      if (num_frames % 2 == 0 && num_frames > 1) {
        double val1 = windowed_vacf[num_frames - 2];
        double val2 = windowed_vacf[num_frames - 1];

        double arg1 = theta * (num_frames - 2);
        double arg2 = theta * (num_frames - 1);

        integral_real +=
            0.5 * dt * (val1 * std::cos(arg1) + val2 * std::cos(arg2));
        integral_imag +=
            0.5 * dt * (val1 * std::sin(arg1) + val2 * std::sin(arg2));
      }
    }

    double damping = std::exp(-0.5 * std::pow(nu / (0.5 * nyquist_thz), 2.0));
    intensities_real[k] = integral_real * damping;
    intensities_imag[k] = integral_imag * damping;
  });

  return {frequencies, intensities_real, intensities_imag};
}
