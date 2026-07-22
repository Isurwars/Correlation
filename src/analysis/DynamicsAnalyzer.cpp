/**
 * @file DynamicsAnalyzer.cpp
 * @brief Implementation of VACF, VDOS, and MSD analysis.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#include "analysis/DynamicsAnalyzer.hpp"
#include "math/Constants.hpp"
#include "math/FFTUtils.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"
#include "physics/PhysicalData.hpp"

#include <algorithm>
#include <numeric>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <utility>

namespace correlation::analysis {

namespace {
std::pair<real_t, real_t> integrate_vdos_frequency(real_t theta, const std::vector<real_t> &windowed_vacf,
                                                   real_t time_step) {
  size_t const num_frames = windowed_vacf.size();
  real_t integral_real = 0.0;
  real_t integral_imag = 0.0;

  if (std::abs(theta) < 1e-6) {
    for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
      real_t const val = windowed_vacf[frame_idx];
      real_t const arg = theta * static_cast<real_t>(frame_idx);

      real_t const term_real = val * std::cos(arg);

      if (frame_idx == 0 || frame_idx == num_frames - 1) {
        integral_real += static_cast<real_t>(0.5) * term_real;
      } else {
        integral_real += term_real;
      }
    }
    integral_real *= time_step;
  } else {
    real_t const theta2 = theta * theta;
    real_t const theta3 = theta2 * theta;
    real_t const sin_theta = std::sin(theta);
    real_t const cos_theta = std::cos(theta);

    real_t const alpha = static_cast<real_t>(1.0) / theta + sin_theta * cos_theta / theta2 -
                         static_cast<real_t>(2.0) * sin_theta * sin_theta / theta3;
    real_t const beta = static_cast<real_t>(2.0) * ((static_cast<real_t>(1.0) + cos_theta * cos_theta) / theta2 -
                                                    static_cast<real_t>(2.0) * sin_theta * cos_theta / theta3);
    real_t const gamma = static_cast<real_t>(4.0) * (sin_theta / theta3 - cos_theta / theta2);

    size_t two_N = (num_frames % 2 == 0) ? num_frames - 2 : num_frames - 1;
    two_N = std::max<size_t>(two_N, 0);

    real_t const f_0 = windowed_vacf[0];
    real_t const f_2N = windowed_vacf[two_N];

    real_t const arg_2N = theta * static_cast<real_t>(two_N);

    real_t sum_odd_cos = 0.0;
    real_t sum_odd_sin = 0.0;
    for (size_t frame_idx = 1; frame_idx < two_N; frame_idx += 2) {
      real_t const arg = theta * static_cast<real_t>(frame_idx);
      sum_odd_cos += windowed_vacf[frame_idx] * std::cos(arg);
      sum_odd_sin += windowed_vacf[frame_idx] * std::sin(arg);
    }

    real_t sum_even_cos = 0.0;
    real_t sum_even_sin = 0.0;
    for (size_t frame_idx = 2; frame_idx < two_N; frame_idx += 2) {
      real_t const arg = theta * static_cast<real_t>(frame_idx);
      sum_even_cos += windowed_vacf[frame_idx] * std::cos(arg);
      sum_even_sin += windowed_vacf[frame_idx] * std::sin(arg);
    }

    real_t const even_cos_trapz =
        static_cast<real_t>(sum_even_cos + static_cast<real_t>(0.5) * (f_0 * std::cos(0.0) + f_2N * std::cos(arg_2N)));
    real_t const even_sin_trapz =
        static_cast<real_t>(sum_even_sin + static_cast<real_t>(0.5) * (f_0 * std::sin(0.0) + f_2N * std::sin(arg_2N)));

    real_t const bound_cos = f_2N * std::sin(arg_2N);
    real_t const bound_sin = static_cast<real_t>(f_0 * std::cos(0.0) - f_2N * std::cos(arg_2N));

    integral_real = time_step * (alpha * bound_cos + beta * even_cos_trapz + gamma * sum_odd_cos);
    integral_imag = time_step * (alpha * bound_sin + beta * even_sin_trapz + gamma * sum_odd_sin);

    if (num_frames % 2 == 0 && num_frames > 1) {
      real_t const val1 = windowed_vacf[num_frames - 2];
      real_t const val2 = windowed_vacf[num_frames - 1];

      real_t const arg1 = theta * static_cast<real_t>(num_frames - 2);
      real_t const arg2 = theta * static_cast<real_t>(num_frames - 1);

      integral_real += static_cast<real_t>(0.5) * time_step * (val1 * std::cos(arg1) + val2 * std::cos(arg2));
      integral_imag += static_cast<real_t>(0.5) * time_step * (val1 * std::sin(arg1) + val2 * std::sin(arg2));
    }
  }

  return {integral_real, integral_imag};
}
} // namespace

std::vector<real_t> DynamicsAnalyzer::calculateVACF(const correlation::core::Trajectory &traj,
                                                    MaxFrames max_correlation_frames_arg, StartFrame start_frame_arg,
                                                    EndFrame end_frame_arg) {
  int max_correlation_frames = max_correlation_frames_arg.value;
  size_t start_frame = start_frame_arg.value;
  size_t end_frame = end_frame_arg.value;
  if (traj.getFrameCount() == 0) {
    return {};
  }
  const auto &frames = traj.getFrames(); // materialises lazy trajectory for computation

  size_t const total_frames = frames.size();
  size_t const num_atoms = frames[0].atoms().size();

  start_frame = std::min(start_frame, total_frames > 0 ? total_frames - 1 : 0);
  end_frame = std::min(end_frame, total_frames);
  if (start_frame >= end_frame) {
    return {};
  }

  size_t num_frames = end_frame - start_frame;

  if (max_correlation_frames < 0 || std::cmp_greater_equal(max_correlation_frames, num_frames)) {
    max_correlation_frames = static_cast<int>(num_frames) - 1;
  }

  // Get masses for COM removal
  const auto &atoms = frames[start_frame].atoms();
  std::vector<real_t> masses(num_atoms);
  real_t total_mass = 0.0;
  for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
    try {
      masses[atom_idx] = correlation::physics::getAtomicMass(atoms[atom_idx].element().symbol);
    } catch (const std::out_of_range &) {
      masses[atom_idx] = 1.0;
    }
    total_mass += masses[atom_idx];
  }

  // Create corrected velocities (COM removed) mapped relative to start_frame
  std::vector<std::vector<correlation::math::Vector3<real_t>>> atom_velocities(
      num_atoms, std::vector<correlation::math::Vector3<real_t>>(num_frames));

  for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
    const size_t traj_t = start_frame + frame_idx;
    correlation::math::Vector3<real_t> momentum_sum = {0.0, 0.0, 0.0};
    for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
      momentum_sum += frames[traj_t].atoms()[atom_idx].velocity() * masses[atom_idx];
    }
    correlation::math::Vector3<real_t> const v_com = momentum_sum / total_mass;

    for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
      atom_velocities[atom_idx][frame_idx] = frames[traj_t].atoms()[atom_idx].velocity() - v_com;
    }
  }

  std::vector<real_t> vacf(max_correlation_frames + 1, 0.0);
  std::vector<int> counts(max_correlation_frames + 1, 0);

  // Count time origins per lag analytically
  for (int lag = 0; lag <= max_correlation_frames; ++lag) {
    counts[lag] = static_cast<int>(num_frames - static_cast<size_t>(lag));
  }

  // Thread-local VACF accumulators: eliminates all mutex acquisitions.
  tbb::enumerable_thread_specific<std::vector<real_t>> ets_vacf(
      [&] { return std::vector<real_t>(max_correlation_frames + 1, 0.0); });

  // Thread-local FFT workspace: one heap allocation per thread, reused
  // per-atom across all three autocorrelate calls.
  tbb::enumerable_thread_specific<std::vector<std::complex<real_t>>> ets_ws;

  tbb::parallel_for(static_cast<size_t>(0), num_atoms, [&](size_t atom_idx) {
    auto &local_vacf = ets_vacf.local();
    auto &workspace = ets_ws.local();
    const auto &v_i = atom_velocities[atom_idx];

    std::vector<real_t> vel_x(num_frames);
    std::vector<real_t> vel_y(num_frames);
    std::vector<real_t> vel_z(num_frames);
    for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
      vel_x[frame_idx] = v_i[frame_idx].x();
      vel_y[frame_idx] = v_i[frame_idx].y();
      vel_z[frame_idx] = v_i[frame_idx].z();
    }

    auto S2_x = correlation::math::autocorrelate(vel_x, workspace);
    auto S2_y = correlation::math::autocorrelate(vel_y, workspace);
    auto S2_z = correlation::math::autocorrelate(vel_z, workspace);

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
      vacf[lag] /= (static_cast<real_t>(counts[lag]) * static_cast<real_t>(num_atoms));
    }
  }

  return vacf;
}

std::vector<real_t> DynamicsAnalyzer::calculateMSD(const correlation::core::Trajectory &traj,
                                                   MaxFrames max_correlation_frames_arg, StartFrame start_frame_arg,
                                                   EndFrame end_frame_arg) {
  int max_correlation_frames = max_correlation_frames_arg.value;
  size_t start_frame = start_frame_arg.value;
  size_t end_frame = end_frame_arg.value;
  if (traj.getFrameCount() == 0) {
    return {};
  }
  const auto &frames = traj.getFrames(); // materialises lazy trajectory for computation

  size_t const total_frames = frames.size();
  size_t const num_atoms = frames[0].atoms().size();

  start_frame = std::min(start_frame, total_frames > 0 ? total_frames - 1 : 0);
  end_frame = std::min(end_frame, total_frames);
  if (start_frame >= end_frame) {
    return {};
  }

  size_t num_frames = end_frame - start_frame;

  // Default: use half the trajectory length to ensure good statistics
  if (max_correlation_frames < 0 || std::cmp_greater_equal(max_correlation_frames, num_frames)) {
    max_correlation_frames = static_cast<int>(num_frames / 2);
  }
  if (max_correlation_frames == 0) {
    return {};
  }

  // --- Build unwrapped trajectories using minimum image convention ---
  // For each atom i and frame frame_idx (relative to start_frame):
  //   unwrapped[i][frame_idx] = sum_{s=0}^{frame_idx-1} min_image( r(s+1) - r(s) )
  // This correctly handles PBC crossings without needing explicit unwrapping.

  std::vector<std::vector<correlation::math::Vector3<real_t>>> unwrapped(
      num_atoms, std::vector<correlation::math::Vector3<real_t>>(num_frames, {0.0, 0.0, 0.0}));

  std::vector<correlation::KahanAccumulator<real_t>> kahan_x(num_atoms);
  std::vector<correlation::KahanAccumulator<real_t>> kahan_y(num_atoms);
  std::vector<correlation::KahanAccumulator<real_t>> kahan_z(num_atoms);

  for (size_t frame_idx = 1; frame_idx < num_frames; ++frame_idx) {
    const size_t traj_frame = start_frame + frame_idx;
    const size_t traj_frame_prev = start_frame + frame_idx - 1;

    // Use Cell::minimumImage() for correct triclinic PBC handling.
    const bool use_pbc = (frames[traj_frame].volume() > 1e-9);

    const auto &curr_atoms = frames[traj_frame].atoms();
    const auto &prev_atoms = frames[traj_frame_prev].atoms();

    for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
      const math::Vector3<real_t> delta_r = curr_atoms[atom_idx].position() - prev_atoms[atom_idx].position();

      // Apply minimum image convention for correct unwrapping across PBC.
      const math::Vector3<real_t> min_delta_r = use_pbc ? frames[traj_frame].minimumImage(delta_r) : delta_r;

      kahan_x[atom_idx].add(min_delta_r.x());
      kahan_y[atom_idx].add(min_delta_r.y());
      kahan_z[atom_idx].add(min_delta_r.z());

      unwrapped[atom_idx][frame_idx] = {kahan_x[atom_idx].value(), kahan_y[atom_idx].value(),
                                        kahan_z[atom_idx].value()};
    }
  }

  // --- Compute MSD using Wiener-Khinchin Fast Correlation Algorithm ---
  // MSD(k) = (S1(k) - 2 * S2(k)) / (N - k)
  // S1(k) = sum_{i=0}^{N-1-k} (r_{i+k}^2 + r_i^2)
  // S2(k) = sum_{i=0}^{N-1-k} r_{i+k} * r_i

  std::vector<real_t> counts(max_correlation_frames + 1, 0.0);
  for (int lag = 0; lag <= max_correlation_frames; ++lag) {
    counts[lag] = static_cast<real_t>(num_frames - lag);
  }

  // Thread-local MSD accumulators to eliminate mutex contention
  tbb::enumerable_thread_specific<std::vector<real_t>> ets_msd(
      [&] { return std::vector<real_t>(max_correlation_frames + 1, 0.0); });

  // Thread-local FFT workspace reused across per-atom autocorrelate calls.
  tbb::enumerable_thread_specific<std::vector<std::complex<real_t>>> ets_ws_msd;

  tbb::parallel_for(static_cast<size_t>(0), num_atoms, [&](size_t atom_idx) {
    auto &local_msd = ets_msd.local();
    auto &workspace = ets_ws_msd.local();
    const auto &u_i = unwrapped[atom_idx];

    std::vector<real_t> pos_x(num_frames);
    std::vector<real_t> pos_y(num_frames);
    std::vector<real_t> pos_z(num_frames);
    std::vector<real_t> r_sq(num_frames);
    for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
      pos_x[frame_idx] = u_i[frame_idx].x();
      pos_y[frame_idx] = u_i[frame_idx].y();
      pos_z[frame_idx] = u_i[frame_idx].z();
      r_sq[frame_idx] = pos_x[frame_idx] * pos_x[frame_idx] + pos_y[frame_idx] * pos_y[frame_idx] +
                        pos_z[frame_idx] * pos_z[frame_idx];
    }

    auto S2_x = correlation::math::autocorrelate(pos_x, workspace);
    auto S2_y = correlation::math::autocorrelate(pos_y, workspace);
    auto S2_z = correlation::math::autocorrelate(pos_z, workspace);

    std::vector<real_t> sum_S1(max_correlation_frames + 1, 0.0);
    sum_S1[0] = static_cast<real_t>(static_cast<real_t>(2.0) * std::accumulate(r_sq.begin(), r_sq.end(), 0.0));
    for (int lag_idx = 1; lag_idx <= max_correlation_frames; ++lag_idx) {
      sum_S1[lag_idx] = sum_S1[lag_idx - 1] - r_sq[lag_idx - 1] - r_sq[num_frames - lag_idx];
    }

    for (int lag = 1; lag <= max_correlation_frames; ++lag) {
      real_t const sum_S2 = S2_x[lag] + S2_y[lag] + S2_z[lag];
      local_msd[lag] += (sum_S1[lag] - static_cast<real_t>(2.0) * sum_S2);
    }
  });

  std::vector<real_t> msd(max_correlation_frames + 1, 0.0);
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
      msd[lag] /= (counts[lag] * static_cast<real_t>(num_atoms));
    }
  }

  return msd;
}

std::vector<real_t> DynamicsAnalyzer::calculateNormalizedVACF(const correlation::core::Trajectory &traj,
                                                              MaxFrames max_correlation_frames, StartFrame start_frame,
                                                              EndFrame end_frame) {
  std::vector<real_t> vacf = calculateVACF(traj, max_correlation_frames, start_frame, end_frame);
  if (!vacf.empty() && vacf[0] != 0.0) {
    real_t const normalization_factor = static_cast<real_t>(1.0) / vacf[0];
    for (real_t &val : vacf) {
      val *= normalization_factor;
    }
  }
  return vacf;
}

std::tuple<std::vector<real_t>, std::vector<real_t>, std::vector<real_t>>
DynamicsAnalyzer::calculateVDOS(const std::vector<real_t> &vacf, real_t time_step) {
  if (vacf.empty()) {
    return {};
  }

  size_t num_frames = vacf.size();
  real_t const t_max = static_cast<real_t>(num_frames - 1) * time_step; // Total time of correlation

  // Define frequency range
  // Max frequency (Nyquist) is 1 / (2 * dt)
  // But we are interested in physical vibrations, usually up to ~100 THz or
  // ~3000 cm^-1 time_step is in fs. 1 fs = 10^-15 s. 1 THz = 10^12 Hz. Nyquist freq in
  // THz = 1 / (2 * time_step * 10^-15) * 10^-12 = 1000 / (2 * time_step)

  real_t nyquist_thz = static_cast<real_t>(1000.0) / (static_cast<real_t>(2.0) * time_step);

  // Number of frequency points. Let's make it high resolution.
  size_t const num_freq_points = 2000;
  real_t d_nu = nyquist_thz / num_freq_points; // Frequency step in THz

  std::vector<real_t> frequencies(num_freq_points);
  std::vector<real_t> intensities_real(num_freq_points);
  std::vector<real_t> intensities_imag(num_freq_points);

  // Pre-calculate windowed VACF
  std::vector<real_t> windowed_vacf(num_frames);
  // Define standard deviation for the Gaussian window (e.g., stopping at 3
  // sigma at t_max)
  real_t const sigma = t_max >= 1e-6 ? t_max / static_cast<real_t>(3.0) : static_cast<real_t>(1.0);
  for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
    real_t const time_val = static_cast<real_t>(frame_idx) * time_step;
    real_t const window = std::exp(-static_cast<real_t>(0.5) * std::pow(time_val / sigma, static_cast<real_t>(2.0)));
    windowed_vacf[frame_idx] = vacf[frame_idx] * window;
  }

  std::vector<size_t> freq_indices(num_freq_points);
  std::iota(freq_indices.begin(), freq_indices.end(), 0);

  tbb::parallel_for(static_cast<size_t>(0), num_freq_points, [&](size_t freq_idx) {
    real_t const freq_val = static_cast<real_t>(freq_idx) * d_nu; // Frequency in THz
    frequencies[freq_idx] = freq_val;

    // $\omega = 2 * \pi * \nu * 0.001$ (to handle THz to fs
    // scale)
    real_t const theta = correlation::math::two_pi * freq_val * time_step * static_cast<real_t>(0.001);

    auto [integral_real, integral_imag] = integrate_vdos_frequency(theta, windowed_vacf, time_step);

    real_t const damping =
        std::exp(-static_cast<real_t>(0.5) *
                 std::pow(freq_val / (static_cast<real_t>(0.5) * nyquist_thz), static_cast<real_t>(2.0)));
    intensities_real[freq_idx] = integral_real * damping;
    intensities_imag[freq_idx] = integral_imag * damping;
  });

  return {frequencies, intensities_real, intensities_imag};
}

real_t DynamicsAnalyzer::computeDiffusionCoefficientMSD(const std::vector<real_t> &time,
                                                        const std::vector<real_t> &msd) {
  if (time.size() < 2 || time.size() != msd.size()) {
    return 0.0;
  }
  // Fit on the second half of the data
  size_t const start_idx = time.size() / 2;
  size_t const num_points = time.size() - start_idx;
  if (num_points < 2) {
    return 0.0;
  }

  real_t sum_x = 0.0;
  real_t sum_y = 0.0;
  real_t sum_xx = 0.0;
  real_t sum_xy = 0.0;

  for (size_t time_idx = start_idx; time_idx < time.size(); ++time_idx) {
    real_t const time_val = time[time_idx];
    real_t const msd_val = msd[time_idx];
    sum_x += time_val;
    sum_y += msd_val;
    sum_xx += time_val * time_val;
    sum_xy += time_val * msd_val;
  }

  real_t const denominator = (static_cast<real_t>(num_points) * sum_xx - sum_x * sum_x);
  if (std::abs(denominator) < 1e-12) {
    return 0.0;
  }

  real_t const slope = (static_cast<real_t>(num_points) * sum_xy - sum_x * sum_y) / denominator;
  if (slope < 0.0) {
    return 0.0; // Self-diffusion coefficient cannot be negative
  }
  // D = slope / 6.0
  return slope / static_cast<real_t>(6.0);
}

real_t DynamicsAnalyzer::computeDiffusionCoefficientVACF(const std::vector<real_t> &time,
                                                         const std::vector<real_t> &vacf) {
  if (time.size() < 2 || time.size() != vacf.size()) {
    return 0.0;
  }
  real_t integral = 0.0;
  for (size_t time_idx = 0; time_idx < time.size() - 1; ++time_idx) {
    real_t const dt_step = time[time_idx + 1] - time[time_idx];
    if (dt_step <= 0.0) {
      return 0.0; // Time must be strictly increasing
    }
    integral += static_cast<real_t>(0.5) * (vacf[time_idx] + vacf[time_idx + 1]) * dt_step;
  }
  if (integral < 0.0) {
    return 0.0; // Self-diffusion coefficient cannot be negative
  }
  // Green-Kubo: D = 1/3 * integral
  return integral / static_cast<real_t>(3.0);
}

real_t DynamicsAnalyzer::computeRelaxationTime(const std::vector<real_t> &time,
                                               const std::vector<real_t> &normalized_vacf) {
  if (time.size() < 2 || time.size() != normalized_vacf.size()) {
    return 0.0;
  }
  real_t integral = 0.0;
  for (size_t time_idx = 0; time_idx < time.size() - 1; ++time_idx) {
    real_t const dt_step = time[time_idx + 1] - time[time_idx];
    if (dt_step <= 0.0) {
      return 0.0; // Time must be strictly increasing
    }
    integral += static_cast<real_t>(0.5) * (normalized_vacf[time_idx] + normalized_vacf[time_idx + 1]) * dt_step;
  }
  if (integral < 0.0) {
    return 0.0; // Relaxation time must be non-negative
  }
  return integral;
}

} // namespace correlation::analysis
