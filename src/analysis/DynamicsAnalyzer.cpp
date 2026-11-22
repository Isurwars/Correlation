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

[[nodiscard]] real_t evaluateVdosWindow(size_t frame_idx, size_t num_frames,
                                        VDOSWindowType window_type) noexcept {
  if (num_frames <= 1 || window_type == VDOSWindowType::None) {
    return static_cast<real_t>(1.0);
  }
  const double ratio = static_cast<double>(frame_idx) / static_cast<double>(num_frames - 1);
  switch (window_type) {
  case VDOSWindowType::Hann:
    return static_cast<real_t>(0.5 * (1.0 + std::cos(correlation::math::PI * ratio)));
  case VDOSWindowType::Blackman:
    return static_cast<real_t>(0.42 + 0.5 * std::cos(correlation::math::PI * ratio) +
                               0.08 * std::cos(correlation::math::TWO_PI * ratio));
  case VDOSWindowType::Gaussian: {
    const double scaled = 3.0 * ratio;
    return static_cast<real_t>(std::exp(-0.5 * scaled * scaled));
  }
  case VDOSWindowType::None:
  default:
    return static_cast<real_t>(1.0);
  }
}

[[nodiscard]] size_t computeEvenPaddedFFTSize(size_t num_frames, size_t zero_pad_factor) noexcept {
  const size_t pad_mult = std::max<size_t>(1, zero_pad_factor);
  const size_t min_sym_len = (num_frames > 1) ? (2 * num_frames - 2) : 2;
  const auto target = std::max<size_t>(
      {static_cast<size_t>(2048), min_sym_len, static_cast<size_t>(2) * pad_mult * num_frames});
  size_t fft_len = correlation::math::findNextGoodFFTSize(target);
  if ((fft_len % 2) != 0) {
    fft_len *= 2;
  }
  return fft_len;
}

std::vector<std::vector<correlation::math::Vector3<real_t>>>
computeUnwrappedTrajectories(const std::vector<correlation::core::Cell> &frames, size_t start_frame,
                             size_t num_frames, size_t num_atoms) {
  std::vector<std::vector<correlation::math::Vector3<real_t>>> unwrapped(
      num_atoms, std::vector<correlation::math::Vector3<real_t>>(num_frames, {0.0, 0.0, 0.0}));

  std::vector<correlation::KahanAccumulator<real_t>> kahan_x(num_atoms);
  std::vector<correlation::KahanAccumulator<real_t>> kahan_y(num_atoms);
  std::vector<correlation::KahanAccumulator<real_t>> kahan_z(num_atoms);

  for (size_t frame_idx = 1; frame_idx < num_frames; ++frame_idx) {
    const size_t traj_frame = start_frame + frame_idx;
    const size_t traj_frame_prev = start_frame + frame_idx - 1;
    const bool use_pbc = (frames[traj_frame].volume() > 1e-9);

    const auto &curr_atoms = frames[traj_frame].atoms();
    const auto &prev_atoms = frames[traj_frame_prev].atoms();

    for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
      const math::Vector3<real_t> delta_r =
          curr_atoms[atom_idx].position() - prev_atoms[atom_idx].position();
      const math::Vector3<real_t> min_delta_r =
          use_pbc ? frames[traj_frame].minimumImage(delta_r) : delta_r;

      kahan_x[atom_idx].add(min_delta_r.x());
      kahan_y[atom_idx].add(min_delta_r.y());
      kahan_z[atom_idx].add(min_delta_r.z());

      unwrapped[atom_idx][frame_idx] = {kahan_x[atom_idx].value(), kahan_y[atom_idx].value(),
                                        kahan_z[atom_idx].value()};
    }
  }
  return unwrapped;
}

std::vector<std::vector<correlation::math::Vector3<real_t>>>
computeCorrectedVelocities(const std::vector<correlation::core::Cell> &frames, size_t start_frame,
                           size_t num_frames, size_t num_atoms) {
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

  if (total_mass <= 0.0) {
    return {};
  }

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

  return atom_velocities;
}
struct ElementMassLayout {
  std::vector<std::string> symbols;
  std::vector<size_t> atom_elem_idx;
  std::vector<size_t> counts;
  std::vector<real_t> masses;
  real_t total_mass{0.0};
};

[[nodiscard]] ElementMassLayout
buildElementMassLayout(const std::vector<correlation::core::Atom> &ref_atoms) {
  ElementMassLayout layout;
  const size_t num_atoms = ref_atoms.size();
  layout.atom_elem_idx.resize(num_atoms);
  std::map<std::string, size_t> elem_to_idx;

  for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
    const std::string &sym = ref_atoms[atom_idx].element().symbol;
    auto [elem_it, inserted] = elem_to_idx.emplace(sym, layout.symbols.size());
    if (inserted) {
      layout.symbols.push_back(sym);
      layout.counts.push_back(0);
      auto mass = static_cast<real_t>(1.0);
      if (const auto m_opt = correlation::physics::findAtomicMass(sym);
          m_opt.has_value() && *m_opt > static_cast<real_t>(0.0)) {
        mass = *m_opt;
      }
      layout.masses.push_back(mass);
    }
    layout.atom_elem_idx[atom_idx] = elem_it->second;
    layout.counts[elem_it->second]++;
    layout.total_mass += layout.masses[elem_it->second];
  }
  return layout;
}

[[nodiscard]] PartialVACFResult
assemblePartialVACFResult(const std::vector<std::vector<real_t>> &elem_sums,
                          const ElementMassLayout &layout, size_t num_frames, size_t num_atoms) {
  const size_t num_elems = layout.symbols.size();
  const size_t lag_len = elem_sums.empty() ? 0 : elem_sums[0].size();
  PartialVACFResult result;
  result.total_unweighted.assign(lag_len, static_cast<real_t>(0.0));
  result.total_mass_weighted.assign(lag_len, static_cast<real_t>(0.0));

  for (size_t elem_idx = 0; elem_idx < num_elems; ++elem_idx) {
    std::vector<real_t> elem_vacf(lag_len, static_cast<real_t>(0.0));
    const auto count_e = static_cast<real_t>(layout.counts[elem_idx]);
    const real_t weight_e = layout.total_mass > static_cast<real_t>(0.0)
                                ? (count_e * layout.masses[elem_idx]) / layout.total_mass
                                : static_cast<real_t>(0.0);
    result.element_mass_weights[layout.symbols[elem_idx]] = weight_e;

    for (size_t lag = 0; lag < lag_len; ++lag) {
      const auto n_origins = static_cast<real_t>(num_frames - lag);
      const real_t val = elem_sums[elem_idx][lag] / (n_origins * count_e);
      elem_vacf[lag] = val;
      result.total_unweighted[lag] +=
          elem_sums[elem_idx][lag] / (n_origins * static_cast<real_t>(num_atoms));
      result.total_mass_weighted[lag] += weight_e * val;
    }
    result.element_vacfs[layout.symbols[elem_idx]] = std::move(elem_vacf);
  }
  return result;
}

} // namespace

PartialVACFResult DynamicsAnalyzer::calculatePartialVACF(const correlation::core::Trajectory &traj,
                                                         MaxFrames max_correlation_frames_arg,
                                                         StartFrame start_frame_arg,
                                                         EndFrame end_frame_arg) {
  int max_correlation_frames = max_correlation_frames_arg.value;
  if (traj.getFrameCount() < 2) {
    return {};
  }
  const auto &frames = traj.getFrames();
  const size_t total_frames = frames.size();
  if (total_frames < 2 || frames[0].atoms().empty()) {
    return {};
  }

  const size_t start_frame = std::min(start_frame_arg.value, total_frames - 1);
  const size_t end_frame = std::min(end_frame_arg.value, total_frames);
  if (start_frame >= end_frame || (end_frame - start_frame) < 2) {
    return {};
  }

  const size_t num_frames = end_frame - start_frame;
  const size_t num_atoms = frames[0].atoms().size();
  if (max_correlation_frames < 0) {
    max_correlation_frames = std::max(1, static_cast<int>(num_frames / 2));
  } else if (std::cmp_greater_equal(max_correlation_frames, num_frames)) {
    max_correlation_frames = static_cast<int>(num_frames) - 1;
  }
  if (max_correlation_frames <= 0) {
    return {};
  }

  const auto atom_velocities =
      computeCorrectedVelocities(frames, start_frame, num_frames, num_atoms);
  if (atom_velocities.empty()) {
    return {};
  }

  const ElementMassLayout layout = buildElementMassLayout(frames[start_frame].atoms());
  const size_t num_elems = layout.symbols.size();
  const size_t lag_len = static_cast<size_t>(max_correlation_frames) + 1;

  tbb::enumerable_thread_specific<std::vector<std::vector<real_t>>> ets_elem_vacf([&] {
    return std::vector<std::vector<real_t>>(num_elems, std::vector<real_t>(lag_len, 0.0));
  });
  tbb::enumerable_thread_specific<std::vector<std::complex<real_t>>> ets_ws;

  tbb::parallel_for(static_cast<size_t>(0), num_atoms, [&](size_t atom_idx) {
    auto &local_elem = ets_elem_vacf.local()[layout.atom_elem_idx[atom_idx]];
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

    const auto s2_x = correlation::math::autocorrelate(vel_x, workspace);
    const auto s2_y = correlation::math::autocorrelate(vel_y, workspace);
    const auto s2_z = correlation::math::autocorrelate(vel_z, workspace);
    for (size_t lag = 0; lag < lag_len; ++lag) {
      local_elem[lag] += s2_x[lag] + s2_y[lag] + s2_z[lag];
    }
  });

  std::vector<std::vector<real_t>> elem_sums(num_elems, std::vector<real_t>(lag_len, 0.0));
  for (const auto &local_all : ets_elem_vacf) {
    for (size_t elem_idx = 0; elem_idx < num_elems; ++elem_idx) {
      for (size_t lag = 0; lag < lag_len; ++lag) {
        elem_sums[elem_idx][lag] += local_all[elem_idx][lag];
      }
    }
  }

  return assemblePartialVACFResult(elem_sums, layout, num_frames, num_atoms);
}

std::vector<real_t> DynamicsAnalyzer::calculateVACF(const correlation::core::Trajectory &traj,
                                                    MaxFrames max_correlation_frames_arg,
                                                    StartFrame start_frame_arg,
                                                    EndFrame end_frame_arg) {
  return calculatePartialVACF(traj, max_correlation_frames_arg, start_frame_arg, end_frame_arg)
      .total_unweighted;
}

std::vector<real_t> DynamicsAnalyzer::calculateMSD(const correlation::core::Trajectory &traj,
                                                   MaxFrames max_correlation_frames_arg,
                                                   StartFrame start_frame_arg,
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
  if (num_atoms == 0) {
    return {};
  }

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
  const auto unwrapped = computeUnwrappedTrajectories(frames, start_frame, num_frames, num_atoms);

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

    auto s2_x = correlation::math::autocorrelate(pos_x, workspace);
    auto s2_y = correlation::math::autocorrelate(pos_y, workspace);
    auto s2_z = correlation::math::autocorrelate(pos_z, workspace);

    std::vector<real_t> sum_s1(max_correlation_frames + 1, 0.0);
    sum_s1[0] = static_cast<real_t>(static_cast<real_t>(2.0) *
                                    std::accumulate(r_sq.begin(), r_sq.end(), 0.0));
    for (int lag_idx = 1; lag_idx <= max_correlation_frames; ++lag_idx) {
      sum_s1[lag_idx] = sum_s1[lag_idx - 1] - r_sq[lag_idx - 1] - r_sq[num_frames - lag_idx];
    }

    for (int lag = 1; lag <= max_correlation_frames; ++lag) {
      real_t const sum_s2 = s2_x[lag] + s2_y[lag] + s2_z[lag];
      local_msd[lag] += (sum_s1[lag] - static_cast<real_t>(2.0) * sum_s2);
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

std::vector<real_t>
DynamicsAnalyzer::calculateNormalizedVACF(const correlation::core::Trajectory &traj,
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

std::pair<std::vector<real_t>, std::vector<real_t>>
DynamicsAnalyzer::calculateVDOS(const std::vector<real_t> &vacf, real_t time_step,
                                VDOSWindowType window_type, size_t zero_pad_factor) {
  if (vacf.empty() || time_step <= static_cast<real_t>(0.0)) {
    return {};
  }

  const size_t num_frames = vacf.size();
  const size_t fft_len = computeEvenPaddedFFTSize(num_frames, zero_pad_factor);
  const size_t num_freq_points = fft_len / 2 + 1;

  // Build exact even-symmetric zero-padded sequence in double precision for numerical stability.
  // Re(DFT[j]) = x_0 + 2 * sum_{k=1}^{M-1} x_k * cos(2*pi*j*k / L), giving the exact trapezoidal
  // cosine transform with 1/2 endpoint weight at t = 0 and zero imaginary component.
  std::vector<std::complex<double>> buffer(fft_len, {0.0, 0.0});
  buffer[0] = {static_cast<double>(vacf[0]) *
                   static_cast<double>(evaluateVdosWindow(0, num_frames, window_type)),
               0.0};

  const size_t half_len = fft_len / 2;
  for (size_t k = 1; k < num_frames; ++k) {
    const double val = static_cast<double>(vacf[k]) *
                       static_cast<double>(evaluateVdosWindow(k, num_frames, window_type));
    buffer[k] = {val, 0.0};
    if (k < half_len) {
      buffer[fft_len - k] = {val, 0.0};
    }
  }

  correlation::math::computeFFT(buffer, false);

  std::vector<real_t> frequencies(num_freq_points);
  std::vector<real_t> intensities(num_freq_points);

  // dt in fs -> 1000.0 / (L * dt) is d_nu in THz.
  const auto dt_d = static_cast<double>(time_step);
  const double d_nu = 1000.0 / (static_cast<double>(fft_len) * dt_d);

  for (size_t j = 0; j < num_freq_points; ++j) {
    frequencies[j] = static_cast<real_t>(static_cast<double>(j) * d_nu);
    // Clamp tiny negative numerical ripples from finite-time truncation to 0.
    const double raw_val = buffer[j].real() * dt_d;
    intensities[j] = static_cast<real_t>(std::max(0.0, raw_val));
  }

  return {std::move(frequencies), std::move(intensities)};
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
    integral += static_cast<real_t>(0.5) *
                (normalized_vacf[time_idx] + normalized_vacf[time_idx + 1]) * dt_step;
  }
  if (integral < 0.0) {
    return 0.0; // Relaxation time must be non-negative
  }
  return integral;
}

} // namespace correlation::analysis
