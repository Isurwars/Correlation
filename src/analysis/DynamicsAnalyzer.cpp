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
#include <cmath>
#include <complex>
#include <span>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <utility>

namespace correlation::analysis {

namespace {

[[nodiscard]] double evaluateVdosWindow(size_t frame_idx, size_t num_frames,
                                        VDOSWindowType window_type) noexcept {
  if (num_frames <= 1 || window_type == VDOSWindowType::None) {
    return 1.0;
  }
  const double ratio = static_cast<double>(frame_idx) / static_cast<double>(num_frames - 1);
  switch (window_type) {
  case VDOSWindowType::Hann:
    return 0.5 * (1.0 + std::cos(correlation::math::PI * ratio));
  case VDOSWindowType::Blackman:
    return 0.42 + 0.5 * std::cos(correlation::math::PI * ratio) +
           0.08 * std::cos(correlation::math::TWO_PI * ratio);
  case VDOSWindowType::Gaussian: {
    const double scaled = 3.0 * ratio;
    return std::exp(-0.5 * scaled * scaled);
  }
  case VDOSWindowType::None:
  default:
    return 1.0;
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

/// Contiguous Structure-of-Arrays trajectory buffer in double precision: [atom * num_frames + t].
struct TrajectorySoA {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  size_t num_atoms{0};
  size_t num_frames{0};

  [[nodiscard]] bool empty() const noexcept { return num_atoms == 0 || num_frames == 0; }

  [[nodiscard]] std::span<const double> atomX(size_t atom_idx) const noexcept {
    return {x.data() + (atom_idx * num_frames), num_frames};
  }
  [[nodiscard]] std::span<const double> atomY(size_t atom_idx) const noexcept {
    return {y.data() + (atom_idx * num_frames), num_frames};
  }
  [[nodiscard]] std::span<const double> atomZ(size_t atom_idx) const noexcept {
    return {z.data() + (atom_idx * num_frames), num_frames};
  }
};

struct ElementMassLayout {
  std::vector<std::string> symbols;
  std::vector<size_t> atom_elem_idx;
  std::vector<size_t> counts;
  std::vector<double> masses;
  double total_mass{0.0};
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
      double mass = 1.0;
      if (const auto m_opt = correlation::physics::findAtomicMass(sym);
          m_opt.has_value() && *m_opt > static_cast<real_t>(0.0)) {
        mass = static_cast<double>(*m_opt);
      }
      layout.masses.push_back(mass);
    }
    layout.atom_elem_idx[atom_idx] = elem_it->second;
    layout.counts[elem_it->second]++;
    layout.total_mass += layout.masses[elem_it->second];
  }
  return layout;
}

/// Zero-allocation double-precision FFT autocorrelation accumulating into out_accum[0..lag_len).
void autocorrelateAccumulate(std::span<const double> signal,
                             std::vector<std::complex<double>> &workspace,
                             std::span<double> out_accum) {
  const size_t num_pts = signal.size();
  if (num_pts == 0 || out_accum.empty()) {
    return;
  }

  const size_t fft_len = correlation::math::findNextGoodFFTSize(2 * num_pts - 1);
  workspace.assign(fft_len, {0.0, 0.0});
  for (size_t i = 0; i < num_pts; ++i) {
    workspace[i] = {signal[i], 0.0};
  }

  correlation::math::computeFFT(workspace, false);
  for (size_t i = 0; i < fft_len; ++i) {
    const double re_val = workspace[i].real();
    const double im_val = workspace[i].imag();
    workspace[i] = {re_val * re_val + im_val * im_val, 0.0};
  }
  correlation::math::computeFFT(workspace, true);

  const size_t count = std::min(out_accum.size(), num_pts);
  for (size_t lag = 0; lag < count; ++lag) {
    out_accum[lag] += workspace[lag].real();
  }
}

void subtractComDrift(TrajectorySoA &unwrapped, const ElementMassLayout &layout) {
  if (unwrapped.num_atoms <= 1 || layout.total_mass <= 0.0) {
    return;
  }
  const size_t num_atoms = unwrapped.num_atoms;
  const size_t num_frames = unwrapped.num_frames;
  const double inv_total_mass = 1.0 / layout.total_mass;

  // Use grain_size >= 64 (8+ 64-byte cache lines of double) to prevent false sharing across frames
  constexpr size_t K_FRAME_GRAIN_SIZE = 64;
  tbb::parallel_for(
      tbb::blocked_range<size_t>(1, num_frames, K_FRAME_GRAIN_SIZE),
      [&](const tbb::blocked_range<size_t> &range) {
        for (size_t frame_idx = range.begin(); frame_idx != range.end(); ++frame_idx) {
          double com_x = 0.0;
          double com_y = 0.0;
          double com_z = 0.0;
          for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
            const double mass = layout.masses[layout.atom_elem_idx[atom_idx]];
            const size_t offset = (atom_idx * num_frames) + frame_idx;
            com_x += unwrapped.x[offset] * mass;
            com_y += unwrapped.y[offset] * mass;
            com_z += unwrapped.z[offset] * mass;
          }
          com_x *= inv_total_mass;
          com_y *= inv_total_mass;
          com_z *= inv_total_mass;

          for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
            const size_t offset = (atom_idx * num_frames) + frame_idx;
            unwrapped.x[offset] -= com_x;
            unwrapped.y[offset] -= com_y;
            unwrapped.z[offset] -= com_z;
          }
        }
      });
}

[[nodiscard]] TrajectorySoA
computeUnwrappedTrajectories(const std::vector<correlation::core::Cell> &frames, size_t start_frame,
                             size_t num_frames, size_t num_atoms, const ElementMassLayout &layout,
                             bool remove_com_drift) {
  TrajectorySoA unwrapped;
  unwrapped.num_atoms = num_atoms;
  unwrapped.num_frames = num_frames;
  const size_t total_elems = num_atoms * num_frames;
  unwrapped.x.assign(total_elems, 0.0);
  unwrapped.y.assign(total_elems, 0.0);
  unwrapped.z.assign(total_elems, 0.0);

  for (size_t frame_idx = 1; frame_idx < num_frames; ++frame_idx) {
    const size_t traj_frame = start_frame + frame_idx;
    const size_t traj_frame_prev = start_frame + frame_idx - 1;
    const bool use_pbc = (frames[traj_frame].volume() > static_cast<real_t>(1e-9));

    const auto &curr_atoms = frames[traj_frame].atoms();
    const auto &prev_atoms = frames[traj_frame_prev].atoms();

    for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
      const math::Vector3<real_t> delta_r =
          curr_atoms[atom_idx].position() - prev_atoms[atom_idx].position();
      const math::Vector3<real_t> min_delta_r =
          use_pbc ? frames[traj_frame].minimumImage(delta_r) : delta_r;

      const size_t offset = (atom_idx * num_frames) + frame_idx;
      unwrapped.x[offset] = unwrapped.x[offset - 1] + static_cast<double>(min_delta_r.x());
      unwrapped.y[offset] = unwrapped.y[offset - 1] + static_cast<double>(min_delta_r.y());
      unwrapped.z[offset] = unwrapped.z[offset - 1] + static_cast<double>(min_delta_r.z());
    }
  }

  if (remove_com_drift) {
    subtractComDrift(unwrapped, layout);
  }
  return unwrapped;
}

[[nodiscard]] TrajectorySoA
computeCorrectedVelocities(const std::vector<correlation::core::Cell> &frames, size_t start_frame,
                           size_t num_frames, size_t num_atoms, const ElementMassLayout &layout) {
  if (layout.total_mass <= 0.0) {
    return {};
  }

  TrajectorySoA velocities;
  velocities.num_atoms = num_atoms;
  velocities.num_frames = num_frames;
  const size_t total_elems = num_atoms * num_frames;
  velocities.x.resize(total_elems);
  velocities.y.resize(total_elems);
  velocities.z.resize(total_elems);

  const double inv_total_mass = 1.0 / layout.total_mass;

  constexpr size_t K_FRAME_GRAIN_SIZE = 64;
  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, num_frames, K_FRAME_GRAIN_SIZE),
      [&](const tbb::blocked_range<size_t> &range) {
        for (size_t frame_idx = range.begin(); frame_idx != range.end(); ++frame_idx) {
          const auto &frame = frames[start_frame + frame_idx];
          double mom_x = 0.0;
          double mom_y = 0.0;
          double mom_z = 0.0;
          for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
            const double mass = layout.masses[layout.atom_elem_idx[atom_idx]];
            const auto vel = frame.atomVelocity(atom_idx);
            mom_x += static_cast<double>(vel.x()) * mass;
            mom_y += static_cast<double>(vel.y()) * mass;
            mom_z += static_cast<double>(vel.z()) * mass;
          }
          const double com_vx = mom_x * inv_total_mass;
          const double com_vy = mom_y * inv_total_mass;
          const double com_vz = mom_z * inv_total_mass;

          for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
            const auto vel = frame.atomVelocity(atom_idx);
            const size_t offset = (atom_idx * num_frames) + frame_idx;
            velocities.x[offset] = static_cast<double>(vel.x()) - com_vx;
            velocities.y[offset] = static_cast<double>(vel.y()) - com_vy;
            velocities.z[offset] = static_cast<double>(vel.z()) - com_vz;
          }
        }
      });

  return velocities;
}

[[nodiscard]] PartialVACFResult
assemblePartialVACFResult(const std::vector<std::vector<double>> &elem_sums,
                          const ElementMassLayout &layout, size_t num_frames, size_t num_atoms) {
  const size_t num_elems = layout.symbols.size();
  const size_t lag_len = elem_sums.empty() ? 0 : elem_sums[0].size();
  PartialVACFResult result;
  std::vector<double> total_unweighted_d(lag_len, 0.0);
  std::vector<double> total_mass_weighted_d(lag_len, 0.0);

  for (size_t elem_idx = 0; elem_idx < num_elems; ++elem_idx) {
    std::vector<real_t> elem_vacf(lag_len, static_cast<real_t>(0.0));
    const auto count_e = static_cast<double>(layout.counts[elem_idx]);
    const double weight_e =
        layout.total_mass > 0.0 ? (count_e * layout.masses[elem_idx]) / layout.total_mass : 0.0;
    result.element_mass_weights[layout.symbols[elem_idx]] = static_cast<real_t>(weight_e);

    for (size_t lag = 0; lag < lag_len; ++lag) {
      const auto n_origins = static_cast<double>(num_frames - lag);
      const double val = elem_sums[elem_idx][lag] / (n_origins * count_e);
      elem_vacf[lag] = static_cast<real_t>(val);
      total_unweighted_d[lag] +=
          elem_sums[elem_idx][lag] / (n_origins * static_cast<double>(num_atoms));
      total_mass_weighted_d[lag] += weight_e * val;
    }
    result.element_vacfs[layout.symbols[elem_idx]] = std::move(elem_vacf);
  }

  result.total_unweighted.resize(lag_len);
  result.total_mass_weighted.resize(lag_len);
  for (size_t lag = 0; lag < lag_len; ++lag) {
    result.total_unweighted[lag] = static_cast<real_t>(total_unweighted_d[lag]);
    result.total_mass_weighted[lag] = static_cast<real_t>(total_mass_weighted_d[lag]);
  }
  return result;
}

[[nodiscard]] PartialMSDResult
assemblePartialMSDResult(const std::vector<std::vector<double>> &elem_sums,
                         const ElementMassLayout &layout, size_t num_frames, size_t num_atoms) {
  const size_t num_elems = layout.symbols.size();
  const size_t lag_len = elem_sums.empty() ? 0 : elem_sums[0].size();
  PartialMSDResult result;
  std::vector<double> total_d(lag_len, 0.0);

  for (size_t elem_idx = 0; elem_idx < num_elems; ++elem_idx) {
    std::vector<real_t> elem_msd(lag_len, static_cast<real_t>(0.0));
    const auto count_e = static_cast<double>(layout.counts[elem_idx]);
    for (size_t lag = 1; lag < lag_len; ++lag) {
      const auto n_origins = static_cast<double>(num_frames - lag);
      const double val = std::max(0.0, elem_sums[elem_idx][lag] / (n_origins * count_e));
      elem_msd[lag] = static_cast<real_t>(val);
      total_d[lag] += elem_sums[elem_idx][lag] / (n_origins * static_cast<double>(num_atoms));
    }
    result.element_msds[layout.symbols[elem_idx]] = std::move(elem_msd);
  }

  result.total.assign(lag_len, static_cast<real_t>(0.0));
  for (size_t lag = 1; lag < lag_len; ++lag) {
    result.total[lag] = static_cast<real_t>(std::max(0.0, total_d[lag]));
  }
  return result;
}

struct MSDThreadScratch {
  std::vector<std::complex<double>> fft_ws;
  std::vector<double> r_sq;
  std::vector<double> s2_sum;
};

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
  const size_t num_atoms = frames[start_frame].atoms().size();
  if (num_atoms == 0) {
    return {};
  }
  for (size_t frame_idx = start_frame + 1; frame_idx < end_frame; ++frame_idx) {
    if (frames[frame_idx].atoms().size() != num_atoms) {
      throw std::invalid_argument(
          "DynamicsAnalyzer requires a constant atom count across analyzed frames (frame " +
          std::to_string(start_frame) + " has " + std::to_string(num_atoms) + " atoms, frame " +
          std::to_string(frame_idx) + " has " + std::to_string(frames[frame_idx].atoms().size()) +
          " atoms).");
    }
  }
  if (max_correlation_frames < 0) {
    max_correlation_frames = std::max(1, static_cast<int>(num_frames / 2));
  } else if (std::cmp_greater_equal(max_correlation_frames, num_frames)) {
    max_correlation_frames = static_cast<int>(num_frames) - 1;
  }
  if (max_correlation_frames <= 0) {
    return {};
  }

  const ElementMassLayout layout = buildElementMassLayout(frames[start_frame].atoms());
  const TrajectorySoA atom_velocities =
      computeCorrectedVelocities(frames, start_frame, num_frames, num_atoms, layout);
  if (atom_velocities.empty()) {
    return {};
  }

  const size_t num_elems = layout.symbols.size();
  const size_t lag_len = static_cast<size_t>(max_correlation_frames) + 1;

  tbb::enumerable_thread_specific<std::vector<std::vector<double>>> ets_elem_vacf([&] {
    return std::vector<std::vector<double>>(num_elems, std::vector<double>(lag_len, 0.0));
  });
  tbb::enumerable_thread_specific<std::vector<std::complex<double>>> ets_ws;

  tbb::parallel_for(static_cast<size_t>(0), num_atoms, [&](size_t atom_idx) {
    auto &local_elem = ets_elem_vacf.local()[layout.atom_elem_idx[atom_idx]];
    auto &workspace = ets_ws.local();
    const std::span<double> accum_span(local_elem.data(), lag_len);
    autocorrelateAccumulate(atom_velocities.atomX(atom_idx), workspace, accum_span);
    autocorrelateAccumulate(atom_velocities.atomY(atom_idx), workspace, accum_span);
    autocorrelateAccumulate(atom_velocities.atomZ(atom_idx), workspace, accum_span);
  });

  std::vector<std::vector<double>> elem_sums(num_elems, std::vector<double>(lag_len, 0.0));
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

PartialMSDResult DynamicsAnalyzer::calculatePartialMSD(const correlation::core::Trajectory &traj,
                                                       MaxFrames max_correlation_frames_arg,
                                                       StartFrame start_frame_arg,
                                                       EndFrame end_frame_arg,
                                                       bool remove_com_drift) {
  int max_correlation_frames = max_correlation_frames_arg.value;
  if (traj.getFrameCount() < 2) {
    return {};
  }
  const auto &frames = traj.getFrames();
  const size_t total_frames = frames.size();
  if (total_frames < 2) {
    return {};
  }

  const size_t start_frame = std::min(start_frame_arg.value, total_frames - 1);
  const size_t end_frame = std::min(end_frame_arg.value, total_frames);
  if (start_frame >= end_frame || (end_frame - start_frame) < 2) {
    return {};
  }

  const size_t num_atoms = frames[start_frame].atoms().size();
  if (num_atoms == 0) {
    return {};
  }
  for (size_t frame_idx = start_frame + 1; frame_idx < end_frame; ++frame_idx) {
    if (frames[frame_idx].atoms().size() != num_atoms) {
      throw std::invalid_argument(
          "DynamicsAnalyzer requires a constant atom count across analyzed frames (frame " +
          std::to_string(start_frame) + " has " + std::to_string(num_atoms) + " atoms, frame " +
          std::to_string(frame_idx) + " has " + std::to_string(frames[frame_idx].atoms().size()) +
          " atoms).");
    }
  }

  const size_t num_frames = end_frame - start_frame;
  if (max_correlation_frames < 0 || std::cmp_greater_equal(max_correlation_frames, num_frames)) {
    max_correlation_frames = static_cast<int>(num_frames / 2);
  }
  if (max_correlation_frames <= 0) {
    return {};
  }

  const ElementMassLayout layout = buildElementMassLayout(frames[start_frame].atoms());
  const TrajectorySoA unwrapped = computeUnwrappedTrajectories(frames, start_frame, num_frames,
                                                               num_atoms, layout, remove_com_drift);

  const size_t num_elems = layout.symbols.size();
  const size_t lag_len = static_cast<size_t>(max_correlation_frames) + 1;

  tbb::enumerable_thread_specific<std::vector<std::vector<double>>> ets_elem_msd([&] {
    return std::vector<std::vector<double>>(num_elems, std::vector<double>(lag_len, 0.0));
  });
  tbb::enumerable_thread_specific<MSDThreadScratch> ets_scratch;

  tbb::parallel_for(static_cast<size_t>(0), num_atoms, [&](size_t atom_idx) {
    auto &local_elem = ets_elem_msd.local()[layout.atom_elem_idx[atom_idx]];
    auto &scratch = ets_scratch.local();
    const auto pos_x = unwrapped.atomX(atom_idx);
    const auto pos_y = unwrapped.atomY(atom_idx);
    const auto pos_z = unwrapped.atomZ(atom_idx);

    scratch.r_sq.resize(num_frames);
    double sum_r_sq = 0.0;
    for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
      const double r2_val = (pos_x[frame_idx] * pos_x[frame_idx]) +
                            (pos_y[frame_idx] * pos_y[frame_idx]) +
                            (pos_z[frame_idx] * pos_z[frame_idx]);
      scratch.r_sq[frame_idx] = r2_val;
      sum_r_sq += r2_val;
    }

    scratch.s2_sum.assign(lag_len, 0.0);
    const std::span<double> s2_span(scratch.s2_sum.data(), lag_len);
    autocorrelateAccumulate(pos_x, scratch.fft_ws, s2_span);
    autocorrelateAccumulate(pos_y, scratch.fft_ws, s2_span);
    autocorrelateAccumulate(pos_z, scratch.fft_ws, s2_span);

    double s1_running = 2.0 * sum_r_sq;
    for (size_t lag = 1; lag < lag_len; ++lag) {
      s1_running -= scratch.r_sq[lag - 1] + scratch.r_sq[num_frames - lag];
      local_elem[lag] += std::max(0.0, s1_running - (2.0 * scratch.s2_sum[lag]));
    }
  });

  std::vector<std::vector<double>> elem_sums(num_elems, std::vector<double>(lag_len, 0.0));
  for (const auto &local_all : ets_elem_msd) {
    for (size_t elem_idx = 0; elem_idx < num_elems; ++elem_idx) {
      for (size_t lag = 1; lag < lag_len; ++lag) {
        elem_sums[elem_idx][lag] += local_all[elem_idx][lag];
      }
    }
  }

  return assemblePartialMSDResult(elem_sums, layout, num_frames, num_atoms);
}

std::vector<real_t> DynamicsAnalyzer::calculateMSD(const correlation::core::Trajectory &traj,
                                                   MaxFrames max_correlation_frames_arg,
                                                   StartFrame start_frame_arg,
                                                   EndFrame end_frame_arg, bool remove_com_drift) {
  return calculatePartialMSD(traj, max_correlation_frames_arg, start_frame_arg, end_frame_arg,
                             remove_com_drift)
      .total;
}

std::vector<real_t>
DynamicsAnalyzer::calculateNormalizedVACF(const correlation::core::Trajectory &traj,
                                          MaxFrames max_correlation_frames, StartFrame start_frame,
                                          EndFrame end_frame) {
  std::vector<real_t> vacf = calculateVACF(traj, max_correlation_frames, start_frame, end_frame);
  if (!vacf.empty() && vacf[0] != static_cast<real_t>(0.0)) {
    const double inv_c0 = 1.0 / static_cast<double>(vacf[0]);
    for (real_t &val : vacf) {
      val = static_cast<real_t>(static_cast<double>(val) * inv_c0);
    }
  }
  return vacf;
}

std::pair<std::vector<real_t>, std::vector<real_t>>
DynamicsAnalyzer::calculateVDOS(const std::vector<real_t> &vacf, real_t time_step,
                                VDOSWindowType window_type, size_t zero_pad_factor) {
  if (vacf.empty() || time_step <= static_cast<real_t>(0.0) || !std::isfinite(time_step)) {
    return {};
  }

  const size_t num_frames = vacf.size();
  const size_t fft_len = computeEvenPaddedFFTSize(num_frames, zero_pad_factor);
  const size_t num_freq_points = fft_len / 2 + 1;

  // Build exact even-symmetric zero-padded sequence in double precision for numerical stability.
  // Re(DFT[j]) = x_0 + 2 * sum_{k=1}^{M-1} x_k * cos(2*pi*j*k / L), giving the exact trapezoidal
  // cosine transform with 1/2 endpoint weight at t = 0 and zero imaginary component.
  std::vector<std::complex<double>> buffer(fft_len, {0.0, 0.0});
  buffer[0] = {static_cast<double>(vacf[0]) * evaluateVdosWindow(0, num_frames, window_type), 0.0};

  const size_t half_len = fft_len / 2;
  for (size_t k = 1; k < num_frames; ++k) {
    const double val =
        static_cast<double>(vacf[k]) * evaluateVdosWindow(k, num_frames, window_type);
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
    intensities[j] = static_cast<real_t>(std::isfinite(raw_val) ? std::max(0.0, raw_val) : 0.0);
  }

  return {std::move(frequencies), std::move(intensities)};
}

real_t DynamicsAnalyzer::computeDiffusionCoefficientMSD(const std::vector<real_t> &time,
                                                        const std::vector<real_t> &msd) {
  if (time.size() < 2 || time.size() != msd.size()) {
    return static_cast<real_t>(0.0);
  }
  // Centered OLS in double precision on the second half to prevent catastrophic cancellation.
  const size_t start_idx = time.size() / 2;
  const size_t num_points = time.size() - start_idx;
  if (num_points < 2) {
    return static_cast<real_t>(0.0);
  }

  double mean_x = 0.0;
  double mean_y = 0.0;
  for (size_t idx = start_idx; idx < time.size(); ++idx) {
    mean_x += static_cast<double>(time[idx]);
    mean_y += static_cast<double>(msd[idx]);
  }
  const auto inv_n = 1.0 / static_cast<double>(num_points);
  mean_x *= inv_n;
  mean_y *= inv_n;

  double var_x = 0.0;
  double cov_xy = 0.0;
  for (size_t idx = start_idx; idx < time.size(); ++idx) {
    const double diff_x = static_cast<double>(time[idx]) - mean_x;
    const double diff_y = static_cast<double>(msd[idx]) - mean_y;
    var_x += diff_x * diff_x;
    cov_xy += diff_x * diff_y;
  }

  if (var_x < 1e-18 || !std::isfinite(var_x) || !std::isfinite(cov_xy)) {
    return static_cast<real_t>(0.0);
  }

  const double slope = cov_xy / var_x;
  if (slope < 0.0 || !std::isfinite(slope)) {
    return static_cast<real_t>(0.0); // Self-diffusion coefficient cannot be negative
  }
  // D = slope / 6.0
  return static_cast<real_t>(slope / 6.0);
}

real_t DynamicsAnalyzer::computeDiffusionCoefficientVACF(const std::vector<real_t> &time,
                                                         const std::vector<real_t> &vacf) {
  if (time.size() < 2 || time.size() != vacf.size()) {
    return static_cast<real_t>(0.0);
  }
  double integral = 0.0;
  for (size_t time_idx = 0; time_idx < time.size() - 1; ++time_idx) {
    const double dt_step =
        static_cast<double>(time[time_idx + 1]) - static_cast<double>(time[time_idx]);
    if (dt_step <= 0.0 || !std::isfinite(dt_step)) {
      return static_cast<real_t>(0.0); // Time must be strictly increasing
    }
    integral += 0.5 *
                (static_cast<double>(vacf[time_idx]) + static_cast<double>(vacf[time_idx + 1])) *
                dt_step;
  }
  if (integral < 0.0 || !std::isfinite(integral)) {
    return static_cast<real_t>(0.0); // Self-diffusion coefficient cannot be negative
  }
  // Green-Kubo: D = 1/3 * integral
  return static_cast<real_t>(integral / 3.0);
}

real_t DynamicsAnalyzer::computeRelaxationTime(const std::vector<real_t> &time,
                                               const std::vector<real_t> &normalized_vacf) {
  if (time.size() < 2 || time.size() != normalized_vacf.size()) {
    return static_cast<real_t>(0.0);
  }
  double integral = 0.0;
  for (size_t time_idx = 0; time_idx < time.size() - 1; ++time_idx) {
    const double dt_step =
        static_cast<double>(time[time_idx + 1]) - static_cast<double>(time[time_idx]);
    if (dt_step <= 0.0 || !std::isfinite(dt_step)) {
      return static_cast<real_t>(0.0); // Time must be strictly increasing
    }
    integral += 0.5 *
                (static_cast<double>(normalized_vacf[time_idx]) +
                 static_cast<double>(normalized_vacf[time_idx + 1])) *
                dt_step;
  }
  if (integral < 0.0 || !std::isfinite(integral)) {
    return static_cast<real_t>(0.0); // Relaxation time must be non-negative
  }
  return static_cast<real_t>(integral);
}

} // namespace correlation::analysis
