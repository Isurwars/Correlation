/**
 * @file VelocityIntegrator.cpp
 * @brief Implementation of finite-difference atomic velocity calculation under PBC.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/VelocityIntegrator.hpp"
#include "core/Trajectory.hpp"

namespace correlation::core {

namespace {

[[nodiscard]] math::Vector3<real_t>
computeDisplacement(const Cell &frame, const math::Vector3<real_t> &target_pos,
                    const math::Vector3<real_t> &source_pos) {
  const math::Vector3<real_t> delta_r = target_pos - source_pos;
  const bool use_pbc = (frame.volume() > static_cast<real_t>(1e-9));
  return use_pbc ? frame.minimumImage(delta_r) : delta_r;
}

void computeForwardBoundary(std::span<Cell> frames, real_t time_step) {
  auto &frame0 = frames[0];
  const auto &frame1 = frames[1];
  const size_t num_atoms = frame0.atomCount();
  for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
    const auto &pos0 = frame0.atoms()[atom_idx].position();
    const auto &pos1 = frame1.atoms()[atom_idx].position();
    frame0.atoms()[atom_idx].setVelocity(computeDisplacement(frame0, pos1, pos0) / time_step);
  }
}

void computeBackwardBoundary(std::span<Cell> frames, real_t time_step) {
  const size_t num_frames = frames.size();
  auto &frame_last = frames[num_frames - 1];
  const auto &frame_prev = frames[num_frames - 2];
  const size_t num_atoms = frame_last.atomCount();
  for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
    const auto &pos_n = frame_last.atoms()[atom_idx].position();
    const auto &pos_n_1 = frame_prev.atoms()[atom_idx].position();
    frame_last.atoms()[atom_idx].setVelocity(computeDisplacement(frame_last, pos_n, pos_n_1) /
                                             time_step);
  }
}

void computeCentralDifferences(std::span<Cell> frames, real_t time_step) {
  const size_t num_frames = frames.size();
  const real_t two_dt = static_cast<real_t>(2.0) * time_step;
  for (size_t frame_idx = 1; frame_idx + 1 < num_frames; ++frame_idx) {
    auto &current_frame = frames[frame_idx];
    const auto &frame_next = frames[frame_idx + 1];
    const auto &frame_prev = frames[frame_idx - 1];
    const size_t num_atoms = current_frame.atomCount();

    for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
      const auto &pos_next = frame_next.atoms()[atom_idx].position();
      const auto &pos_prev = frame_prev.atoms()[atom_idx].position();
      current_frame.atoms()[atom_idx].setVelocity(
          computeDisplacement(current_frame, pos_next, pos_prev) / two_dt);
    }
  }
}

} // namespace

void VelocityIntegrator::calculateVelocities(std::span<Cell> frames, real_t time_step) {
  if (frames.size() < 2 || time_step <= static_cast<real_t>(0.0)) {
    return;
  }

  computeForwardBoundary(frames, time_step);
  computeBackwardBoundary(frames, time_step);
  computeCentralDifferences(frames, time_step);
}

void VelocityIntegrator::calculateVelocities(Trajectory &trajectory) {
  trajectory.calculateVelocities();
}

} // namespace correlation::core
