// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "Trajectory.hpp"
#include "PhysicalData.hpp"
#include <cmath>
#include <stdexcept>

//---------------------------------------------------------------------------//
//----------------------------- Constructors --------------------------------//
//---------------------------------------------------------------------------//

Trajectory::Trajectory() : time_step_(1.0) {}

Trajectory::Trajectory(std::vector<Cell> frames, double time_step)
    : frames_(std::move(frames)), time_step_(time_step) {
  // Validate all frames against the first one
  if (!frames_.empty()) {
    for (size_t i = 1; i < frames_.size(); ++i) {
      validateFrame(frames_[i]);
    }
    removeDuplicatedFrames();
  }
}

//---------------------------------------------------------------------------//
//------------------------------- Accessors ---------------------------------//
//---------------------------------------------------------------------------//

double Trajectory::getBondCutoff(int type1, int type2) {
  if (bond_cutoffs_sq_.empty())
    precomputeBondCutoffs();

  if (bond_cutoffs_sq_.empty() || type1 >= bond_cutoffs_sq_.size() ||
      type2 >= bond_cutoffs_sq_.size())
    return 0.0;

  return std::sqrt(bond_cutoffs_sq_[type1][type2]);
}

//---------------------------------------------------------------------------//
//-------------------------------- Methods ----------------------------------//
//---------------------------------------------------------------------------//

void Trajectory::addFrame(const Cell &frame) {
  if (!frames_.empty()) {
    validateFrame(frame);
  }
  frames_.push_back(frame);
}

void Trajectory::precomputeBondCutoffs() {
  if (frames_.empty())
    return;

  const auto &elements = frames_[0].elements();
  const size_t num_elements = elements.size();
  bond_cutoffs_sq_.resize(num_elements, std::vector<double>(num_elements));

  for (size_t i = 0; i < num_elements; ++i) {
    const double radius_A = CovalentRadii::get(elements[i].symbol);
    for (size_t j = i; j < num_elements; ++j) {
      const double radius_B = CovalentRadii::get(elements[j].symbol);
      const double max_bond_dist = (radius_A + radius_B) * 1.3;
      const double max_bond_dist_sq = max_bond_dist * max_bond_dist;
      bond_cutoffs_sq_[i][j] = max_bond_dist_sq;
      bond_cutoffs_sq_[j][i] = max_bond_dist_sq;
    }
  }
}

void Trajectory::removeDuplicatedFrames() {
  if (frames_.size() < 2) {
    return;
  }

  std::vector<Cell> unique_frames;
  unique_frames.reserve(frames_.size());
  unique_frames.push_back(frames_[0]);

  const double epsilon = 1e-5; // Tolerance for position comparison

  for (size_t i = 1; i < frames_.size(); ++i) {
    const auto &current_frame = frames_[i];
    const auto &last_unique_frame = unique_frames.back();

    bool is_duplicate = true;

    // Check atom count
    if (current_frame.atomCount() != last_unique_frame.atomCount()) {
      is_duplicate = false;
    } else {
      // Check positions
      const auto &current_atoms = current_frame.atoms();
      const auto &last_atoms = last_unique_frame.atoms();

      for (size_t j = 0; j < current_atoms.size(); ++j) {
        if (linalg::norm(current_atoms[j].position() -
                         last_atoms[j].position()) > epsilon) {
          is_duplicate = false;
          break;
        }
      }
    }

    if (!is_duplicate) {
      unique_frames.push_back(current_frame);
    }
  }

  if (unique_frames.size() < frames_.size()) {
    removed_frames_count_ = frames_.size() - unique_frames.size();
    frames_ = std::move(unique_frames);
  }
}

void Trajectory::calculateVelocities() {
  if (frames_.size() < 2)
    return;
  size_t num_frames = frames_.size();
  size_t num_atoms = frames_[0].atoms().size();

  velocities_.assign(num_frames,
                     std::vector<linalg::Vector3<double>>(num_atoms));

  if (time_step_ <= 0.0)
    return; // Cannot calculate valid velocities

  for (size_t t = 0; t < num_frames; ++t) {
    // Determine simulation box for PBC (using current frame)
    const auto &lattice = frames_[t].latticeVectors();
    linalg::Vector3<double> box = {lattice[0][0], lattice[1][1], lattice[2][2]};
    // Check if box is valid (not zero), otherwise disable PBC correction
    bool use_pbc = (box[0] > 0.0 && box[1] > 0.0 && box[2] > 0.0);

    // Helper lambda to get minimum image displacement
    auto displacement = [&](const linalg::Vector3<double> &r2,
                            const linalg::Vector3<double> &r1) {
      linalg::Vector3<double> dr = r2 - r1;
      if (use_pbc) {
        if (dr[0] > box[0] * 0.5)
          dr[0] -= box[0];
        if (dr[0] < -box[0] * 0.5)
          dr[0] += box[0];
        if (dr[1] > box[1] * 0.5)
          dr[1] -= box[1];
        if (dr[1] < -box[1] * 0.5)
          dr[1] += box[1];
        if (dr[2] > box[2] * 0.5)
          dr[2] -= box[2];
        if (dr[2] < -box[2] * 0.5)
          dr[2] += box[2];
      }
      return dr;
    };

    for (size_t i = 0; i < num_atoms; ++i) {
      if (t == 0) {
        // Forward difference for the first frame
        // v(0) = (r(1) - r(0)) / dt
        const auto &r0 = frames_[0].atoms()[i].position();
        const auto &r1 = frames_[1].atoms()[i].position();
        velocities_[t][i] = displacement(r1, r0) / time_step_;
      } else if (t == num_frames - 1) {
        // Backward difference for the last frame
        // v(N) = (r(N) - r(N-1)) / dt
        const auto &rN = frames_[num_frames - 1].atoms()[i].position();
        const auto &rN_1 = frames_[num_frames - 2].atoms()[i].position();
        velocities_[t][i] = displacement(rN, rN_1) / time_step_;
      } else {
        // Central difference for internal frames
        // v(t) = (r(t+1) - r(t-1)) / (2 * dt)
        const auto &r_next = frames_[t + 1].atoms()[i].position();
        const auto &r_prev = frames_[t - 1].atoms()[i].position();
        velocities_[t][i] = displacement(r_next, r_prev) / (2.0 * time_step_);
      }
    }
  }
}

//---------------------------------------------------------------------------//
//--------------------------- Private Methods -------------------------------//
//---------------------------------------------------------------------------//

void Trajectory::validateFrame(const Cell &new_frame) const {
  if (frames_.empty()) {
    return;
  }

  const Cell &reference = frames_[0];

  if (new_frame.atomCount() != reference.atomCount()) {
    throw std::runtime_error(
        "Frame validation failed: Atom count mismatch. Expected " +
        std::to_string(reference.atomCount()) + ", but got " +
        std::to_string(new_frame.atomCount()));
  }

  // We should also ideally check if the atom types match in strict order,
  // but for now, we rely on count consistency.
  // Warning: If element registration order differs, elementIDs might differ
  // even if symbols are same. Check symbols.
  const auto &ref_atoms = reference.atoms();
  const auto &new_atoms = new_frame.atoms();
  const auto &ref_elements = reference.elements();
  const auto &new_elements = new_frame.elements();

  for (size_t i = 0; i < ref_atoms.size(); ++i) {
    const std::string &ref_sym = ref_elements[ref_atoms[i].element_id()].symbol;
    const std::string &new_sym = new_elements[new_atoms[i].element_id()].symbol;

    if (ref_sym != new_sym) {
      throw std::runtime_error(
          "Frame validation failed: Atom symbol mismatch at index " +
          std::to_string(i) + ". Expected " + ref_sym + ", but got " + new_sym);
    }
  }
}
