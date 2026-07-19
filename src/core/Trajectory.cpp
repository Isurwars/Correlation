/**
 * @file Trajectory.cpp
 * @brief Implementation of the trajectory data container.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/Trajectory.hpp"
#include "core/MappedFile.hpp"
#include "math/LinearAlgebra.hpp"
#include "physics/PhysicalData.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace correlation::core {



Trajectory::Trajectory() : time_step_(1.0) {}

Trajectory::Trajectory(std::vector<Cell> frames, double time_step) : frames_(std::move(frames)), time_step_(time_step) {
  // Validate all frames against the first one
  if (!frames_.empty()) {
    for (size_t i = 1; i < frames_.size(); ++i) {
      validateFrame(frames_[i]);
    }
    removeDuplicatedFrames();
  }
}

Trajectory::Trajectory(std::shared_ptr<MappedFile> mapped_file, std::vector<size_t> frame_offsets, FrameParser parser,
                       double time_step)
    : mapped_file_(std::move(mapped_file)), frame_offsets_(std::move(frame_offsets)), parser_(std::move(parser)),
      time_step_(time_step) {
  if (getFrameCount() > 0) {
    precomputeBondCutoffs();
  }
}



std::vector<Cell> &Trajectory::getFrames() {
  ensureMaterialized();
  return frames_;
}

const std::vector<Cell> &Trajectory::getFrames() const {
  ensureMaterialized();
  return frames_;
}

size_t Trajectory::getFrameCount() const {
  if (mapped_file_) {
    return frame_offsets_.size() - 1;
  }
  return frames_.size();
}

Cell Trajectory::getFrame(size_t index) const {
  if (index >= getFrameCount()) {
    throw std::out_of_range("Trajectory::getFrame: index out of range");
  }
  if (!frames_.empty()) {
    return frames_[index];
  }
  if (mapped_file_) {
    return parser_(mapped_file_->data() + frame_offsets_[index], frame_offsets_[index + 1] - frame_offsets_[index]);
  }
  throw std::runtime_error("Trajectory::getFrame: trajectory is empty");
}

const Cell &Trajectory::firstFrame() const {
  if (!frames_.empty()) {
    return frames_[0];
  }
  if (!first_frame_) {
    if (getFrameCount() == 0) {
      throw std::runtime_error("Trajectory::firstFrame: trajectory is empty");
    }
    first_frame_ = getFrame(0);
  }
  return *first_frame_;
}

void Trajectory::ensureMaterialized() const {
  if (mapped_file_) {
    std::scoped_lock const lock(*init_mutex_);
    if (!mapped_file_ || !frames_.empty()) {
      return;
    }
    size_t const count = getFrameCount();
    frames_.resize(count);
    for (size_t i = 0; i < count; ++i) {
      frames_[i] = parser_(mapped_file_->data() + frame_offsets_[i], frame_offsets_[i + 1] - frame_offsets_[i]);
    }
    first_frame_.reset();
  }
}

double Trajectory::getBondCutoffSQ(size_t type1, size_t type2) const {
  if (bond_cutoffs_sq_.empty()) {
    precomputeBondCutoffs();
  }

  if (type1 >= bond_cutoffs_sq_.size() || type2 >= bond_cutoffs_sq_.size()) {
    return 0.0;
  }

  return bond_cutoffs_sq_[type1][type2];
}

double Trajectory::getBondCutoff(size_t type1, size_t type2) const { return std::sqrt(getBondCutoffSQ(type1, type2)); }



void Trajectory::addFrame(const Cell &frame) {
  ensureMaterialized();
  mapped_file_.reset();
  frame_offsets_.clear();
  parser_ = nullptr;

  if (!frames_.empty()) {
    validateFrame(frame);
  }
  frames_.push_back(frame);
  if (frames_.size() == 1) {
    precomputeBondCutoffs();
  }
}

void Trajectory::precomputeBondCutoffs() const {
  if (getFrameCount() == 0) {
    return;
  }

  Cell const first_frame = getFrame(0);
  const auto &elements = first_frame.elements();
  const size_t num_elements = elements.size();
  bond_cutoffs_sq_.resize(num_elements, std::vector<double>(num_elements));

  auto safeGetRadius = [](const std::string &symbol) -> double {
    try {
      return physics::getCovalentRadius(symbol);
    } catch (const std::out_of_range &) {
      return 1.5; // Default covalent radius for unknown elements
    }
  };

  for (size_t i = 0; i < num_elements; ++i) {
    const double radius_A = safeGetRadius(elements[i].symbol);
    for (size_t j = i; j < num_elements; ++j) {
      const double radius_B = safeGetRadius(elements[j].symbol);
      const double max_bond_dist = (radius_A + radius_B) * 1.3;
      const double max_bond_dist_sq = max_bond_dist * max_bond_dist;
      bond_cutoffs_sq_[i][j] = max_bond_dist_sq;
      bond_cutoffs_sq_[j][i] = max_bond_dist_sq;
    }
  }
}

void Trajectory::removeDuplicatedFrames() {
  ensureMaterialized();
  mapped_file_.reset();
  frame_offsets_.clear();
  parser_ = nullptr;

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
        if (math::norm(current_atoms[j].position() - last_atoms[j].position()) > epsilon) {
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
  ensureMaterialized();
  mapped_file_.reset();
  frame_offsets_.clear();
  parser_ = nullptr;

  if (frames_.size() < 2) {
    return;
  }
  size_t const num_frames = frames_.size();
  size_t const num_atoms = frames_[0].atoms().size();

  if (time_step_ <= 0.0) {
    return; // Cannot calculate valid velocities
  }

  for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
    // Check if the frame has a valid periodic cell (volume > 0).
    const bool use_pbc = (frames_[frame_idx].volume() > 1e-9);

    // Minimum-image displacement: delegates to Cell::minimumImage() which
    // correctly handles both orthogonal and triclinic cells.
    auto displacement = [&](const math::Vector3<real_t> &target_pos,
                            const math::Vector3<real_t> &source_pos) -> math::Vector3<real_t> {
      const math::Vector3<real_t> delta_r = target_pos - source_pos;
      return use_pbc ? frames_[frame_idx].minimumImage(delta_r) : delta_r;
    };
    if (frame_idx == 0) {
      for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
        const auto &pos0 = frames_[0].atoms()[atom_idx].position();
        const auto &pos1 = frames_[1].atoms()[atom_idx].position();
        frames_[frame_idx].atoms()[atom_idx].setVelocity(displacement(pos1, pos0) / time_step_);
      }
    } else if (frame_idx == num_frames - 1) {
      for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
        const auto &posN = frames_[num_frames - 1].atoms()[atom_idx].position();
        const auto &posN_1 = frames_[num_frames - 2].atoms()[atom_idx].position();
        frames_[frame_idx].atoms()[atom_idx].setVelocity(displacement(posN, posN_1) / time_step_);
      }
    } else {
      for (size_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
        // Central difference for internal frames
        // v(t) = (r(t+1) - r(t-1)) / (2 * dt)
        const auto &pos_next = frames_[frame_idx + 1].atoms()[atom_idx].position();
        const auto &pos_prev = frames_[frame_idx - 1].atoms()[atom_idx].position();
        frames_[frame_idx].atoms()[atom_idx].setVelocity(displacement(pos_next, pos_prev) / (2.0 * time_step_));
      }
    }
  }
}



void Trajectory::validateFrame(const Cell &new_frame) const {
  if (getFrameCount() == 0) {
    return;
  }

  Cell reference = getFrame(0);

  if (new_frame.atomCount() != reference.atomCount()) {
    throw std::runtime_error("Frame validation failed: Atom count mismatch. Expected " +
                             std::to_string(reference.atomCount()) + ", but got " +
                             std::to_string(new_frame.atomCount()));
  }

  // Check the element in the cell match the reference frame
  const auto &ref_elements = reference.elements();
  const auto &new_elements = new_frame.elements();
  if (ref_elements.size() != new_elements.size()) {
    throw std::runtime_error("Frame validation failed: Element count mismatch. Expected " +
                             std::to_string(ref_elements.size()) + ", but got " + std::to_string(new_elements.size()));
  }

  // Map new element IDs to reference element IDs for fast comparison
  std::vector<int> new_to_ref(new_elements.size(), -1);
  for (size_t i = 0; i < new_elements.size(); ++i) {
    auto element_iterator = std::ranges::find(ref_elements, new_elements[i]);
    if (element_iterator == ref_elements.end()) {
      throw std::runtime_error("Frame validation failed: Element symbol mismatch at index " + std::to_string(i) +
                               ". Expected " + ref_elements[i].symbol + ", but got " + new_elements[i].symbol);
    }
    new_to_ref[i] = static_cast<int>(std::distance(ref_elements.begin(), element_iterator));
  }

  // Check if the element IDs match in strict order
  // This ensures that the atom types match in the same order as the reference
  // frame
  const auto &ref_atoms = reference.atoms();
  const auto &new_atoms = new_frame.atoms();

  for (size_t i = 0; i < ref_atoms.size(); ++i) {
    const int &ref_id = ref_atoms[i].element_id();
    const int &new_id = new_to_ref[new_atoms[i].element_id()];

    if (ref_id != new_id) {
      // For a helpful error message, we get the mapped original new_id
      const int original_new_id = new_atoms[i].element_id();
      throw std::runtime_error("Frame validation failed: Atom symbol mismatch at index " + std::to_string(i) +
                               ". Expected " + ref_elements[ref_id].symbol + ", but got " +
                               new_elements[original_new_id].symbol);
    }
  }
}
} // namespace correlation::core
