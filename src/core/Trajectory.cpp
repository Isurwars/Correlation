/**
 * @file Trajectory.cpp
 * @brief Implementation of the trajectory data container.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/Trajectory.hpp"
#include "core/BondCutoffEngine.hpp"
#include "core/MappedFile.hpp"
#include "core/TrajectoryFilter.hpp"
#include "core/VelocityIntegrator.hpp"
#include "math/LinearAlgebra.hpp"
#include "physics/PhysicalData.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace correlation::core {

Trajectory::Trajectory() : time_step_(1.0) {}

Trajectory::Trajectory(std::vector<Cell> frames, real_t time_step)
    : frames_(std::move(frames)), time_step_(time_step) {
  // Validate all frames against the first one
  if (!frames_.empty()) {
    for (size_t i = 1; i < frames_.size(); ++i) {
      validateFrame(frames_[i]);
    }
    removeDuplicatedFrames();
  }
}

Trajectory::Trajectory(std::shared_ptr<MappedFile> mapped_file, std::vector<size_t> frame_offsets,
                       FrameParser parser, real_t time_step)
    : time_step_(time_step), parser_(std::move(parser)) {
  if (mapped_file && frame_offsets.size() > 1) {
    mapped_files_.push_back(std::move(mapped_file));
    frame_slices_.reserve(frame_offsets.size() - 1);
    for (size_t i = 0; i + 1 < frame_offsets.size(); ++i) {
      frame_slices_.push_back({.file_index = 0,
                               .offset = frame_offsets[i],
                               .size = frame_offsets[i + 1] - frame_offsets[i],});
    }
  }
  if (getFrameCount() > 0) {
    precomputeBondCutoffs();
  }
}

Trajectory::Trajectory(std::vector<std::shared_ptr<MappedFile>> mapped_files,
                       std::vector<MappedFrameSlice> frame_slices, FrameParser parser,
                       real_t time_step)
    : time_step_(time_step), mapped_files_(std::move(mapped_files)),
      frame_slices_(std::move(frame_slices)), parser_(std::move(parser)) {
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

std::span<Cell> Trajectory::frameSpan() {
  ensureMaterialized();
  return frames_;
}

std::span<const Cell> Trajectory::frameSpan() const {
  ensureMaterialized();
  return frames_;
}

size_t Trajectory::getFrameCount() const noexcept {
  if (!mapped_files_.empty()) {
    return frame_slices_.size();
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
  if (!mapped_files_.empty()) {
    const auto &slice = frame_slices_[index];
    if (slice.file_index >= mapped_files_.size() || !mapped_files_[slice.file_index]) {
      throw std::runtime_error("Trajectory::getFrame: invalid mapped file index");
    }
    return parser_(mapped_files_[slice.file_index]->data() + slice.offset, slice.size);
  }
  throw std::runtime_error("Trajectory::getFrame: trajectory is empty");
}

const Cell &Trajectory::firstFrame() const {
  if (!frames_.empty()) {
    return frames_[0];
  }
  std::scoped_lock const lock(*init_mutex_);
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
  if (!mapped_files_.empty()) {
    std::scoped_lock const lock(*init_mutex_);
    if (mapped_files_.empty() || !frames_.empty()) {
      return;
    }
    size_t const count = getFrameCount();
    frames_.resize(count);
    for (size_t i = 0; i < count; ++i) {
      const auto &slice = frame_slices_[i];
      frames_[i] = parser_(mapped_files_[slice.file_index]->data() + slice.offset, slice.size);
    }
    first_frame_.reset();
  }
}

const correlation::analysis::BondCutoffRange &Trajectory::getBondCutoffRange(size_t type1,
                                                                             size_t type2) const {
  if (bond_cutoffs_.empty()) {
    precomputeBondCutoffs();
  }
  return BondCutoffEngine::getBondCutoffRange(bond_cutoffs_, type1, type2);
}

real_t Trajectory::getBondCutoffSQ(size_t type1, size_t type2) const {
  return BondCutoffEngine::getBondCutoffSQ(bond_cutoffs_, type1, type2);
}

real_t Trajectory::getMinBondCutoffSQ(size_t type1, size_t type2) const {
  return BondCutoffEngine::getMinBondCutoffSQ(bond_cutoffs_, type1, type2);
}

real_t Trajectory::getBondCutoff(size_t type1, size_t type2) const {
  return BondCutoffEngine::getBondCutoff(bond_cutoffs_, type1, type2);
}

real_t Trajectory::getMinBondCutoff(size_t type1, size_t type2) const {
  return BondCutoffEngine::getMinBondCutoff(bond_cutoffs_, type1, type2);
}

void Trajectory::setBondCutoffs(const correlation::analysis::BondCutoffMatrix &cutoffs) {
  BondCutoffEngine::validate(cutoffs);
  bond_cutoffs_ = cutoffs;
}

void Trajectory::addFrame(const Cell &frame) {
  ensureMaterialized();
  mapped_files_.clear();
  frame_slices_.clear();
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
  bond_cutoffs_ = BondCutoffEngine::compute(getFrame(0));
}

void Trajectory::removeDuplicatedFrames(real_t tolerance) {
  ensureMaterialized();
  mapped_files_.clear();
  frame_slices_.clear();
  parser_ = nullptr;

  if (frames_.size() < 2) {
    return;
  }

  auto unique_frames = TrajectoryFilter::filterConsecutiveDuplicates(frames_, tolerance);
  if (unique_frames.size() < frames_.size()) {
    removed_frames_count_ = frames_.size() - unique_frames.size();
    frames_ = std::move(unique_frames);
  }
}

void Trajectory::calculateVelocities() {
  ensureMaterialized();
  mapped_files_.clear();
  frame_slices_.clear();
  parser_ = nullptr;

  VelocityIntegrator::calculateVelocities(frames_, time_step_);
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
                             std::to_string(ref_elements.size()) + ", but got " +
                             std::to_string(new_elements.size()));
  }

  // Map new element IDs to reference element IDs for fast comparison
  std::vector<int> new_to_ref(new_elements.size(), -1);
  for (size_t i = 0; i < new_elements.size(); ++i) {
    const auto element_iterator = std::ranges::find(ref_elements, new_elements[i]);
    if (element_iterator == ref_elements.end()) {
      throw std::runtime_error("Frame validation failed: Element symbol mismatch at index " +
                               std::to_string(i) + ". Expected " + ref_elements[i].symbol +
                               ", but got " + new_elements[i].symbol);
    }
    new_to_ref[i] = static_cast<int>(std::distance(ref_elements.begin(), element_iterator));
  }

  // Check if the element IDs match in strict order
  // This ensures that the atom types match in the same order as the reference
  // frame
  const auto &ref_atoms = reference.atoms();
  const auto &new_atoms = new_frame.atoms();

  for (size_t i = 0; i < ref_atoms.size(); ++i) {
    const int &ref_id = ref_atoms[i].elementId();
    const int &new_id = new_to_ref[new_atoms[i].elementId()];

    if (ref_id != new_id) {
      // For a helpful error message, we get the mapped original new_id
      const int original_new_id = new_atoms[i].elementId();
      throw std::runtime_error("Frame validation failed: Atom symbol mismatch at index " +
                               std::to_string(i) + ". Expected " + ref_elements[ref_id].symbol +
                               ", but got " + new_elements[original_new_id].symbol);
    }
  }
}
} // namespace correlation::core
