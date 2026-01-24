// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/Trajectory.hpp"
#include <stdexcept>

Trajectory::Trajectory() : time_step_(1.0) {}

Trajectory::Trajectory(std::vector<Cell> frames, double time_step)
    : frames_(std::move(frames)), time_step_(time_step) {
  // Validate all frames against the first one
  if (!frames_.empty()) {
    for (size_t i = 1; i < frames_.size(); ++i) {
      validateFrame(frames_[i]);
    }
  }
}

void Trajectory::addFrame(const Cell &frame) {
  if (!frames_.empty()) {
    validateFrame(frame);
  }
  frames_.push_back(frame);
}

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
