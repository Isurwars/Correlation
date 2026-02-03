// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "Trajectory.hpp"
#include <stdexcept>
#include <cmath>
#include "PhysicalData.hpp"

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

void Trajectory::precomputeBondCutoffs() {
  if (frames_.empty()) return;
  
  const auto& elements = frames_[0].elements();
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

double Trajectory::getBondCutoff(int type1, int type2) {
  if (bond_cutoffs_sq_.empty())
    precomputeBondCutoffs();
    
  if (bond_cutoffs_sq_.empty() || type1 >= bond_cutoffs_sq_.size() || type2 >= bond_cutoffs_sq_.size())
      return 0.0;
      
  return std::sqrt(bond_cutoffs_sq_[type1][type2]);
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
        if (linalg::norm(current_atoms[j].position() - last_atoms[j].position()) > epsilon) {
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
      frames_ = std::move(unique_frames);
  }
}
