/**
 * @file StructureAnalyzer.cpp
 * @brief Implementation of structural analysis and neighbor search.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/AngleCalculator.hpp"
#include "calculators/DihedralCalculator.hpp"
#include "calculators/DistanceCalculator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_group.h>
#include <vector>

namespace correlation::analysis {

StructureAnalyzer::StructureAnalyzer(const correlation::core::Cell &cell, real_t cutoff,
                                     const BondCutoffMatrix &bond_cutoffs,
                                     bool ignore_periodic_self_interactions)
    // Use the member initializer list for all members for correctness and
    // efficiency.
    : cell_(cell), cutoff_sq_(cutoff * cutoff), bond_cutoffs_(bond_cutoffs),
      ignore_periodic_self_interactions_(ignore_periodic_self_interactions) {
  if (cutoff <= 0) {
    throw std::invalid_argument("Cutoff distance must be positive.");
  }

  // Validate bond cutoff ranges
  for (size_t i = 0; i < bond_cutoffs_.size(); ++i) {
    for (size_t j = 0; j < bond_cutoffs_[i].size(); ++j) {
      const auto &range = bond_cutoffs_[i][j];
      if (range.min_sq < 0.0) {
        throw std::invalid_argument("Minimum bond cutoff squared must be non-negative.");
      }
      if (range.max_sq > 0.0 && range.max_sq < range.min_sq) {
        throw std::invalid_argument("Maximum bond cutoff squared cannot be less than minimum bond cutoff squared.");
      }
    }
  }

  // Ensure cutoff covers the largest bond distance
  const auto &elements = cell.elements();
  real_t max_bond_dist = 0.0;
  max_bond_dist = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, elements.size()), static_cast<real_t>(0.0),
      [&](const tbb::blocked_range<size_t> &range, real_t init) {
        for (size_t i = range.begin(); i != range.end(); ++i) {
          for (size_t j = i; j < elements.size(); ++j) {
            // Find element indices in the cutoff matrix
            // Assuming bond_cutoffs indices match element indices in frame
            // This assumption holds if trajectory validation works.
            if (i < bond_cutoffs_.size() && j < bond_cutoffs_[i].size()) {
              init = std::max(init, static_cast<real_t>(std::sqrt(bond_cutoffs_[i][j].max_sq)));
            }
          }
        }
        return init;
      },
      [](real_t lhs, real_t rhs) { return std::max(lhs, rhs); });
  if (cutoff < max_bond_dist) {
    cutoff = max_bond_dist;
    cutoff_sq_ = cutoff * cutoff;
  }

  if (cell.isEmpty()) {
    return; // Nothing to compute for an empty cell
  }

  // Initialize the tensors and bond list with the correct dimensions
  const size_t num_elements = cell.elements().size();
  distance_tensor_.resize(num_elements, std::vector<std::vector<real_t>>(num_elements));
  angle_tensor_.resize(num_elements,
                       std::vector<std::vector<std::vector<real_t>>>(
                           num_elements, std::vector<std::vector<real_t>>(num_elements, std::vector<real_t>())));

  dihedral_tensor_.resize(
      num_elements,
      std::vector<std::vector<std::vector<std::vector<real_t>>>>(
          num_elements, std::vector<std::vector<std::vector<real_t>>>(
                            num_elements, std::vector<std::vector<real_t>>(num_elements, std::vector<real_t>()))));

  neighbor_graph_ = correlation::core::NeighborGraph(cell.atomCount());

  // The constructor orchestrates the computation by delegating to dedicated
  // calculators
  correlation::calculators::DistanceCalculator::compute(cell_, cutoff_sq_, bond_cutoffs_,
                                                        ignore_periodic_self_interactions_, distance_tensor_,
                                                        neighbor_graph_);

  tbb::task_group task_group;
  task_group.run([&]() { correlation::calculators::AngleCalculator::compute(cell_, neighbor_graph_, angle_tensor_); });
  task_group.run(
      [&]() { correlation::calculators::DihedralCalculator::compute(cell_, neighbor_graph_, dihedral_tensor_); });
  task_group.wait();
}

} // namespace correlation::analysis
