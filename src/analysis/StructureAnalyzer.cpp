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
#include <mutex>
#include <stdexcept>
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <utility>
#include <vector>

namespace correlation::analysis {

namespace {

void validateBondCutoffs(const BondCutoffMatrix &bond_cutoffs) {
  for (const auto &row : bond_cutoffs) {
    for (const auto &range : row) {
      if (range.min_sq < 0.0) {
        throw std::invalid_argument("Minimum bond cutoff squared must be non-negative.");
      }
      if (range.max_sq > 0.0 && range.max_sq < range.min_sq) {
        throw std::invalid_argument("Maximum bond cutoff squared cannot be less than minimum bond cutoff squared.");
      }
    }
  }
}

real_t computeMaxBondDistance(const BondCutoffMatrix &bond_cutoffs, size_t num_elements) {
  return tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, num_elements), static_cast<real_t>(0.0),
      [&](const tbb::blocked_range<size_t> &range, real_t init) {
        for (size_t i = range.begin(); i != range.end(); ++i) {
          for (size_t j = i; j < num_elements; ++j) {
            if (i < bond_cutoffs.size() && j < bond_cutoffs[i].size()) {
              init = std::max(init, static_cast<real_t>(std::sqrt(bond_cutoffs[i][j].max_sq)));
            }
          }
        }
        return init;
      },
      [](real_t lhs, real_t rhs) { return std::max(lhs, rhs); });
}

} // namespace

StructureAnalyzer::StructureAnalyzer(std::shared_ptr<const correlation::core::Cell> cell, real_t cutoff,
                                     BondCutoffMatrix bond_cutoffs, bool ignore_periodic_self_interactions,
                                     real_t r_bin_width)
    : cell_(std::move(cell)), cutoff_sq_(cutoff * cutoff), bond_cutoffs_(std::move(bond_cutoffs)),
      ignore_periodic_self_interactions_(ignore_periodic_self_interactions) {
  if (!cell_) {
    throw std::invalid_argument("Cell pointer cannot be null.");
  }
  if (cutoff <= 0) {
    throw std::invalid_argument("Cutoff distance must be positive.");
  }

  validateBondCutoffs(bond_cutoffs_);

  const real_t max_bond_dist = computeMaxBondDistance(bond_cutoffs_, cell_->elements().size());
  if (cutoff < max_bond_dist) {
    cutoff = max_bond_dist;
    cutoff_sq_ = cutoff * cutoff;
  }

  if (cell_->isEmpty()) {
    return; // Nothing to compute for an empty cell
  }

  // Initialize the raw histogram tensor and neighbor graph
  const size_t num_elements = cell_->elements().size();
  const size_t num_bins =
      (r_bin_width > 0.0 && cutoff > 0.0) ? static_cast<size_t>(std::ceil(cutoff / r_bin_width)) : 0;
  raw_histograms_.resize(num_elements,
                         std::vector<std::vector<real_t>>(num_elements, std::vector<real_t>(num_bins, 0.0)));
  neighbor_graph_ = correlation::core::NeighborGraph(cell_->atomCount());

  calculators::DistanceCalculationConfig const hist_config{
      .r_max = cutoff,
      .r_bin_width = r_bin_width,
      .num_bins = num_bins,
  };

  // Delegate distance and neighbor graph computation to DistanceCalculator
  correlation::calculators::DistanceCalculator::compute(*cell_, cutoff_sq_, bond_cutoffs_,
                                                        ignore_periodic_self_interactions_, neighbor_graph_,
                                                        &raw_histograms_, hist_config);
}

StructureAnalyzer::StructureAnalyzer(const correlation::core::Cell &cell, real_t cutoff, BondCutoffMatrix bond_cutoffs,
                                     bool ignore_periodic_self_interactions, real_t r_bin_width)
    : StructureAnalyzer(std::make_shared<correlation::core::Cell>(cell), cutoff, std::move(bond_cutoffs),
                        ignore_periodic_self_interactions, r_bin_width) {}

StructureAnalyzer::StructureAnalyzer(StructureAnalyzer &&other) noexcept
    : cell_(std::move(other.cell_)), cutoff_sq_(other.cutoff_sq_), bond_cutoffs_(std::move(other.bond_cutoffs_)),
      ignore_periodic_self_interactions_(other.ignore_periodic_self_interactions_),
      neighbor_graph_(std::move(other.neighbor_graph_)), raw_histograms_(std::move(other.raw_histograms_)),
      angle_tensor_(std::move(other.angle_tensor_)), dihedral_tensor_(std::move(other.dihedral_tensor_)),
      angles_computed_(other.angles_computed_.load(std::memory_order_relaxed)),
      dihedrals_computed_(other.dihedrals_computed_.load(std::memory_order_relaxed)) {}

StructureAnalyzer &StructureAnalyzer::operator=(StructureAnalyzer &&other) noexcept {
  if (this != &other) {
    const std::scoped_lock lock(compute_mutex_, other.compute_mutex_);
    cell_ = std::move(other.cell_);
    cutoff_sq_ = other.cutoff_sq_;
    bond_cutoffs_ = std::move(other.bond_cutoffs_);
    ignore_periodic_self_interactions_ = other.ignore_periodic_self_interactions_;
    neighbor_graph_ = std::move(other.neighbor_graph_);
    raw_histograms_ = std::move(other.raw_histograms_);
    angle_tensor_ = std::move(other.angle_tensor_);
    dihedral_tensor_ = std::move(other.dihedral_tensor_);
    angles_computed_.store(other.angles_computed_.load(std::memory_order_relaxed), std::memory_order_relaxed);
    dihedrals_computed_.store(other.dihedrals_computed_.load(std::memory_order_relaxed), std::memory_order_relaxed);
  }
  return *this;
}

const StructureAnalyzer::AngleTensor &StructureAnalyzer::angles() const {
  ensureAnglesComputed();
  return angle_tensor_;
}

const StructureAnalyzer::DihedralTensor &StructureAnalyzer::dihedrals() const {
  ensureDihedralsComputed();
  return dihedral_tensor_;
}

void StructureAnalyzer::ensureAnglesComputed() const {
  if (angles_computed_.load(std::memory_order_acquire)) {
    return;
  }
  const std::scoped_lock lock(compute_mutex_);
  if (angles_computed_.load(std::memory_order_relaxed)) {
    return;
  }
  const size_t num_elements = cell_->elements().size();
  angle_tensor_.resize(num_elements, std::vector<std::vector<std::vector<real_t>>>(
                                         num_elements, std::vector<std::vector<real_t>>(num_elements)));

  correlation::calculators::AngleCalculator::compute(*cell_, neighbor_graph_, angle_tensor_);
  angles_computed_.store(true, std::memory_order_release);
}

void StructureAnalyzer::ensureDihedralsComputed() const {
  if (dihedrals_computed_.load(std::memory_order_acquire)) {
    return;
  }
  const std::scoped_lock lock(compute_mutex_);
  if (dihedrals_computed_.load(std::memory_order_relaxed)) {
    return;
  }
  const size_t num_elements = cell_->elements().size();
  dihedral_tensor_.resize(num_elements,
                          std::vector<std::vector<std::vector<std::vector<real_t>>>>(
                              num_elements, std::vector<std::vector<std::vector<real_t>>>(
                                                num_elements, std::vector<std::vector<real_t>>(num_elements))));

  correlation::calculators::DihedralCalculator::compute(*cell_, neighbor_graph_, dihedral_tensor_);
  dihedrals_computed_.store(true, std::memory_order_release);
}

} // namespace correlation::analysis
