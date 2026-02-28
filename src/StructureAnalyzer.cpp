// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
#include "StructureAnalyzer.hpp"
#include "calculators/AngleCalculator.hpp"
#include "calculators/DihedralCalculator.hpp"
#include "calculators/DistanceCalculator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

//---------------------------------------------------------------------------//
//------------------------------- Constructors ------------------------------//
//---------------------------------------------------------------------------//
StructureAnalyzer::StructureAnalyzer(
    Cell &cell, double cutoff,
    const std::vector<std::vector<double>> &bond_cutoffs_sq,
    bool ignore_periodic_self_interactions)
    // Use the member initializer list for all members for correctness and
    // efficiency.
    : cell_(cell), cutoff_sq_(cutoff * cutoff),
      bond_cutoffs_sq_(bond_cutoffs_sq),
      ignore_periodic_self_interactions_(ignore_periodic_self_interactions) {
  if (cutoff <= 0) {
    throw std::invalid_argument("Cutoff distance must be positive.");
  }

  // Ensure cutoff covers the largest bond distance
  const auto &elements = cell.elements();
  double max_bond_dist = 0.0;
  for (size_t i = 0; i < elements.size(); ++i) {
    for (size_t j = i; j < elements.size(); ++j) {
      // Find element indices in the cutoff matrix
      // Assuming bond_cutoffs_sq indices match element indices in frame
      // This assumption holds if trajectory validation works.
      if (i < bond_cutoffs_sq.size() && j < bond_cutoffs_sq[i].size()) {
        max_bond_dist =
            std::max(max_bond_dist, std::sqrt(bond_cutoffs_sq[i][j]));
      }
    }
  }
  if (cutoff < max_bond_dist) {
    cutoff = max_bond_dist;
    cutoff_sq_ = cutoff * cutoff;
  }

  if (cell.isEmpty()) {
    return; // Nothing to compute for an empty cell
  }

  // Initialize the tensors and bond list with the correct dimensions
  const size_t num_elements = cell.elements().size();
  distance_tensor_.resize(num_elements,
                          std::vector<std::vector<double>>(num_elements));
  angle_tensor_.resize(
      num_elements,
      std::vector<std::vector<std::vector<double>>>(
          num_elements, std::vector<std::vector<double>>(
                            num_elements, std::vector<double>())));

  dihedral_tensor_.resize(
      num_elements,
      std::vector<std::vector<std::vector<std::vector<double>>>>(
          num_elements,
          std::vector<std::vector<std::vector<double>>>(
              num_elements, std::vector<std::vector<double>>(
                                num_elements, std::vector<double>()))));

  neighbor_graph_ = NeighborGraph(cell.atomCount());

  // The constructor orchestrates the computation by delegating to dedicated
  // calculators
  calculators::DistanceCalculator::compute(cell_, cutoff_sq_, bond_cutoffs_sq_,
                                           ignore_periodic_self_interactions_,
                                           distance_tensor_, neighbor_graph_);
  calculators::AngleCalculator::compute(cell_, neighbor_graph_, angle_tensor_);
  calculators::DihedralCalculator::compute(cell_, neighbor_graph_,
                                           dihedral_tensor_);
}
