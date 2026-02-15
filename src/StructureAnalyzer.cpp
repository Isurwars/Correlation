// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "StructureAnalyzer.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for_each.h>
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
  neighbor_tensor_.resize(cell.atomCount());

  // The constructor orchestrates the computation
  computeDistances();
  computeAngles();
}

//---------------------------------------------------------------------------//
//--------------------------------- Methods ---------------------------------//
//---------------------------------------------------------------------------//

void StructureAnalyzer::computeDistances() {
  const auto &atoms = cell_.atoms();
  const size_t atom_count = atoms.size();
  const size_t num_elements = cell_.elements().size();
  const auto &lattice = cell_.latticeVectors();

  linalg::Vector3<double> box_sidelengths = {linalg::norm(lattice[0]),
                                             linalg::norm(lattice[1]),
                                             linalg::norm(lattice[2])};
  int nx =
      static_cast<int>(std::ceil(std::sqrt(cutoff_sq_) / box_sidelengths.x()));
  int ny =
      static_cast<int>(std::ceil(std::sqrt(cutoff_sq_) / box_sidelengths.y()));
  int nz =
      static_cast<int>(std::ceil(std::sqrt(cutoff_sq_) / box_sidelengths.z()));

  if (nx + ny + nz > 8) {
    ignore_periodic_self_interactions_ = false;
  }

  std::vector<linalg::Vector3<double>> displacements;
  for (int i = -nx; i <= nx; ++i) {
    for (int j = -ny; j <= ny; ++j) {
      for (int k = -nz; k <= nz; ++k) {
        displacements.push_back(lattice[0] * i + lattice[1] * j +
                                lattice[2] * k);
      }
    }
  }

  std::vector<size_t> atom_indices(atom_count);
  std::iota(atom_indices.begin(), atom_indices.end(), 0);

  // lock-free calculation phase using TBB.
  // We use enumerable_thread_specific to give each thread its own local results
  // storage, avoiding mutex contention.
  tbb::enumerable_thread_specific<ThreadLocalResults> ets(num_elements,
                                                          atom_count);

  tbb::parallel_for_each(
      atom_indices.begin(), atom_indices.end(), [&](size_t i) {
        // Each thread gets its own private copy of the results structure.
        ThreadLocalResults &local_results = ets.local();
        auto &distance_local = local_results.distance_tensor_local;
        auto &neighbor_local = local_results.neighbor_tensor_local;

        const auto &atom_A = atoms[i];
        const int type_A = atom_A.element_id();

        for (size_t j = i; j < atom_count; ++j) {
          if (i == j && ignore_periodic_self_interactions_) {
            continue;
          }
          const auto &atom_B = atoms[j];
          const int type_B = atom_B.element_id();
          double max_bond_dist_sq = 0.0;
          if (type_A < bond_cutoffs_sq_.size() &&
              type_B < bond_cutoffs_sq_[type_A].size()) {
            max_bond_dist_sq = bond_cutoffs_sq_[type_A][type_B];
          }
          for (const auto &disp : displacements) {
            if (i == j && linalg::norm_sq(disp) < 1e-9) {
              continue;
            }

            linalg::Vector3<double> r_ij =
                atom_B.position() + disp - atom_A.position();
            double d_sq = linalg::norm_sq(r_ij);

            if (d_sq < cutoff_sq_) {
              double dist = std::sqrt(d_sq);
              distance_local[type_A][type_B].push_back(dist);
              if (i != j && type_A != type_B) {
                distance_local[type_B][type_A].push_back(dist);
              }

              if (d_sq <= max_bond_dist_sq) {
                neighbor_local[i].push_back({atom_B.id(), dist, r_ij});
                if (i != j) {
                  neighbor_local[j].push_back({atom_A.id(), dist, -1.0 * r_ij});
                }
              }
            }
          }
        }
      });

  for (auto &local_results : ets) {
    for (size_t i = 0; i < num_elements; ++i) {
      for (size_t j = 0; j < num_elements; ++j) {
        distance_tensor_[i][j].insert(
            distance_tensor_[i][j].end(),
            local_results.distance_tensor_local[i][j].begin(),
            local_results.distance_tensor_local[i][j].end());
      }
    }
    for (size_t i = 0; i < atom_count; ++i) {
      neighbor_tensor_[i].insert(neighbor_tensor_[i].end(),
                                 local_results.neighbor_tensor_local[i].begin(),
                                 local_results.neighbor_tensor_local[i].end());
    }
  }
}

void StructureAnalyzer::computeAngles() {
  const auto &atoms = cell_.atoms();
  const size_t atom_count = atoms.size();

  const size_t num_elements = cell_.elements().size();
  using AngleTensor =
      std::vector<std::vector<std::vector<std::vector<double>>>>;

  // Initialize thread-local storage
  tbb::enumerable_thread_specific<AngleTensor> ets([&]() {
    return AngleTensor(
        num_elements,
        std::vector<std::vector<std::vector<double>>>(
            num_elements, std::vector<std::vector<double>>(num_elements)));
  });

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, atom_count),
      [&](const tbb::blocked_range<size_t> &r) {
        auto &local_tensor = ets.local();
        for (size_t i = r.begin(); i != r.end(); ++i) {
          const auto &central_atom_neighbors = neighbor_tensor_[i];

          // Skip if the central atom has fewer than two neighbors
          if (central_atom_neighbors.size() < 2)
            continue;

          const int type_central = atoms[i].element_id();

          // Loop over all unique pairs of neighbors (j, k)
          for (size_t j = 0; j < central_atom_neighbors.size(); ++j) {
            for (size_t k = j + 1; k < central_atom_neighbors.size(); ++k) {
              const auto &neighbor1 = central_atom_neighbors[j];
              const auto &neighbor2 = central_atom_neighbors[k];

              const int type1 = atoms[neighbor1.index].element_id();
              const int type2 = atoms[neighbor2.index].element_id();

              const auto &v1 = neighbor1.r_ij;
              const auto &v2 = neighbor2.r_ij;

              // Calculate angle
              double cos_theta = linalg::dot(v1, v2) /
                                 (neighbor1.distance * neighbor2.distance);
              cos_theta = std::clamp(cos_theta, -1.0, 1.0);
              double angle_rad = std::acos(cos_theta);
              local_tensor[type1][type_central][type2].push_back(angle_rad);
              if (type1 != type2) {
                local_tensor[type2][type_central][type1].push_back(angle_rad);
              }
            }
          }
        }
      });

  // Merge results
  for (const auto &local_tensor : ets) {
    for (size_t i = 0; i < num_elements; ++i) {
      for (size_t j = 0; j < num_elements; ++j) {
        for (size_t k = 0; k < num_elements; ++k) {
          angle_tensor_[i][j][k].insert(angle_tensor_[i][j][k].end(),
                                        local_tensor[i][j][k].begin(),
                                        local_tensor[i][j][k].end());
        }
      }
    }
  }
}
