// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/StructureAnalyzer.hpp"

#include <algorithm>
#include <cmath>
#include <execution>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "../include/PhysicalData.hpp"

//---------------------------------------------------------------------------//
//------------------------------- Constructors ------------------------------//
//---------------------------------------------------------------------------//
StructureAnalyzer::StructureAnalyzer(const Cell &cell, double cutoff,
                                     double bond_factor)
    : cell_(cell), cutoff_sq_(cutoff * cutoff), bond_factor_(bond_factor) {
  if (cutoff <= 0) {
    throw std::invalid_argument("Cutoff distance must be positive.");
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
          num_elements, std::vector<std::vector<double>>(num_elements)));
  neighbor_tensor_.resize(cell.atomCount());

  // The constructor orchestrates the computation
  precomputeBondCutoffs();
  computeDistances();
  computeAngles();
}

//---------------------------------------------------------------------------//
//--------------------------------- Methods ---------------------------------//
//---------------------------------------------------------------------------//

void StructureAnalyzer::precomputeBondCutoffs() {
  const auto &elements = cell_.elements();
  const size_t num_elements = elements.size();
  auto placeholder = std::vector<double>(num_elements);
  bond_cutoffs_sq_.resize(num_elements, placeholder);

  for (size_t i = 0; i < num_elements; ++i) {
    const double radius_A = CovalentRadii::get(elements[i].symbol);
    for (size_t j = i; j < num_elements; ++j) {
      const double radius_B = CovalentRadii::get(elements[j].symbol);
      const double max_bond_dist = (radius_A + radius_B) * bond_factor_;
      const double max_bond_dist_sq = max_bond_dist * max_bond_dist;
      bond_cutoffs_sq_[i][j] = max_bond_dist_sq;
      bond_cutoffs_sq_[j][i] = max_bond_dist_sq;
    }
  }
}

void StructureAnalyzer::computeDistances() {
  const std::vector<Atom> &atoms = cell_.atoms();
  const size_t atom_count = atoms.size();
  const size_t num_elements = cell_.elements().size();
  const linalg::Matrix3<double> &lattice = cell_.latticeVectors();

  // Determine the number of image cells to check based on the cutoff
  linalg::Vector3<double> box_sidelengths = {linalg::norm(lattice[0]),
                                             linalg::norm(lattice[1]),
                                             linalg::norm(lattice[2])};
  int nx =
      static_cast<int>(std::ceil(std::sqrt(cutoff_sq_) / box_sidelengths.x()));
  int ny =
      static_cast<int>(std::ceil(std::sqrt(cutoff_sq_) / box_sidelengths.y()));
  int nz =
      static_cast<int>(std::ceil(std::sqrt(cutoff_sq_) / box_sidelengths.z()));

  // Pre-calculate displacement vectors for periodic images
  std::vector<linalg::Vector3<double>> displacements;
  for (int i = -nx; i <= nx; ++i) {
    for (int j = -ny; j <= ny; ++j) {
      for (int k = -nz; k <= nz; ++k) {
        displacements.push_back(lattice[0] * i + lattice[1] * j +
                                lattice[2] * k);
      }
    }
  }
  // A mutex is used to safely merge results from multiple threads.
  // This is the direct replacement for the '#pragma omp critical' directive.
  std::mutex results_mutex;

  // Create a vector of indices [0, 1, 2, ..., atom_count-1] to iterate over.
  // This is a common pattern for using standard algorithms on an index-based
  // loop.
  std::vector<size_t> atom_indices(atom_count);
  std::iota(atom_indices.begin(), atom_indices.end(), 0);

  // Use std::for_each with a parallel execution policy. The loop body is a
  // lambda.
  std::for_each(
      std::execution::par, atom_indices.begin(), atom_indices.end(),
      [&](size_t i) {
        // These containers are created within the lambda, making them naturally
        // "thread-local" for each task executed by the parallel algorithm.
        DistanceTensor distance_tensor_local(
            num_elements, std::vector<std::vector<double>>(num_elements));
        NeighborTensor neighbor_tensor_local(atom_count);

        const Atom &atom_A = atoms[i];
        const int type_A = atom_A.element_id();

        // Iterate over unique pairs (j>=i)
        for (size_t j = i; j < atom_count; ++j) {
          const Atom &atom_B = atoms[j];
          const int type_B = atom_B.element_id();
          const double max_bond_dist_sq = bond_cutoffs_sq_[type_A][type_B];

          // Check against all periodic images of atom B
          for (const auto &disp : displacements) {
            linalg::Vector3<double> r_ij =
                atom_B.position() + disp - atom_A.position();
            double d_sq = linalg::norm_sq(r_ij);

            if (d_sq > 1e-9 && d_sq < cutoff_sq_) {
              double dist = std::sqrt(d_sq);

              distance_tensor_local[type_A][type_B].push_back(dist);
              if (i != j && type_A != type_B) {
                distance_tensor_local[type_B][type_A].push_back(dist);
              } else if (i == j && type_A != type_B) {
                distance_tensor_local[type_B][type_A].push_back(dist);
              }

              if (d_sq <= max_bond_dist_sq) {
                if (i != j || linalg::norm_sq(disp) > 1e-9) {
                  neighbor_tensor_local[i].push_back({atom_B.id(), dist, r_ij});
                  if (i != j) {
                    neighbor_tensor_local[j].push_back(
                        {atom_A.id(), dist, -1.0 * r_ij});
                  }
                }
              }
            }
          }
        }

        // --- Reduction Step ---
        // Lock the mutex to ensure only one thread at a time modifies the
        // shared member variables (distance_tensor_ and neighbor_tensor_).
        // std::lock_guard automatically unlocks the mutex when it goes out of
        // scope.
        std::lock_guard<std::mutex> lock(results_mutex);
        for (size_t type1 = 0; type1 < num_elements; ++type1) {
          for (size_t type2 = 0; type2 < num_elements; ++type2) {
            if (!distance_tensor_local[type1][type2].empty()) {
              distance_tensor_[type1][type2].insert(
                  distance_tensor_[type1][type2].end(),
                  distance_tensor_local[type1][type2].begin(),
                  distance_tensor_local[type1][type2].end());
            }
          }
        }
        for (size_t atom_idx = 0; atom_idx < atom_count; ++atom_idx) {
          if (!neighbor_tensor_local[atom_idx].empty()) {
            neighbor_tensor_[atom_idx].insert(
                neighbor_tensor_[atom_idx].end(),
                neighbor_tensor_local[atom_idx].begin(),
                neighbor_tensor_local[atom_idx].end());
          }
        }
      }); // End of parallel std::for_each
}
void StructureAnalyzer::computeAngles() {
  const auto &all_neighbors = neighbor_tensor_;
  const auto &atoms = cell_.atoms();
  const size_t atom_count = atoms.size();
  const size_t num_elements = cell_.elements().size();

#pragma omp parallel
  {
    // Thread-local storage for angles
    AngleTensor angle_tensor_local(
        num_elements,
        std::vector<std::vector<std::vector<double>>>(
            num_elements, std::vector<std::vector<double>>(num_elements)));

#pragma omp for schedule(dynamic)
    for (size_t i = 0; i < atom_count; ++i) {
      const auto &central_atom_neighbors = all_neighbors[i];
      if (central_atom_neighbors.size() < 2)
        continue;

      const int type_central = atoms[i].element_id();

      // Iterate over all unique pairs of neighbors for the central atom
      for (size_t j = 0; j < central_atom_neighbors.size(); ++j) {
        for (size_t k = j + 1; k < central_atom_neighbors.size(); ++k) {
          const auto &neighbor1 = central_atom_neighbors[j];
          const auto &neighbor2 = central_atom_neighbors[k];

          const int type1 = atoms[neighbor1.index.value].element_id();
          const int type2 = atoms[neighbor2.index.value].element_id();

          const auto &v1 = neighbor1.r_ij;
          const auto &v2 = neighbor2.r_ij;

          double cos_theta =
              linalg::dot(v1, v2) / (neighbor1.distance * neighbor2.distance);
          cos_theta = std::clamp(cos_theta, -1.0, 1.0);
          double angle_rad = std::acos(cos_theta);
          // Add to the local angle tensor
          angle_tensor_local[type1][type_central][type2].push_back(angle_rad);
          if (type1 != type2) {
            angle_tensor_local[type2][type_central][type1].push_back(angle_rad);
          }
        }
      }
    }

// Merge thread-local angle tensors
#pragma omp critical
    {
      for (size_t i = 0; i < num_elements; ++i) {
        for (size_t j = 0; j < num_elements; ++j) {
          for (size_t k = 0; k < num_elements; ++k) {
            angle_tensor_[i][j][k].insert(angle_tensor_[i][j][k].end(),
                                          angle_tensor_local[i][j][k].begin(),
                                          angle_tensor_local[i][j][k].end());
          }
        }
      }
    }
  }
}
