// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/DistanceCalculator.hpp"

#include <cmath>
#include <numeric>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for_each.h>

namespace calculators {

struct ThreadLocalDistances {
  DistanceTensor distance_tensor_local;
  std::vector<std::vector<Neighbor>> neighbor_list_local;

  ThreadLocalDistances(size_t num_elements, size_t atom_count)
      : distance_tensor_local(num_elements,
                              std::vector<std::vector<double>>(num_elements)),
        neighbor_list_local(atom_count) {}
};

void DistanceCalculator::compute(
    const Cell &cell, double cutoff_sq,
    const std::vector<std::vector<double>> &bond_cutoffs_sq,
    bool ignore_periodic_self_interactions, DistanceTensor &out_distances,
    NeighborGraph &out_graph) {

  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const size_t num_elements = cell.elements().size();
  const auto &lattice = cell.latticeVectors();

  linalg::Vector3<double> box_sidelengths = {linalg::norm(lattice[0]),
                                             linalg::norm(lattice[1]),
                                             linalg::norm(lattice[2])};
  int nx =
      static_cast<int>(std::ceil(std::sqrt(cutoff_sq) / box_sidelengths.x()));
  int ny =
      static_cast<int>(std::ceil(std::sqrt(cutoff_sq) / box_sidelengths.y()));
  int nz =
      static_cast<int>(std::ceil(std::sqrt(cutoff_sq) / box_sidelengths.z()));

  if (nx + ny + nz > 8) {
    ignore_periodic_self_interactions = false;
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

  tbb::enumerable_thread_specific<ThreadLocalDistances> ets(num_elements,
                                                            atom_count);

  tbb::parallel_for_each(
      atom_indices.begin(), atom_indices.end(), [&](size_t i) {
        ThreadLocalDistances &local_results = ets.local();
        auto &distance_local = local_results.distance_tensor_local;
        auto &neighbor_local = local_results.neighbor_list_local;

        const auto &atom_A = atoms[i];
        const int type_A = atom_A.element_id();

        for (size_t j = i; j < atom_count; ++j) {
          if (i == j && ignore_periodic_self_interactions) {
            continue;
          }
          const auto &atom_B = atoms[j];
          const int type_B = atom_B.element_id();
          double max_bond_dist_sq = 0.0;
          if (static_cast<size_t>(type_A) < bond_cutoffs_sq.size() &&
              static_cast<size_t>(type_B) < bond_cutoffs_sq[type_A].size()) {
            max_bond_dist_sq = bond_cutoffs_sq[type_A][type_B];
          }
          for (const auto &disp : displacements) {
            if (i == j && linalg::norm_sq(disp) < 1e-9) {
              continue;
            }

            linalg::Vector3<double> r_ij =
                atom_B.position() + disp - atom_A.position();
            double d_sq = linalg::norm_sq(r_ij);

            if (d_sq < cutoff_sq) {
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
        out_distances[i][j].insert(
            out_distances[i][j].end(),
            local_results.distance_tensor_local[i][j].begin(),
            local_results.distance_tensor_local[i][j].end());
      }
    }
    for (size_t i = 0; i < atom_count; ++i) {
      for (const auto &neighbor : local_results.neighbor_list_local[i]) {
        out_graph.addDirectedEdge(i, neighbor.index, neighbor.distance,
                                  neighbor.r_ij);
      }
    }
  }
}

} // namespace calculators
