// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/AngleCalculator.hpp"

#include <algorithm>
#include <cmath>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace calculators {

void AngleCalculator::compute(const Cell &cell, const NeighborGraph &graph,
                              AngleTensor &out_angles) {
  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const size_t num_elements = cell.elements().size();

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
          const auto &central_atom_neighbors = graph.getNeighbors(i);

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
          out_angles[i][j][k].insert(out_angles[i][j][k].end(),
                                     local_tensor[i][j][k].begin(),
                                     local_tensor[i][j][k].end());
        }
      }
    }
  }
}

} // namespace calculators
