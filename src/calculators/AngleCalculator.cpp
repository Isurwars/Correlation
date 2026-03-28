// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/AngleCalculator.hpp"
#include "DistributionFunctions.hpp"

#include <algorithm>
#include <cmath>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include "math/SIMDUtils.hpp"

namespace calculators {

void AngleCalculator::calculateFrame(DistributionFunctions &df,
                                     const AnalysisSettings &settings) const {
  // AngleCalculator is a foundational calculator. It is currently
  // called by StructureAnalyzer during its construction.
}

void AngleCalculator::compute(const Cell &cell, const NeighborGraph &graph,
                              AngleTensor &out_angles) {
  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const size_t num_elements = cell.elements().size();

  // Thread-local AngleTensor accumulator
  tbb::enumerable_thread_specific<AngleTensor> ets([&]() {
    return AngleTensor(
        num_elements,
        std::vector<std::vector<std::vector<double>>>(
            num_elements, std::vector<std::vector<double>>(num_elements)));
  });

  // Thread-local SoA scratch: pre-built per central atom and reused across
  // atoms on the same thread.  Grows to max(CN) and stays there – zero
  // heap allocs after the first high-CN atom is processed.
  struct AngleScratch {
    std::vector<double> nb_x, nb_y, nb_z, nb_dist, dots;
    void ensure(size_t cn) {
      if (nb_x.size() < cn) {
        nb_x.resize(cn);
        nb_y.resize(cn);
        nb_z.resize(cn);
        nb_dist.resize(cn);
        dots.resize(cn);
      }
    }
  };
  tbb::enumerable_thread_specific<AngleScratch> scratch_ets;

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, atom_count),
      [&](const tbb::blocked_range<size_t> &r) {
        auto &local_tensor = ets.local();
        auto &sc = scratch_ets.local();

        for (size_t i = r.begin(); i != r.end(); ++i) {
          const auto &neighbors = graph.getNeighbors(i);
          const size_t cn = neighbors.size();
          if (cn < 2)
            continue;

          const int type_central = atoms[i].element_id();

          // Build SoA for all neighbors of atom i (once per atom)
          sc.ensure(cn);
          for (size_t a = 0; a < cn; ++a) {
            const auto &r_ij = neighbors[a].r_ij;
            sc.nb_x[a] = r_ij.x();
            sc.nb_y[a] = r_ij.y();
            sc.nb_z[a] = r_ij.z();
            sc.nb_dist[a] = neighbors[a].distance;
          }

          // For each j: batch-compute dot(v_j, v_k) for all k > j via SIMD,
          // then apply scalar acos per (j,k) pair.
          for (size_t j = 0; j < cn - 1; ++j) {
            const size_t k_count = cn - j - 1; // number of k's for this j

            // SIMD dot products: dots[m] = dot(v_j, v_{j+1+m})
            correlation::math::simd::dot_block(
                sc.nb_x[j], sc.nb_y[j], sc.nb_z[j], sc.nb_x.data() + j + 1,
                sc.nb_y.data() + j + 1, sc.nb_z.data() + j + 1, sc.dots.data(),
                k_count);

            const int type1 = atoms[neighbors[j].index].element_id();
            const double d1 = sc.nb_dist[j];

            for (size_t m = 0; m < k_count; ++m) {
              const size_t k = j + 1 + m;
              const int type2 = atoms[neighbors[k].index].element_id();

              double cos_theta =
                  std::clamp(sc.dots[m] / (d1 * sc.nb_dist[k]), -1.0, 1.0);
              double angle_rad = std::acos(cos_theta);

              local_tensor[type1][type_central][type2].push_back(angle_rad);
              if (type1 != type2)
                local_tensor[type2][type_central][type1].push_back(angle_rad);
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
