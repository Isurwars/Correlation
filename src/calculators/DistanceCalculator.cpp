/**
 * @file DistanceCalculator.cpp
 * @brief Implementation of pairwise distance calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "calculators/DistanceCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/SIMDUtils.hpp"

#include <cmath>
#include <numeric>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for_each.h>

namespace correlation::calculators {

void DistanceCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &df,
    const correlation::analysis::AnalysisSettings &settings) const {
  // DistanceCalculator is a foundational calculator. In the current
  // architecture, it's called by StructureAnalyzer, which df already has.
  // We provide this implementation for completeness within the BaseCalculator
  // framework.
  // Note: Since StructureAnalyzer already runs this in its constructor,
  // calling it again here might be redundant if df.neighbors() is already
  // populated.
}

struct ThreadLocalDistances {
  DistanceTensor distance_tensor_local;
  std::vector<std::vector<correlation::core::Neighbor>> neighbor_list_local;

  // Per-thread SoA scratch buffers (reused across iterations to avoid allocs)
  std::vector<double> soa_x;
  std::vector<double> soa_y;
  std::vector<double> soa_z;
  std::vector<double> dsq_scratch;

  ThreadLocalDistances(size_t num_elements, size_t atom_count)
      : distance_tensor_local(num_elements,
                              std::vector<std::vector<double>>(num_elements)),
        neighbor_list_local(atom_count) {}
};

void DistanceCalculator::compute(
    const correlation::core::Cell &cell, double cutoff_sq,
    const std::vector<std::vector<double>> &bond_cutoffs_sq,
    bool ignore_periodic_self_interactions, DistanceTensor &out_distances,
    correlation::core::NeighborGraph &out_graph) {

  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const size_t num_elements = cell.elements().size();
  const auto &lattice = cell.latticeVectors();

  correlation::math::Vector3<double> box_sidelengths = {
      correlation::math::norm(lattice[0]), correlation::math::norm(lattice[1]),
      correlation::math::norm(lattice[2])};
  int nx =
      static_cast<int>(std::ceil(std::sqrt(cutoff_sq) / box_sidelengths.x()));
  int ny =
      static_cast<int>(std::ceil(std::sqrt(cutoff_sq) / box_sidelengths.y()));
  int nz =
      static_cast<int>(std::ceil(std::sqrt(cutoff_sq) / box_sidelengths.z()));

  if (nx + ny + nz > 8) {
    ignore_periodic_self_interactions = false;
  }

  std::vector<correlation::math::Vector3<double>> displacements;
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

        auto &soa_x = local_results.soa_x;
        auto &soa_y = local_results.soa_y;
        auto &soa_z = local_results.soa_z;
        auto &dsq_buf = local_results.dsq_scratch;

        const auto &atom_A = atoms[i];
        const int type_A = atom_A.element_id();
        const double ax = atom_A.position().x();
        const double ay = atom_A.position().y();
        const double az = atom_A.position().z();

        // j ranges from i to atom_count-1 (upper triangle + self)
        const size_t j_count = atom_count - i;

        // Ensure scratch buffers are large enough
        if (soa_x.size() < j_count) {
          soa_x.resize(j_count);
          soa_y.resize(j_count);
          soa_z.resize(j_count);
          dsq_buf.resize(j_count);
        }

        // -----------------------------------------------------------------------
        // For each periodic image displacement, build a shifted SoA block for
        // all j in [i, atom_count) and run the SIMD squared-distance kernel.
        // -----------------------------------------------------------------------
        for (const auto &disp : displacements) {
          const double disp_sq = correlation::math::norm_sq(disp);
          const bool zero_disp = (disp_sq < 1e-9);

          // Build SoA: shifted positions of atom_B candidates
          for (size_t jj = 0; jj < j_count; ++jj) {
            const auto &pos_B = atoms[i + jj].position();
            soa_x[jj] = pos_B.x() + disp.x();
            soa_y[jj] = pos_B.y() + disp.y();
            soa_z[jj] = pos_B.z() + disp.z();
          }

          // SIMD pass: compute dsq[jj] = ||atom_A - shifted_atom_B[jj]||²
          correlation::math::PositionBlock block{soa_x.data(), soa_y.data(),
                                                 soa_z.data(), j_count};
          correlation::math::compute_dsq_block(ax, ay, az, block,
                                               dsq_buf.data());

          // -----------------------------------------------------------------------
          // Scalar post-processing: apply cutoff and record output
          // Mirrors the original logic exactly:
          //   - i == j with zero-displacement → skip always (same atom, same
          //   image)
          //   - i == j with non-zero displacement → periodic self-image:
          //       skip if ignore_periodic_self_interactions is true
          //   - i != j → normal pair
          // -----------------------------------------------------------------------
          for (size_t jj = 0; jj < j_count; ++jj) {
            const size_t j = i + jj;

            // Skip self with zero displacement (exact same image of same atom)
            if (i == j && zero_disp) {
              continue;
            }
            // Skip periodic self-images if requested
            if (i == j && ignore_periodic_self_interactions) {
              continue;
            }

            const double d_sq = dsq_buf[jj];
            if (d_sq >= cutoff_sq)
              continue;

            const double dist = std::sqrt(d_sq);

            distance_local[type_A][atoms[j].element_id()].push_back(dist);
            if (i != j && type_A != atoms[j].element_id()) {
              distance_local[atoms[j].element_id()][type_A].push_back(dist);
            }

            // Neighbor graph
            const int type_B = atoms[j].element_id();
            double max_bond_dist_sq = 0.0;
            if (static_cast<size_t>(type_A) < bond_cutoffs_sq.size() &&
                static_cast<size_t>(type_B) < bond_cutoffs_sq[type_A].size()) {
              max_bond_dist_sq = bond_cutoffs_sq[type_A][type_B];
            }

            if (d_sq <= max_bond_dist_sq) {
              // r_ij = (pos_B + disp) - pos_A  — already stored in soa arrays
              correlation::math::Vector3<double> r_ij = {
                  soa_x[jj] - ax, soa_y[jj] - ay, soa_z[jj] - az};
              neighbor_local[i].push_back({atoms[j].id(), dist, r_ij});
              if (i != j) {
                neighbor_local[j].push_back({atom_A.id(), dist, -1.0 * r_ij});
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

} // namespace correlation::calculators
