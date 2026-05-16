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

/**
 * @struct ThreadLocalDistances
 * @brief Thread-local storage for parallel distance calculations.
 * 
 * Contains local tensor results, neighbor lists, and scratch buffers for SoA
 * (Structure of Arrays) SIMD operations to avoid frequent memory allocations.
 */
struct ThreadLocalDistances {
  DistanceTensor distance_tensor_local; ///< Local partial distance tensor.
  std::vector<std::vector<correlation::core::Neighbor>> neighbor_list_local; ///< Local partial neighbor list.

  // Per-thread SoA scratch buffers (reused across iterations to avoid allocs)
  std::vector<double> soa_x;      ///< Scratch for x-coordinates.
  std::vector<double> soa_y;      ///< Scratch for y-coordinates.
  std::vector<double> soa_z;      ///< Scratch for z-coordinates.
  std::vector<double> dsq_scratch; ///< Scratch for squared distance results.
  std::vector<size_t> candidate_j; ///< Scratch for candidate indices.

  /**
   * @brief Constructs thread-local storage.
   * @param num_elements Number of unique elements in the cell.
   * @param atom_count Total number of atoms.
   */
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

  // New Cell-List logic
  math::Vector3<double> a = lattice[0];
  math::Vector3<double> b = lattice[1];
  math::Vector3<double> c = lattice[2];
  
  double vol = cell.volume();
  double width_x = vol / correlation::math::norm(correlation::math::cross(b, c));
  double width_y = vol / correlation::math::norm(correlation::math::cross(a, c));
  double width_z = vol / correlation::math::norm(correlation::math::cross(a, b));
  
  double cutoff = std::sqrt(cutoff_sq);
  
  int Kx = std::max(1, static_cast<int>(std::floor(width_x / cutoff)));
  int Ky = std::max(1, static_cast<int>(std::floor(width_y / cutoff)));
  int Kz = std::max(1, static_cast<int>(std::floor(width_z / cutoff)));
  
  int max_dx = (Kx == 1) ? static_cast<int>(std::ceil(cutoff / width_x)) : 1;
  int max_dy = (Ky == 1) ? static_cast<int>(std::ceil(cutoff / width_y)) : 1;
  int max_dz = (Kz == 1) ? static_cast<int>(std::ceil(cutoff / width_z)) : 1;
  
  if (max_dx + max_dy + max_dz > 8) {
      ignore_periodic_self_interactions = false;
  }

  size_t num_bins = static_cast<size_t>(Kx * Ky * Kz);
  std::vector<std::vector<size_t>> bins(num_bins);
  
  std::vector<int> atom_bin(atom_count);
  std::vector<correlation::math::Vector3<double>> wrapped_positions(atom_count);

  for (size_t i = 0; i < atom_count; ++i) {
    correlation::math::Vector3<double> frac =
        cell.inverseLatticeVectors() * atoms[i].position();
    double fx = frac.x() - std::floor(frac.x());
    double fy = frac.y() - std::floor(frac.y());
    double fz = frac.z() - std::floor(frac.z());

    wrapped_positions[i] = lattice[0] * fx + lattice[1] * fy + lattice[2] * fz;

    int cx = std::clamp(static_cast<int>(std::floor(fx * Kx)), 0, Kx - 1);
    int cy = std::clamp(static_cast<int>(std::floor(fy * Ky)), 0, Ky - 1);
    int cz = std::clamp(static_cast<int>(std::floor(fz * Kz)), 0, Kz - 1);

    int bin_idx = cx * (Ky * Kz) + cy * Kz + cz;
    bins[bin_idx].push_back(i);
    atom_bin[i] = bin_idx;
  }

  tbb::enumerable_thread_specific<ThreadLocalDistances> ets(num_elements,
                                                            atom_count);

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, atom_count),
      [&](const tbb::blocked_range<size_t> &r) {
        ThreadLocalDistances &local_results = ets.local();
        auto &distance_local = local_results.distance_tensor_local;
        auto &neighbor_local = local_results.neighbor_list_local;

        auto &soa_x = local_results.soa_x;
        auto &soa_y = local_results.soa_y;
        auto &soa_z = local_results.soa_z;
        auto &dsq_buf = local_results.dsq_scratch;
        auto &candidate_j = local_results.candidate_j;

        for (size_t i = r.begin(); i != r.end(); ++i) {
          const auto &atom_A = atoms[i];
          const int type_A = atom_A.element_id();
          const double ax = wrapped_positions[i].x();
          const double ay = wrapped_positions[i].y();
          const double az = wrapped_positions[i].z();

          int cx = atom_bin[i] / (Ky * Kz);
          int cy = (atom_bin[i] / Kz) % Ky;
          int cz = atom_bin[i] % Kz;
          
          size_t c_count = 0;

          for (int dx = -max_dx; dx <= max_dx; ++dx) {
            int nx_bin = cx + dx;
            int shift_x = (nx_bin >= 0) ? (nx_bin / Kx) : ((nx_bin - Kx + 1) / Kx);
            int wrap_x = nx_bin - shift_x * Kx;
            
            for (int dy = -max_dy; dy <= max_dy; ++dy) {
                int ny_bin = cy + dy;
                int shift_y = (ny_bin >= 0) ? (ny_bin / Ky) : ((ny_bin - Ky + 1) / Ky);
                int wrap_y = ny_bin - shift_y * Ky;
                
                for (int dz = -max_dz; dz <= max_dz; ++dz) {
                    int nz_bin = cz + dz;
                    int shift_z = (nz_bin >= 0) ? (nz_bin / Kz) : ((nz_bin - Kz + 1) / Kz);
                    int wrap_z = nz_bin - shift_z * Kz;
                    
                    int n_bin_idx = wrap_x * (Ky * Kz) + wrap_y * Kz + wrap_z;
                    
                    correlation::math::Vector3<double> disp = lattice[0] * shift_x + lattice[1] * shift_y + lattice[2] * shift_z;
                    bool zero_disp = (shift_x == 0 && shift_y == 0 && shift_z == 0);
                    
                    for (size_t j : bins[n_bin_idx]) {
                        if (j < i) continue;
                        
                        if (i == j && zero_disp) continue;
                        if (i == j && ignore_periodic_self_interactions) continue;
                        
                        auto pos_B = wrapped_positions[j];
                        double shifted_x = pos_B.x() + disp.x();
                        double shifted_y = pos_B.y() + disp.y();
                        double shifted_z = pos_B.z() + disp.z();
                        
                        if (soa_x.size() <= c_count) {
                            soa_x.push_back(shifted_x);
                            soa_y.push_back(shifted_y);
                            soa_z.push_back(shifted_z);
                            candidate_j.push_back(j);
                            dsq_buf.push_back(0.0);
                        } else {
                            soa_x[c_count] = shifted_x;
                            soa_y[c_count] = shifted_y;
                            soa_z[c_count] = shifted_z;
                            candidate_j[c_count] = j;
                        }
                        c_count++;
                    }
                }
            }
        }

        if (c_count > 0) {
            correlation::math::PositionBlock block{soa_x.data(), soa_y.data(), soa_z.data(), c_count};
            correlation::math::compute_dsq_block(ax, ay, az, block, dsq_buf.data());

            for (size_t k = 0; k < c_count; ++k) {
                double d_sq = dsq_buf[k];
                if (d_sq >= cutoff_sq) continue;
                
                size_t j = candidate_j[k];
                double dist = std::sqrt(d_sq);
                int type_B = atoms[j].element_id();
                
                distance_local[type_A][type_B].push_back(dist);
                if (i != j && type_A != type_B) {
                    distance_local[type_B][type_A].push_back(dist);
                }
                
                double max_bond_dist_sq = 0.0;
                if (static_cast<size_t>(type_A) < bond_cutoffs_sq.size() &&
                    static_cast<size_t>(type_B) < bond_cutoffs_sq[type_A].size()) {
                    max_bond_dist_sq = bond_cutoffs_sq[type_A][type_B];
                }
                
                if (d_sq <= max_bond_dist_sq) {
                    correlation::math::Vector3<double> r_ij = {
                        soa_x[k] - ax, soa_y[k] - ay, soa_z[k] - az
                    };
                    neighbor_local[i].push_back({atoms[j].id(), dist, r_ij});
                    if (i != j) {
                        neighbor_local[j].push_back({atom_A.id(), dist, -1.0 * r_ij});
                    }
                }
            }
        }
        } // close for (size_t i...)
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
