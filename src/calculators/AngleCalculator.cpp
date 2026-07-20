/**
 * @file AngleCalculator.cpp
 * @brief Implementation of angular distribution calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/AngleCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "math/Precision.hpp"
#include "math/SIMDUtils.hpp"

#include <algorithm>
#include <cmath>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <utility>

namespace correlation::calculators {

void AngleCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                     const correlation::analysis::AnalysisSettings &settings) const {
  // AngleCalculator is a foundational calculator. It is currently
  // called by StructureAnalyzer during its construction.
}

namespace {

struct AngleScratch {
  std::vector<double> nb_x, nb_y, nb_z, nb_dist, dots;
  void ensure(size_t neighbor_count) {
    if (nb_x.size() < neighbor_count) {
      nb_x.resize(neighbor_count);
      nb_y.resize(neighbor_count);
      nb_z.resize(neighbor_count);
      nb_dist.resize(neighbor_count);
      dots.resize(neighbor_count);
    }
  }
};

void computeTriadAngles(int type_central, const std::vector<correlation::core::Atom> &atoms, size_t neighbor_count,
                        const std::vector<correlation::core::Neighbor> &neighbors, size_t num_elements,
                        AngleScratch &scratch, AngleTensor &local_tensor) {
  for (size_t j_idx = 0; j_idx < neighbor_count - 1; ++j_idx) {
    const size_t k_count = neighbor_count - j_idx - 1; // number of k's for this j

    // SIMD dot products: dots[m_idx] = dot(v_j, v_{j+1+m_idx})
    correlation::math::dot_block(scratch.nb_x[j_idx], scratch.nb_y[j_idx], scratch.nb_z[j_idx],
                                 scratch.nb_x.data() + j_idx + 1, scratch.nb_y.data() + j_idx + 1,
                                 scratch.nb_z.data() + j_idx + 1, scratch.dots.data(), k_count);

    const int type1 = atoms[neighbors[j_idx].index].element_id();
    const double dist1 = scratch.nb_dist[j_idx];
    if (dist1 < 1e-6) {
      continue;
    }
    if (type1 < 0 || std::cmp_greater_equal(type1, num_elements)) {
      continue;
    }

    for (size_t m_idx = 0; m_idx < k_count; ++m_idx) {
      const size_t k_idx = j_idx + 1 + m_idx;
      const double dist2 = scratch.nb_dist[k_idx];
      if (dist2 < 1e-6) {
        continue;
      }

      const int type2 = atoms[neighbors[k_idx].index].element_id();
      if (type2 < 0 || std::cmp_greater_equal(type2, num_elements)) {
        continue;
      }

      const real_t cos_theta = static_cast<real_t>(std::clamp(scratch.dots[m_idx] / (dist1 * dist2), -1.0, 1.0));
      const real_t angle_rad = std::acos(cos_theta);

      local_tensor[type1][type_central][type2].push_back(angle_rad);
      if (type1 != type2) {
        local_tensor[type2][type_central][type1].push_back(angle_rad);
      }
    }
  }
}

void processCentralAtom(size_t atom_idx, const correlation::core::Cell &cell,
                        const correlation::core::NeighborGraph &graph, size_t num_elements, AngleScratch &scratch,
                        AngleTensor &local_tensor) {
  const auto &atoms = cell.atoms();
  const auto &neighbors = graph.getNeighbors(atom_idx);
  const size_t neighbor_count = neighbors.size();
  if (neighbor_count < 2) {
    return;
  }

  const int type_central = atoms[atom_idx].element_id();
  if (type_central < 0 || std::cmp_greater_equal(type_central, num_elements)) {
    return;
  }

  // Build SoA for all neighbors of atom i (once per atom)
  scratch.ensure(neighbor_count);
  for (size_t nb_idx = 0; nb_idx < neighbor_count; ++nb_idx) {
    const auto &r_ij = neighbors[nb_idx].r_ij;
    scratch.nb_x[nb_idx] = r_ij.x();
    scratch.nb_y[nb_idx] = r_ij.y();
    scratch.nb_z[nb_idx] = r_ij.z();
    scratch.nb_dist[nb_idx] = neighbors[nb_idx].distance;
  }

  computeTriadAngles(type_central, atoms, neighbor_count, neighbors, num_elements, scratch, local_tensor);
}

} // namespace

void AngleCalculator::compute(const correlation::core::Cell &cell, const correlation::core::NeighborGraph &graph,
                              AngleTensor &out_angles) {
  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const size_t num_elements = cell.elements().size();

  // Thread-local AngleTensor accumulator
  tbb::enumerable_thread_specific<AngleTensor> ets([&]() {
    return AngleTensor(num_elements, std::vector<std::vector<std::vector<real_t>>>(
                                         num_elements, std::vector<std::vector<real_t>>(num_elements)));
  });

  // Thread-local SoA scratch: pre-built per central atom and reused across
  // atoms on the same thread.  Grows to max(CN) and stays there – zero
  // heap allocs after the first high-CN atom is processed.
  tbb::enumerable_thread_specific<AngleScratch> scratch_ets;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, atom_count), [&](const tbb::blocked_range<size_t> &range) {
    auto &local_tensor = ets.local();
    auto &scratch = scratch_ets.local();

    for (size_t atom_idx = range.begin(); atom_idx != range.end(); ++atom_idx) {
      processCentralAtom(atom_idx, cell, graph, num_elements, scratch, local_tensor);
    }
  });

  // Merge results
  for (const auto &local_tensor : ets) {
    for (size_t i = 0; i < num_elements; ++i) {
      for (size_t j = 0; j < num_elements; ++j) {
        for (size_t k = 0; k < num_elements; ++k) {
          out_angles[i][j][k].insert(out_angles[i][j][k].end(), local_tensor[i][j][k].begin(),
                                     local_tensor[i][j][k].end());
        }
      }
    }
  }
}

} // namespace correlation::calculators
