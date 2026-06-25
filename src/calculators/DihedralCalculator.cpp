/**
 * @file DihedralCalculator.cpp
 * @brief Implementation of dihedral angle geometry.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/DihedralCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "math/LinearAlgebra.hpp"

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <optional>

namespace correlation::calculators {

void DihedralCalculator::calculateFrame(correlation::analysis::DistributionFunctions & /*dists*/,
                                        const correlation::analysis::AnalysisSettings & /*settings*/) const {
  // DihedralCalculator is a foundational calculator — called by
  // StructureAnalyzer during its construction. Nothing to do here.
}

namespace {

/**
 * @brief Bundled displacement vectors for calculating a dihedral angle.
 */
struct DihedralVectors {
  correlation::math::Vector3<double> r_ba;
  correlation::math::Vector3<double> r_bc;
  correlation::math::Vector3<double> r_cd;
};

/**
 * @brief Helper to calculate a dihedral angle from displacement vectors.
 */
std::optional<double> calculateDihedralAngle(const DihedralVectors &vectors) {
  correlation::math::Vector3<double> const vec_b1 = -1.0 * vectors.r_ba;
  correlation::math::Vector3<double> const vec_b2 = vectors.r_bc;
  correlation::math::Vector3<double> const vec_b3 = vectors.r_cd;

  correlation::math::Vector3<double> normal_n1 = correlation::math::cross(vec_b1, vec_b2);
  correlation::math::Vector3<double> normal_n2 = correlation::math::cross(vec_b2, vec_b3);

  double const n1_norm = correlation::math::norm(normal_n1);
  double const n2_norm = correlation::math::norm(normal_n2);

  if (n1_norm < 1e-8 || n2_norm < 1e-8) {
    return std::nullopt; // Collinear atoms, dihedral undefined.
  }

  normal_n1 = correlation::math::normalize(normal_n1);
  normal_n2 = correlation::math::normalize(normal_n2);

  double const b2_norm = correlation::math::norm(vec_b2);
  if (b2_norm < 1e-8) {
    return std::nullopt; // Coincident central bond, dihedral undefined.
  }
  correlation::math::Vector3<double> const b2_hat = vec_b2 / b2_norm;
  correlation::math::Vector3<double> const vec_m = correlation::math::cross(normal_n1, b2_hat);

  double const dot_x = correlation::math::dot(normal_n1, normal_n2);
  double const dot_y = correlation::math::dot(vec_m, normal_n2);

  return std::atan2(dot_y, dot_x);
}

/**
 * @brief Helper to find and bin dihedral angles for a specific central pair B-C.
 */
void findAndProcessDihedralsForPair(size_t idx_b, size_t idx_c, const correlation::core::NeighborGraph &graph,
                                    const std::vector<correlation::core::Atom> &atoms,
                                    const correlation::core::Neighbor &neighbor_c,
                                    correlation::analysis::StructureAnalyzer::DihedralTensor &local_tensor) {
  const auto &neighbors_b = graph.getNeighbors(idx_b);
  const auto &neighbors_c = graph.getNeighbors(idx_c);

  const int type_b = atoms[idx_b].element_id();
  const int type_c = atoms[idx_c].element_id();
  const correlation::math::Vector3<double> &r_bc = neighbor_c.r_ij;

  for (const auto &neighbor_a : neighbors_b) {
    size_t const idx_a = neighbor_a.index;
    if (idx_a == idx_c) {
      continue;
    }

    const int type_a = atoms[idx_a].element_id();
    const correlation::math::Vector3<double> &r_ba = neighbor_a.r_ij;

    for (const auto &neighbor_d : neighbors_c) {
      size_t const idx_d = neighbor_d.index;
      if (idx_d == idx_b || idx_d == idx_a) {
        continue;
      }

      const int type_d = atoms[idx_d].element_id();
      const correlation::math::Vector3<double> &r_cd = neighbor_d.r_ij;

      auto const dihedral_angle_opt = calculateDihedralAngle({.r_ba = r_ba, .r_bc = r_bc, .r_cd = r_cd});
      if (!dihedral_angle_opt.has_value()) {
        continue;
      }
      double const dihedral_angle = dihedral_angle_opt.value();

      local_tensor[type_a][type_b][type_c][type_d].push_back(dihedral_angle);

      if (type_a != type_d || type_b != type_c) {
        local_tensor[type_d][type_c][type_b][type_a].push_back(dihedral_angle);
      }
    }
  }
}

/**
 * @brief Helper to merge thread-local tensors into the final dihedral tensor.
 */
void mergeThreadLocalTensors(
    correlation::analysis::StructureAnalyzer::DihedralTensor &out_dihedrals,
    const tbb::enumerable_thread_specific<correlation::analysis::StructureAnalyzer::DihedralTensor> &ets,
    size_t num_elements) {
  for (const auto &local_tensor : ets) {
    for (size_t idx_a = 0; idx_a < num_elements; ++idx_a) {
      for (size_t idx_b = 0; idx_b < num_elements; ++idx_b) {
        for (size_t idx_c = 0; idx_c < num_elements; ++idx_c) {
          for (size_t idx_d = 0; idx_d < num_elements; ++idx_d) {
            out_dihedrals[idx_a][idx_b][idx_c][idx_d].insert(out_dihedrals[idx_a][idx_b][idx_c][idx_d].end(),
                                                             local_tensor[idx_a][idx_b][idx_c][idx_d].begin(),
                                                             local_tensor[idx_a][idx_b][idx_c][idx_d].end());
          }
        }
      }
    }
  }
}

} // namespace

void DihedralCalculator::compute(const correlation::core::Cell &cell, const correlation::core::NeighborGraph &graph,
                                 correlation::analysis::StructureAnalyzer::DihedralTensor &out_dihedrals) {
  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const size_t num_elements = cell.elements().size();

  // Initialize thread-local storage
  tbb::enumerable_thread_specific<correlation::analysis::StructureAnalyzer::DihedralTensor> ets([&]() {
    return correlation::analysis::StructureAnalyzer::DihedralTensor(
        num_elements, std::vector<std::vector<std::vector<std::vector<double>>>>(
                          num_elements, std::vector<std::vector<std::vector<double>>>(
                                            num_elements, std::vector<std::vector<double>>(num_elements))));
  });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, atom_count), [&](const tbb::blocked_range<size_t> &blocked_range) {
    auto &local_tensor = ets.local();

    // Let atom `idx_b` be atom `B` in the A-B-C-D sequence
    for (size_t idx_b = blocked_range.begin(); idx_b != blocked_range.end(); ++idx_b) {
      const auto &neighbors_b = graph.getNeighbors(idx_b);
      if (neighbors_b.size() < 2) {
        continue;
      }

      // Loop over all neighbors of B, which we consider as `C`
      for (const auto &neighbor_c : neighbors_b) {
        size_t const idx_c = neighbor_c.index;

        // To prevent double counting the identical bond B-C as C-B,
        // we enforce an ordering constraint: B < C
        if (idx_b >= idx_c) {
          continue;
        }

        const auto &neighbors_c = graph.getNeighbors(idx_c);
        if (neighbors_c.size() < 2) {
          continue;
        }

        findAndProcessDihedralsForPair(idx_b, idx_c, graph, atoms, neighbor_c, local_tensor);
      }
    }
  });

  mergeThreadLocalTensors(out_dihedrals, ets, num_elements);
}

} // namespace correlation::calculators
