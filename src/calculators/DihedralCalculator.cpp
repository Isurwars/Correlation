// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/DihedralCalculator.hpp"

#include <cmath>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace calculators {

void DihedralCalculator::compute(const Cell &cell, const NeighborGraph &graph,
                                 DihedralTensor &out_dihedrals) {
  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const size_t num_elements = cell.elements().size();

  // Initialize thread-local storage
  tbb::enumerable_thread_specific<DihedralTensor> ets([&]() {
    return DihedralTensor(
        num_elements,
        std::vector<std::vector<std::vector<std::vector<double>>>>(
            num_elements,
            std::vector<std::vector<std::vector<double>>>(
                num_elements, std::vector<std::vector<double>>(num_elements))));
  });

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, atom_count),
      [&](const tbb::blocked_range<size_t> &r) {
        auto &local_tensor = ets.local();

        // Let atom `i` be atom `B` in the A-B-C-D sequence
        for (size_t B = r.begin(); B != r.end(); ++B) {
          const auto &B_neighbors = graph.getNeighbors(B);
          if (B_neighbors.size() < 2)
            continue;

          const int type_B = atoms[B].element_id();

          // Loop over all neighbors of B, which we consider as `C`
          for (const auto &neighbor_C : B_neighbors) {
            size_t C = neighbor_C.index;

            // To prevent double counting the identical bond B-C as C-B,
            // we enforce an ordering constraint: B < C
            if (B >= C)
              continue;

            const auto &C_neighbors = graph.getNeighbors(C);
            if (C_neighbors.size() < 2)
              continue;

            const int type_C = atoms[C].element_id();
            const linalg::Vector3<double> &r_BC = neighbor_C.r_ij; // B -> C

            // Now, find all A bonded to B (where A != C)
            for (const auto &neighbor_A : B_neighbors) {
              size_t A = neighbor_A.index;
              if (A == C)
                continue;

              const int type_A = atoms[A].element_id();
              const linalg::Vector3<double> &r_BA = neighbor_A.r_ij; // B -> A

              // And find all D bonded to C (where D != B and D != A)
              for (const auto &neighbor_D : C_neighbors) {
                size_t D = neighbor_D.index;
                if (D == B || D == A)
                  continue;

                const int type_D = atoms[D].element_id();

                // Vector C -> D.
                // Wait, r_ij in Neighbor array is Central -> Neighbor.
                // So neighbor_D.r_ij is C -> D.
                const linalg::Vector3<double> &r_CD = neighbor_D.r_ij;

                // We have vectors:
                // b1 = r_BA  (B to A)
                // b2 = r_BC  (B to C)
                // b3 = r_CD  (C to D)
                //
                // The planes are defined by B-C-A and C-B-D? No, standard
                // IUPAC: b1 = r_AB (A to B) = -r_BA b2 = r_BC (B to C) b3 =
                // r_CD (C to D)

                linalg::Vector3<double> b1 = -1.0 * r_BA;
                linalg::Vector3<double> b2 = r_BC;
                linalg::Vector3<double> b3 = r_CD;

                // n1 = b1 x b2
                // n2 = b2 x b3
                linalg::Vector3<double> n1 = linalg::cross(b1, b2);
                linalg::Vector3<double> n2 = linalg::cross(b2, b3);

                double n1_norm = linalg::norm(n1);
                double n2_norm = linalg::norm(n2);

                if (n1_norm < 1e-8 || n2_norm < 1e-8) {
                  continue; // Collinear atoms, dihedral undefined.
                }

                // Normalizing
                n1 = n1 * (1.0 / n1_norm);
                n2 = n2 * (1.0 / n2_norm);

                // m = n1 x (b2 / |b2|)
                linalg::Vector3<double> b2_hat = b2 * (1.0 / linalg::norm(b2));
                linalg::Vector3<double> m = linalg::cross(n1, b2_hat);

                // cos(phi) = n1 . n2
                // sin(phi) = m . n2
                double x = linalg::dot(n1, n2);
                double y = linalg::dot(m, n2);

                double dihedral_angle =
                    std::atan2(y, x); // Returns range [-pi, pi]

                local_tensor[type_A][type_B][type_C][type_D].push_back(
                    dihedral_angle);

                // For unordered element types, depending on symmetry you might
                // also store reverse D-C-B-A. However to keep it clean, we
                // store A-B-C-D only, and the parser can query both later.
                if (type_A != type_D || type_B != type_C) {
                  local_tensor[type_D][type_C][type_B][type_A].push_back(
                      dihedral_angle);
                }
              }
            }
          }
        }
      });

  // Merge results
  for (const auto &local_tensor : ets) {
    for (size_t a = 0; a < num_elements; ++a) {
      for (size_t b = 0; b < num_elements; ++b) {
        for (size_t c = 0; c < num_elements; ++c) {
          for (size_t d = 0; d < num_elements; ++d) {
            out_dihedrals[a][b][c][d].insert(out_dihedrals[a][b][c][d].end(),
                                             local_tensor[a][b][c][d].begin(),
                                             local_tensor[a][b][c][d].end());
          }
        }
      }
    }
  }
}

} // namespace calculators
