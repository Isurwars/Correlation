/**
 * @file DistanceCalculator.cpp
 * @brief Implementation of pairwise distance calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/DistanceCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/SIMDUtils.hpp"

#include <cmath>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for_each.h>

namespace correlation::calculators {

void DistanceCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                        const correlation::analysis::AnalysisSettings &settings) const {}

namespace {
struct SearchGridConfig {
  int K_x;
  int K_y;
  int K_z;
  int max_dx;
  int max_dy;
  int max_dz;
};

/**
 * @struct SoACoordinates
 * @brief Structure of Arrays (SoA) view of wrapped Cartesian positions to prevent parameter swapping.
 */
struct SoACoordinates {
  const double *x;
  const double *y;
  const double *z;
};

/**
 * @struct FlatCellList
 * @brief Reference view of the contiguous flat cell list offsets and indices.
 */
struct FlatCellList {
  const size_t *offsets;
  const size_t *indices;
};

/**
 * @struct ThreadLocalBond
 * @brief Represents a bond/directed edge in flat thread-local storage.
 */
struct ThreadLocalBond {
  size_t from;
  size_t to;
  double distance;
  correlation::math::Vector3<double> r_ij;
};

/**
 * @struct ThreadLocalConfig
 * @brief Thread-local configuration parameters.
 */
struct ThreadLocalConfig {
  size_t num_elements;
  size_t atom_count;
};

/**
 * @struct ThreadLocalDistances
 * @brief Thread-local storage for parallel distance calculations.
 *
 * Contains local tensor results, flat bond list, and pre-allocated scratch buffers
 * for SoA (Structure of Arrays) SIMD operations.
 */
struct ThreadLocalDistances {
  DistanceTensor distance_tensor_local; ///< Local partial distance tensor.
  std::vector<ThreadLocalBond> bonds;   ///< Local flat bond/neighbor lists.

  // Per-thread SoA scratch buffers (reused and pre-reserved to avoid allocations)
  std::vector<double> soa_x;       ///< Scratch for x-coordinates.
  std::vector<double> soa_y;       ///< Scratch for y-coordinates.
  std::vector<double> soa_z;       ///< Scratch for z-coordinates.
  std::vector<double> dsq_scratch; ///< Scratch for squared distance results.
  std::vector<size_t> candidate_j; ///< Scratch for candidate indices.
  size_t candidate_count = 0;      ///< Count of candidates populated in scratch buffers.

  /**
   * @brief Constructs thread-local storage.
   * @param config The configuration parameters.
   */
  explicit ThreadLocalDistances(ThreadLocalConfig config)
      : distance_tensor_local(config.num_elements, std::vector<std::vector<double>>(config.num_elements)) {
    // Pre-reserve capacities to avoid resizing checks in hot loops
    soa_x.reserve(1024);
    soa_y.reserve(1024);
    soa_z.reserve(1024);
    candidate_j.reserve(1024);
    dsq_scratch.reserve(1024);
    bonds.reserve(2048);
  }

  /**
   * @brief Gathers candidate atoms from neighboring bins.
   */
  void collectCandidates(size_t atom_idx, const SoACoordinates &coords, const std::vector<int> &atom_bin,
                         const FlatCellList &cell_list, const correlation::math::Matrix3<double> &lattice,
                         SearchGridConfig grid_config, bool ignore_periodic_self_interactions);

  /**
   * @brief Gathers candidate atoms from a specific bin.
   */
  void collectCandidatesFromBin(size_t atom_idx, const correlation::math::Vector3<double> &disp, int n_bin_idx,
                                const FlatCellList &cell_list, bool zero_disp, const SoACoordinates &coords,
                                bool ignore_periodic_self_interactions, size_t &c_count);

  /**
   * @brief Computes SIMD pairwise distances and updates local structures.
   */
  void computeDistances(size_t atom_idx, const std::vector<correlation::core::Atom> &atoms, const SoACoordinates &coords,
                        double cutoff_sq, const std::vector<std::vector<double>> &bond_cutoffs_sq);
};

void ThreadLocalDistances::collectCandidates(size_t atom_idx, const SoACoordinates &coords,
                                             const std::vector<int> &atom_bin, const FlatCellList &cell_list,
                                             const correlation::math::Matrix3<double> &lattice,
                                             SearchGridConfig grid_config, bool ignore_periodic_self_interactions) {

  size_t c_count = 0;
  int c_x = atom_bin[atom_idx] / (grid_config.K_y * grid_config.K_z);
  int c_y = (atom_bin[atom_idx] / grid_config.K_z) % grid_config.K_y;
  int c_z = atom_bin[atom_idx] % grid_config.K_z;

  for (int dx = -grid_config.max_dx; dx <= grid_config.max_dx; ++dx) {
    int nx_bin = c_x + dx;
    int shift_x = (nx_bin >= 0) ? (nx_bin / grid_config.K_x) : ((nx_bin - grid_config.K_x + 1) / grid_config.K_x);
    int wrap_x = nx_bin - shift_x * grid_config.K_x;

    for (int dy = -grid_config.max_dy; dy <= grid_config.max_dy; ++dy) {
      int ny_bin = c_y + dy;
      int shift_y = (ny_bin >= 0) ? (ny_bin / grid_config.K_y) : ((ny_bin - grid_config.K_y + 1) / grid_config.K_y);
      int wrap_y = ny_bin - shift_y * grid_config.K_y;

      for (int dz = -grid_config.max_dz; dz <= grid_config.max_dz; ++dz) {
        int nz_bin = c_z + dz;
        int shift_z = (nz_bin >= 0) ? (nz_bin / grid_config.K_z) : ((nz_bin - grid_config.K_z + 1) / grid_config.K_z);
        int wrap_z = nz_bin - shift_z * grid_config.K_z;

        int n_bin_idx = wrap_x * (grid_config.K_y * grid_config.K_z) + wrap_y * grid_config.K_z + wrap_z;

        // Use FMA for precise displacement coordinate calculation
        double disp_x = std::fma(shift_z, lattice[2].x(), std::fma(shift_y, lattice[1].x(), shift_x * lattice[0].x()));
        double disp_y = std::fma(shift_z, lattice[2].y(), std::fma(shift_y, lattice[1].y(), shift_x * lattice[0].y()));
        double disp_z = std::fma(shift_z, lattice[2].z(), std::fma(shift_y, lattice[1].z(), shift_x * lattice[0].z()));
        correlation::math::Vector3<double> disp = {disp_x, disp_y, disp_z};

        bool zero_disp = (shift_x == 0 && shift_y == 0 && shift_z == 0);

        collectCandidatesFromBin(atom_idx, disp, n_bin_idx, cell_list, zero_disp, coords,
                                 ignore_periodic_self_interactions, c_count);
      }
    }
  }
  candidate_count = c_count;
}

void ThreadLocalDistances::collectCandidatesFromBin(size_t atom_idx, const correlation::math::Vector3<double> &disp,
                                                    int n_bin_idx, const FlatCellList &cell_list, bool zero_disp,
                                                    const SoACoordinates &coords, bool ignore_periodic_self_interactions,
                                                    size_t &c_count) {

  size_t const start = cell_list.offsets[n_bin_idx];
  size_t const end = cell_list.offsets[n_bin_idx + 1];

  for (size_t offset = start; offset < end; ++offset) {
    size_t const j_idx = cell_list.indices[offset];
    if (j_idx < atom_idx) {
      continue;
    }
    if (atom_idx == j_idx && zero_disp) {
      continue;
    }
    if (atom_idx == j_idx && ignore_periodic_self_interactions) {
      continue;
    }

    double const shifted_x = coords.x[j_idx] + disp.x();
    double const shifted_y = coords.y[j_idx] + disp.y();
    double const shifted_z = coords.z[j_idx] + disp.z();

    if (soa_x.size() <= c_count) {
      soa_x.push_back(shifted_x);
      soa_y.push_back(shifted_y);
      soa_z.push_back(shifted_z);
      candidate_j.push_back(j_idx);
      dsq_scratch.push_back(0.0);
    } else {
      soa_x[c_count] = shifted_x;
      soa_y[c_count] = shifted_y;
      soa_z[c_count] = shifted_z;
      candidate_j[c_count] = j_idx;
    }
    c_count++;
  }
}

void ThreadLocalDistances::computeDistances(size_t atom_idx, const std::vector<correlation::core::Atom> &atoms,
                                            const SoACoordinates &coords, double cutoff_sq,
                                            const std::vector<std::vector<double>> &bond_cutoffs_sq) {

  size_t const c_count = candidate_count;
  if (c_count > 0) {
    const auto &atom_A = atoms[atom_idx];
    const int type_A = atom_A.element_id();
    const double a_x = coords.x[atom_idx];
    const double a_y = coords.y[atom_idx];
    const double a_z = coords.z[atom_idx];

    correlation::math::PositionBlock block{.x = soa_x.data(), .y = soa_y.data(), .z = soa_z.data(), .count = c_count};
    correlation::math::compute_dsq_block(a_x, a_y, a_z, block, dsq_scratch.data());

    for (size_t k = 0; k < c_count; ++k) {
      double const d_sq = dsq_scratch[k];
      if (d_sq >= cutoff_sq) {
        continue;
      }

      size_t const j_idx = candidate_j[k];
      double const dist = std::sqrt(d_sq);
      int const type_B = atoms[j_idx].element_id();

      distance_tensor_local[type_A][type_B].push_back(dist);
      if (atom_idx != j_idx && type_A != type_B) {
        distance_tensor_local[type_B][type_A].push_back(dist);
      }

      double max_bond_dist_sq = 0.0;
      if (static_cast<size_t>(type_A) < bond_cutoffs_sq.size() &&
          static_cast<size_t>(type_B) < bond_cutoffs_sq[type_A].size()) {
        max_bond_dist_sq = bond_cutoffs_sq[type_A][type_B];
      }

      if (d_sq <= max_bond_dist_sq) {
        // Recover displacement vector coordinates precisely to minimize rounding/cancellation errors
        double const disp_x = soa_x[k] - coords.x[j_idx];
        double const disp_y = soa_y[k] - coords.y[j_idx];
        double const disp_z = soa_z[k] - coords.z[j_idx];

        double const r_x = (coords.x[j_idx] - a_x) + disp_x;
        double const r_y = (coords.y[j_idx] - a_y) + disp_y;
        double const r_z = (coords.z[j_idx] - a_z) + disp_z;

        // Recompute consistent distance directly from the exact relative vector
        double const precise_dist = std::sqrt(r_x * r_x + r_y * r_y + r_z * r_z);

        correlation::math::Vector3<double> r_ij = {r_x, r_y, r_z};
        bonds.push_back({.from = atom_idx, .to = j_idx, .distance = precise_dist, .r_ij = r_ij});
        if (atom_idx != j_idx) {
          bonds.push_back({.from = j_idx, .to = atom_idx, .distance = precise_dist, .r_ij = -1.0 * r_ij});
        }
      }
    }
  }
}
} // namespace

void DistanceCalculator::compute(const correlation::core::Cell &cell, double cutoff_sq,
                                 const std::vector<std::vector<double>> &bond_cutoffs_sq,
                                 bool ignore_periodic_self_interactions, DistanceTensor &out_distances,
                                 correlation::core::NeighborGraph &out_graph) {

  if (cutoff_sq <= 0.0) {
    throw std::invalid_argument("Cutoff squared must be strictly positive.");
  }

  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const size_t num_elements = cell.elements().size();
  const auto &lattice = cell.latticeVectors();

  // New Cell-List logic
  math::Vector3<double> lattice_a = lattice[0];
  math::Vector3<double> lattice_b = lattice[1];
  math::Vector3<double> lattice_c = lattice[2];

  double vol = cell.volume();
  double width_x = vol / correlation::math::norm(correlation::math::cross(lattice_b, lattice_c));
  double width_y = vol / correlation::math::norm(correlation::math::cross(lattice_a, lattice_c));
  double width_z = vol / correlation::math::norm(correlation::math::cross(lattice_a, lattice_b));

  double cutoff = std::sqrt(cutoff_sq);

  int K_x = std::max(1, static_cast<int>(std::floor(width_x / cutoff)));
  int K_y = std::max(1, static_cast<int>(std::floor(width_y / cutoff)));
  int K_z = std::max(1, static_cast<int>(std::floor(width_z / cutoff)));

  int max_dx = (K_x == 1) ? static_cast<int>(std::ceil(cutoff / width_x)) : 1;
  int max_dy = (K_y == 1) ? static_cast<int>(std::ceil(cutoff / width_y)) : 1;
  int max_dz = (K_z == 1) ? static_cast<int>(std::ceil(cutoff / width_z)) : 1;

  if (max_dx + max_dy + max_dz > 8) {
    ignore_periodic_self_interactions = false;
  }

  auto num_bins = static_cast<size_t>(K_x) * static_cast<size_t>(K_y) * static_cast<size_t>(K_z);
  std::vector<size_t> bin_counts(num_bins, 0);
  std::vector<int> atom_bin(atom_count);

  std::vector<double> wrapped_x(atom_count);
  std::vector<double> wrapped_y(atom_count);
  std::vector<double> wrapped_z(atom_count);

  for (size_t i = 0; i < atom_count; ++i) {
    correlation::math::Vector3<double> frac = cell.inverseLatticeVectors() * atoms[i].position();
    double f_x = frac.x() - std::floor(frac.x());
    double f_y = frac.y() - std::floor(frac.y());
    double f_z = frac.z() - std::floor(frac.z());

    // Use std::fma for extremely precise coordinate mapping
    wrapped_x[i] = std::fma(f_z, lattice[2].x(), std::fma(f_y, lattice[1].x(), f_x * lattice[0].x()));
    wrapped_y[i] = std::fma(f_z, lattice[2].y(), std::fma(f_y, lattice[1].y(), f_x * lattice[0].y()));
    wrapped_z[i] = std::fma(f_z, lattice[2].z(), std::fma(f_y, lattice[1].z(), f_x * lattice[0].z()));

    int c_x = std::clamp(static_cast<int>(std::floor(f_x * K_x)), 0, K_x - 1);
    int c_y = std::clamp(static_cast<int>(std::floor(f_y * K_y)), 0, K_y - 1);
    int c_z = std::clamp(static_cast<int>(std::floor(f_z * K_z)), 0, K_z - 1);

    int bin_idx = c_x * (K_y * K_z) + c_y * K_z + c_z;
    atom_bin[i] = bin_idx;
    bin_counts[bin_idx]++;
  }

  // Build prefix-sum offsets for flat cell list
  std::vector<size_t> bin_offsets(num_bins + 1, 0);
  for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
    bin_offsets[bin_idx + 1] = bin_offsets[bin_idx] + bin_counts[bin_idx];
  }

  // Populate flat bin indices
  std::vector<size_t> bin_indices(atom_count);
  std::vector<size_t> insertion_cursors = bin_offsets;
  for (size_t i = 0; i < atom_count; ++i) {
    int bin_idx = atom_bin[i];
    bin_indices[insertion_cursors[bin_idx]++] = i;
  }

  SearchGridConfig grid_config{
      .K_x = K_x, .K_y = K_y, .K_z = K_z, .max_dx = max_dx, .max_dy = max_dy, .max_dz = max_dz};

  SoACoordinates coords{
      .x = wrapped_x.data(),
      .y = wrapped_y.data(),
      .z = wrapped_z.data()
  };

  FlatCellList cell_list{
      .offsets = bin_offsets.data(),
      .indices = bin_indices.data()
  };

  tbb::enumerable_thread_specific<ThreadLocalDistances> ets(
      ThreadLocalConfig{.num_elements = num_elements, .atom_count = atom_count});

  tbb::parallel_for(tbb::blocked_range<size_t>(0, atom_count), [&](const tbb::blocked_range<size_t> &range) {
    ThreadLocalDistances &local_results = ets.local();
    for (size_t i = range.begin(); i != range.end(); ++i) {
      local_results.collectCandidates(i, coords, atom_bin, cell_list, lattice, grid_config, ignore_periodic_self_interactions);
      local_results.computeDistances(i, atoms, coords, cutoff_sq, bond_cutoffs_sq);
    }
  });

  for (auto &local_results : ets) {
    for (size_t i = 0; i < num_elements; ++i) {
      for (size_t j_idx = 0; j_idx < num_elements; ++j_idx) {
        out_distances[i][j_idx].insert(out_distances[i][j_idx].end(),
                                       local_results.distance_tensor_local[i][j_idx].begin(),
                                       local_results.distance_tensor_local[i][j_idx].end());
      }
    }
    for (const auto &bond : local_results.bonds) {
      out_graph.addDirectedEdge(bond.from, bond.to, bond.distance, bond.r_ij);
    }
  }
}

} // namespace correlation::calculators
