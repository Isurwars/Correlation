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
#include "math/Precision.hpp"
#include "math/SIMDUtils.hpp"

#if defined(CORRELATION_USE_CUDA) || defined(CORRELATION_USE_HIP)
#include "calculators/GPUDistanceCalculator.hpp"
#endif

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
  const real_t *x;
  const real_t *y;
  const real_t *z;
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
  real_t distance;
  correlation::math::Vector3<real_t> r_ij;
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
  std::vector<real_t> soa_x;       ///< Scratch for x-coordinates.
  std::vector<real_t> soa_y;       ///< Scratch for y-coordinates.
  std::vector<real_t> soa_z;       ///< Scratch for z-coordinates.
  std::vector<real_t> dsq_scratch; ///< Scratch for squared distance results.
  std::vector<size_t> candidate_j; ///< Scratch for candidate indices.
  size_t candidate_count = 0;      ///< Count of candidates populated in scratch buffers.

  /**
   * @brief Constructs thread-local storage.
   * @param config The configuration parameters.
   */
  explicit ThreadLocalDistances(ThreadLocalConfig config)
      : distance_tensor_local(config.num_elements, std::vector<std::vector<real_t>>(config.num_elements)) {
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
                         const FlatCellList &cell_list, const correlation::math::Matrix3<real_t> &lattice,
                         SearchGridConfig grid_config, bool ignore_periodic_self_interactions);

  /**
   * @brief Gathers candidate atoms from a specific bin.
   */
  void collectCandidatesFromBin(size_t atom_idx, const correlation::math::Vector3<real_t> &disp, int n_bin_idx,
                                const FlatCellList &cell_list, bool zero_disp, const SoACoordinates &coords,
                                bool ignore_periodic_self_interactions, size_t &c_count);

  /**
   * @brief Computes SIMD pairwise distances and updates local structures.
   */
  void computeDistances(size_t atom_idx, const std::vector<correlation::core::Atom> &atoms,
                        const SoACoordinates &coords, real_t cutoff_sq,
                        const correlation::analysis::BondCutoffMatrix &bond_cutoffs);
};

void ThreadLocalDistances::collectCandidates(size_t atom_idx, const SoACoordinates &coords,
                                             const std::vector<int> &atom_bin, const FlatCellList &cell_list,
                                             const correlation::math::Matrix3<real_t> &lattice,
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
        real_t disp_x = std::fma(
            static_cast<real_t>(shift_z), lattice[2].x(),
            std::fma(static_cast<real_t>(shift_y), lattice[1].x(), static_cast<real_t>(shift_x) * lattice[0].x()));
        real_t disp_y = std::fma(
            static_cast<real_t>(shift_z), lattice[2].y(),
            std::fma(static_cast<real_t>(shift_y), lattice[1].y(), static_cast<real_t>(shift_x) * lattice[0].y()));
        real_t disp_z = std::fma(
            static_cast<real_t>(shift_z), lattice[2].z(),
            std::fma(static_cast<real_t>(shift_y), lattice[1].z(), static_cast<real_t>(shift_x) * lattice[0].z()));
        correlation::math::Vector3<real_t> disp = {disp_x, disp_y, disp_z};

        bool zero_disp = (shift_x == 0 && shift_y == 0 && shift_z == 0);

        collectCandidatesFromBin(atom_idx, disp, n_bin_idx, cell_list, zero_disp, coords,
                                 ignore_periodic_self_interactions, c_count);
      }
    }
  }
  candidate_count = c_count;
}

void ThreadLocalDistances::collectCandidatesFromBin(size_t atom_idx, const correlation::math::Vector3<real_t> &disp,
                                                    int n_bin_idx, const FlatCellList &cell_list, bool zero_disp,
                                                    const SoACoordinates &coords,
                                                    bool ignore_periodic_self_interactions, size_t &c_count) {

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

    real_t const shifted_x = coords.x[j_idx] + disp.x();
    real_t const shifted_y = coords.y[j_idx] + disp.y();
    real_t const shifted_z = coords.z[j_idx] + disp.z();

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
                                            const SoACoordinates &coords, real_t cutoff_sq,
                                            const correlation::analysis::BondCutoffMatrix &bond_cutoffs) {

  size_t const c_count = candidate_count;
  if (c_count > 0) {
    const auto &atom_A = atoms[atom_idx];
    const int type_A = atom_A.element_id();
    const real_t a_x = coords.x[atom_idx];
    const real_t a_y = coords.y[atom_idx];
    const real_t a_z = coords.z[atom_idx];

    correlation::math::PositionBlock block{
        .x = soa_x.data(),
        .y = soa_y.data(),
        .z = soa_z.data(),
        .count = c_count,
    };
    correlation::math::compute_dsq_block(a_x, a_y, a_z, block, dsq_scratch.data());

    for (size_t k = 0; k < c_count; ++k) {
      real_t const d_sq = dsq_scratch[k];
      if (d_sq >= cutoff_sq) {
        continue;
      }

      size_t const j_idx = candidate_j[k];
      auto const dist = static_cast<real_t>(std::sqrt(d_sq));
      int const type_B = atoms[j_idx].element_id();

      distance_tensor_local[type_A][type_B].push_back(dist);
      if (atom_idx != j_idx && type_A != type_B) {
        distance_tensor_local[type_B][type_A].push_back(dist);
      }

      real_t min_bond_dist_sq = 0.0;
      real_t max_bond_dist_sq = 0.0;
      if (static_cast<size_t>(type_A) < bond_cutoffs.size() &&
          static_cast<size_t>(type_B) < bond_cutoffs[type_A].size()) {
        const auto &range = bond_cutoffs[type_A][type_B];
        min_bond_dist_sq = range.min_sq;
        max_bond_dist_sq = range.max_sq;
      }

      if (d_sq >= min_bond_dist_sq && d_sq <= max_bond_dist_sq) {
        // Recover displacement vector coordinates precisely to minimize rounding/cancellation errors
        real_t const disp_x = soa_x[k] - coords.x[j_idx];
        real_t const disp_y = soa_y[k] - coords.y[j_idx];
        real_t const disp_z = soa_z[k] - coords.z[j_idx];

        real_t const r_x = (coords.x[j_idx] - a_x) + disp_x;
        real_t const r_y = (coords.y[j_idx] - a_y) + disp_y;
        real_t const r_z = (coords.z[j_idx] - a_z) + disp_z;

        // Recompute consistent distance directly from the exact relative vector
        real_t const precise_dist = std::sqrt(r_x * r_x + r_y * r_y + r_z * r_z);

        correlation::math::Vector3<real_t> r_ij = {r_x, r_y, r_z};

        bonds.push_back(ThreadLocalBond{
            .from = atom_idx,
            .to = j_idx,
            .distance = precise_dist,
            .r_ij = r_ij,
        });

        if (atom_idx != j_idx) {
          bonds.push_back(ThreadLocalBond{
              .from = j_idx,
              .to = atom_idx,
              .distance = precise_dist,
              .r_ij = -1.0 * r_ij,
          });
        }
      }
    }
  }
}
/**
 * @brief Validates input parameters for distance computation.
 */
void validateInputs(real_t cutoff_sq, const correlation::analysis::BondCutoffMatrix &bond_cutoffs) {
  if (cutoff_sq <= 0.0) {
    throw std::invalid_argument("Cutoff squared must be strictly positive.");
  }

  for (const auto &row : bond_cutoffs) {
    for (const auto &range : row) {
      if (range.min_sq < 0.0 || range.max_sq < 0.0) {
        throw std::invalid_argument("Bond cutoff bounds cannot be negative.");
      }
      if (range.min_sq > range.max_sq) {
        throw std::invalid_argument("Minimum bond cutoff cannot be greater than maximum bond cutoff.");
      }
    }
  }
}

#if defined(CORRELATION_USE_CUDA) || defined(CORRELATION_USE_HIP)
/**
 * @brief Dispatches calculation to GPU if a compatible device is present.
 */
bool tryComputeGpu(const correlation::core::Cell &cell, real_t cutoff_sq,
                   const correlation::analysis::BondCutoffMatrix &bond_cutoffs, bool ignore_periodic_self_interactions,
                   DistanceTensor &out_distances, correlation::core::NeighborGraph &out_graph) {
  if (gpu::has_gpu_device()) {
    std::vector<std::vector<real_t>> max_cutoffs_sq(bond_cutoffs.size());
    for (size_t i = 0; i < bond_cutoffs.size(); ++i) {
      max_cutoffs_sq[i].resize(bond_cutoffs[i].size());
      for (size_t j = 0; j < bond_cutoffs[i].size(); ++j) {
        max_cutoffs_sq[i][j] = bond_cutoffs[i][j].max_sq;
      }
    }
    gpu::compute_distances_gpu(cell, cutoff_sq, max_cutoffs_sq, ignore_periodic_self_interactions, out_distances,
                               out_graph);
    return true;
  }
  return false;
}
#endif

/**
 * @brief Computes search grid dimensions and cutoff parameters.
 */
SearchGridConfig buildSearchGridConfig(const correlation::core::Cell &cell, real_t cutoff_sq,
                                       bool &ignore_periodic_self_interactions) {
  const auto &lattice = cell.latticeVectors();
  const math::Vector3<real_t> lattice_a = lattice[0];
  const math::Vector3<real_t> lattice_b = lattice[1];
  const math::Vector3<real_t> lattice_c = lattice[2];

  const real_t vol = cell.volume();
  const real_t width_x = vol / correlation::math::norm(correlation::math::cross(lattice_b, lattice_c));
  const real_t width_y = vol / correlation::math::norm(correlation::math::cross(lattice_a, lattice_c));
  const real_t width_z = vol / correlation::math::norm(correlation::math::cross(lattice_a, lattice_b));

  const real_t cutoff = std::sqrt(cutoff_sq);

  const int K_x = std::max(1, static_cast<int>(std::floor(width_x / cutoff)));
  const int K_y = std::max(1, static_cast<int>(std::floor(width_y / cutoff)));
  const int K_z = std::max(1, static_cast<int>(std::floor(width_z / cutoff)));

  const int max_dx = (K_x == 1) ? static_cast<int>(std::ceil(cutoff / width_x)) : 1;
  const int max_dy = (K_y == 1) ? static_cast<int>(std::ceil(cutoff / width_y)) : 1;
  const int max_dz = (K_z == 1) ? static_cast<int>(std::ceil(cutoff / width_z)) : 1;

  if (max_dx + max_dy + max_dz > 8) {
    ignore_periodic_self_interactions = false;
  }

  return SearchGridConfig{
      .K_x = K_x,
      .K_y = K_y,
      .K_z = K_z,
      .max_dx = max_dx,
      .max_dy = max_dy,
      .max_dz = max_dz,
  };
}

/**
 * @struct WrappedPositions
 * @brief Holds wrapped Cartesian coordinates and cell-bin metadata.
 */
struct WrappedPositions {
  std::vector<real_t> wrapped_x;
  std::vector<real_t> wrapped_y;
  std::vector<real_t> wrapped_z;
  std::vector<int> atom_bin;
  std::vector<size_t> bin_counts;
};

/**
 * @brief Computes wrapped coordinates and bin distributions for cell atoms.
 */
WrappedPositions buildWrappedPositionsAndBins(const correlation::core::Cell &cell,
                                              const SearchGridConfig &grid_config) {
  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const auto &lattice = cell.latticeVectors();
  const auto num_bins = static_cast<size_t>(grid_config.K_x) * static_cast<size_t>(grid_config.K_y) *
                        static_cast<size_t>(grid_config.K_z);

  WrappedPositions pos{
      .wrapped_x = std::vector<real_t>(atom_count),
      .wrapped_y = std::vector<real_t>(atom_count),
      .wrapped_z = std::vector<real_t>(atom_count),
      .atom_bin = std::vector<int>(atom_count),
      .bin_counts = std::vector<size_t>(num_bins, 0),
  };

  for (size_t i = 0; i < atom_count; ++i) {
    const correlation::math::Vector3<real_t> frac = cell.inverseLatticeVectors() * atoms[i].position();
    const real_t f_x = frac.x() - std::floor(frac.x());
    const real_t f_y = frac.y() - std::floor(frac.y());
    const real_t f_z = frac.z() - std::floor(frac.z());

    pos.wrapped_x[i] = std::fma(f_z, lattice[2].x(), std::fma(f_y, lattice[1].x(), f_x * lattice[0].x()));
    pos.wrapped_y[i] = std::fma(f_z, lattice[2].y(), std::fma(f_y, lattice[1].y(), f_x * lattice[0].y()));
    pos.wrapped_z[i] = std::fma(f_z, lattice[2].z(), std::fma(f_y, lattice[1].z(), f_x * lattice[0].z()));

    const int c_x =
        std::clamp(static_cast<int>(std::floor(f_x * static_cast<real_t>(grid_config.K_x))), 0, grid_config.K_x - 1);
    const int c_y =
        std::clamp(static_cast<int>(std::floor(f_y * static_cast<real_t>(grid_config.K_y))), 0, grid_config.K_y - 1);
    const int c_z =
        std::clamp(static_cast<int>(std::floor(f_z * static_cast<real_t>(grid_config.K_z))), 0, grid_config.K_z - 1);

    const int bin_idx = c_x * (grid_config.K_y * grid_config.K_z) + c_y * grid_config.K_z + c_z;
    pos.atom_bin[i] = bin_idx;
    pos.bin_counts[bin_idx]++;
  }

  return pos;
}

/**
 * @struct FlatCellListData
 * @brief Storage for flat cell list index offsets and flattened indices.
 */
struct FlatCellListData {
  std::vector<size_t> bin_offsets;
  std::vector<size_t> bin_indices;
};

/**
 * @brief Constructs flattened cell list prefix sum and index buffer.
 */
FlatCellListData buildFlatCellList(const std::vector<size_t> &bin_counts, const std::vector<int> &atom_bin,
                                   size_t atom_count) {
  const size_t num_bins = bin_counts.size();
  std::vector<size_t> bin_offsets(num_bins + 1, 0);
  for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
    bin_offsets[bin_idx + 1] = bin_offsets[bin_idx] + bin_counts[bin_idx];
  }

  std::vector<size_t> bin_indices(atom_count);
  std::vector<size_t> insertion_cursors = bin_offsets;
  for (size_t i = 0; i < atom_count; ++i) {
    const int bin_idx = atom_bin[i];
    bin_indices[insertion_cursors[bin_idx]++] = i;
  }

  return {
      .bin_offsets = std::move(bin_offsets),
      .bin_indices = std::move(bin_indices),
  };
}

/**
 * @brief Merges thread-local distance calculations into output containers.
 */
void mergeThreadResults(const tbb::enumerable_thread_specific<ThreadLocalDistances> &ets, size_t num_elements,
                        DistanceTensor &out_distances, correlation::core::NeighborGraph &out_graph) {
  for (const auto &local_results : ets) {
    for (size_t i = 0; i < num_elements; ++i) {
      for (size_t j_idx = 0; j_idx < num_elements; ++j_idx) {
        out_distances[i][j_idx].insert(out_distances[i][j_idx].end(),
                                       local_results.distance_tensor_local[i][j_idx].begin(),
                                       local_results.distance_tensor_local[i][j_idx].end());
      }
    }
    for (const auto &bond : local_results.bonds) {
      out_graph.addDirectedEdge({
          .source = static_cast<correlation::core::AtomID>(bond.from),
          .target = static_cast<correlation::core::AtomID>(bond.to),
          .distance = bond.distance,
          .r_ij = bond.r_ij,
      });
    }
  }
}

} // namespace

void DistanceCalculator::compute(const correlation::core::Cell &cell, real_t cutoff_sq,
                                 const correlation::analysis::BondCutoffMatrix &bond_cutoffs,
                                 bool ignore_periodic_self_interactions, DistanceTensor &out_distances,
                                 correlation::core::NeighborGraph &out_graph) {
  validateInputs(cutoff_sq, bond_cutoffs);

#if defined(CORRELATION_USE_CUDA) || defined(CORRELATION_USE_HIP)
  if (tryComputeGpu(cell, cutoff_sq, bond_cutoffs, ignore_periodic_self_interactions, out_distances, out_graph)) {
    return;
  }
#endif

  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const size_t num_elements = cell.elements().size();
  const auto &lattice = cell.latticeVectors();

  const SearchGridConfig grid_config = buildSearchGridConfig(cell, cutoff_sq, ignore_periodic_self_interactions);
  const WrappedPositions wrapped = buildWrappedPositionsAndBins(cell, grid_config);
  const FlatCellListData cell_list_data = buildFlatCellList(wrapped.bin_counts, wrapped.atom_bin, atom_count);

  const SoACoordinates coords{
      .x = wrapped.wrapped_x.data(),
      .y = wrapped.wrapped_y.data(),
      .z = wrapped.wrapped_z.data(),
  };

  const FlatCellList cell_list{
      .offsets = cell_list_data.bin_offsets.data(),
      .indices = cell_list_data.bin_indices.data(),
  };

  tbb::enumerable_thread_specific<ThreadLocalDistances> ets(ThreadLocalConfig{
      .num_elements = num_elements,
      .atom_count = atom_count,
  });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, atom_count), [&](const tbb::blocked_range<size_t> &range) {
    ThreadLocalDistances &local_results = ets.local();
    for (size_t i = range.begin(); i != range.end(); ++i) {
      local_results.collectCandidates(i, coords, wrapped.atom_bin, cell_list, lattice, grid_config,
                                      ignore_periodic_self_interactions);
      local_results.computeDistances(i, atoms, coords, cutoff_sq, bond_cutoffs);
    }
  });

  mergeThreadResults(ets, num_elements, out_distances, out_graph);
}

} // namespace correlation::calculators
