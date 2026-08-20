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
enum class CandidateGatherMode : uint8_t {
  Partitioned,
  PrimaryCellOnly,
  PeriodicImageShift,
};

struct SearchGridConfig {
  int K_x;
  int K_y;
  int K_z;
  int max_dx;
  int max_dy;
  int max_dz;
  bool is_small_cell;
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
 * @struct BondCandidate
 * @brief Represents a candidate pair evaluated for bonding.
 */
struct BondCandidate {
  size_t atom_idx;
  size_t j_idx;
  size_t scratch_idx;
  int type_A;
  int type_B;
  real_t d_sq;
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
  size_t num_elements = 0;
  size_t atom_count = 0;
  real_t r_max = 0.0;
  real_t r_bin_width = 0.0;
  size_t num_bins = 0;
};

/**
 * @brief Evaluates whether a candidate pair should be filtered out.
 */
[[nodiscard]] constexpr bool shouldSkipCandidate(CandidateGatherMode mode, size_t atom_idx, size_t j_idx,
                                                 const correlation::math::Vector3<real_t> &disp,
                                                 bool ignore_periodic_self_interactions) noexcept {
  switch (mode) {
  case CandidateGatherMode::PrimaryCellOnly:
    return j_idx <= atom_idx;
  case CandidateGatherMode::Partitioned: {
    if (j_idx < atom_idx) {
      return true;
    }
    if (atom_idx == j_idx) {
      bool const is_zero_disp = (disp.x() == 0.0 && disp.y() == 0.0 && disp.z() == 0.0);
      return is_zero_disp || ignore_periodic_self_interactions;
    }
    return false;
  }
  case CandidateGatherMode::PeriodicImageShift:
    return atom_idx == j_idx && ignore_periodic_self_interactions;
  }
  return false;
}

/**
 * @struct ThreadLocalDistances
 * @brief Thread-local storage for parallel distance calculations and inline histogramming.
 *
 * Contains local raw histogram results, flat bond list, and pre-allocated scratch buffers
 * for SoA (Structure of Arrays) SIMD operations.
 */
struct ThreadLocalDistances {
  RawHistogramTensor local_histograms; ///< Local fixed-size raw histogram bins [e1][e2][bin].
  std::vector<ThreadLocalBond> bonds;  ///< Local flat bond/neighbor lists.
  real_t r_max = 0.0;
  real_t r_bin_width = 0.0;
  size_t num_bins = 0;

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
      : r_max(config.r_max), r_bin_width(config.r_bin_width), num_bins(config.num_bins) {
    if (config.num_bins > 0) {
      local_histograms.resize(config.num_elements, std::vector<std::vector<real_t>>(
                                                       config.num_elements, std::vector<real_t>(config.num_bins, 0.0)));
    }
    // Pre-reserve capacities to avoid resizing checks in hot loops
    soa_x.reserve(1024);
    soa_y.reserve(1024);
    soa_z.reserve(1024);
    candidate_j.reserve(1024);
    dsq_scratch.reserve(1024);
    bonds.reserve(2048);
  }

  /**
   * @brief Gathers candidate atoms using multi-image expansion for small crystal unit cells.
   */
  void collectSmallCellCandidates(size_t atom_idx, const SoACoordinates &coords, const FlatCellList &cell_list,
                                  const correlation::math::Matrix3<real_t> &lattice, SearchGridConfig grid_config,
                                  bool ignore_periodic_self_interactions, size_t &c_count);

  /**
   * @brief Gathers candidate atoms using 27-cell spatial stencils for partitioned cells.
   */
  void collectPartitionedCandidates(size_t atom_idx, const SoACoordinates &coords, const std::vector<int> &atom_bin,
                                    const FlatCellList &cell_list, const correlation::math::Matrix3<real_t> &lattice,
                                    SearchGridConfig grid_config, bool ignore_periodic_self_interactions,
                                    size_t &c_count);

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
                                const FlatCellList &cell_list, CandidateGatherMode mode, const SoACoordinates &coords,
                                bool ignore_periodic_self_interactions, size_t &c_count);

  /**
   * @brief Appends or updates a candidate in the thread-local scratch buffers.
   * @param shifted_pos Shifted 3D coordinate vector of the candidate atom.
   * @param j_idx Atom index of the candidate.
   * @param c_count Current candidate count accumulator reference.
   */
  void appendCandidate(const correlation::math::Vector3<real_t> &shifted_pos, size_t j_idx, size_t &c_count);

  /**
   * @brief Updates pair distribution histograms for a candidate pair.
   */
  void updateHistogram(const BondCandidate &cand);

  /**
   * @brief Evaluates and records bond connections for a candidate pair within cutoff ranges.
   */
  void updateBonds(const BondCandidate &cand, const SoACoordinates &coords,
                   const correlation::analysis::BondCutoffMatrix &bond_cutoffs,
                   const correlation::math::Vector3<real_t> &atom_pos);

  /**
   * @brief Computes SIMD pairwise distances and updates local structures.
   */
  void computeDistances(size_t atom_idx, const std::vector<correlation::core::Atom> &atoms,
                        const SoACoordinates &coords, real_t cutoff_sq,
                        const correlation::analysis::BondCutoffMatrix &bond_cutoffs);
};

void ThreadLocalDistances::collectSmallCellCandidates(size_t atom_idx, const SoACoordinates &coords,
                                                      const FlatCellList &cell_list,
                                                      const correlation::math::Matrix3<real_t> &lattice,
                                                      SearchGridConfig grid_config,
                                                      bool ignore_periodic_self_interactions, size_t &c_count) {
  // 1. Primary unit cell (zero displacement)
  collectCandidatesFromBin(atom_idx, {0.0, 0.0, 0.0}, 0, cell_list, CandidateGatherMode::PrimaryCellOnly, coords,
                           ignore_periodic_self_interactions, c_count);

  // 2. Lexicographically positive multi-image shifts (n > 0)
  for (int dx = -grid_config.max_dx; dx <= grid_config.max_dx; ++dx) {
    for (int dy = -grid_config.max_dy; dy <= grid_config.max_dy; ++dy) {
      for (int dz = -grid_config.max_dz; dz <= grid_config.max_dz; ++dz) {
        const bool is_positive_shift = (dx > 0) || (dx == 0 && dy > 0) || (dx == 0 && dy == 0 && dz > 0);
        if (!is_positive_shift) {
          continue;
        }

        real_t const disp_x =
            std::fma(static_cast<real_t>(dz), lattice[2].x(),
                     std::fma(static_cast<real_t>(dy), lattice[1].x(), static_cast<real_t>(dx) * lattice[0].x()));
        real_t const disp_y =
            std::fma(static_cast<real_t>(dz), lattice[2].y(),
                     std::fma(static_cast<real_t>(dy), lattice[1].y(), static_cast<real_t>(dx) * lattice[0].y()));
        real_t const disp_z =
            std::fma(static_cast<real_t>(dz), lattice[2].z(),
                     std::fma(static_cast<real_t>(dy), lattice[1].z(), static_cast<real_t>(dx) * lattice[0].z()));
        correlation::math::Vector3<real_t> const disp = {disp_x, disp_y, disp_z};

        collectCandidatesFromBin(atom_idx, disp, 0, cell_list, CandidateGatherMode::PeriodicImageShift, coords,
                                 ignore_periodic_self_interactions, c_count);
      }
    }
  }
}

void ThreadLocalDistances::collectPartitionedCandidates(size_t atom_idx, const SoACoordinates &coords,
                                                        const std::vector<int> &atom_bin, const FlatCellList &cell_list,
                                                        const correlation::math::Matrix3<real_t> &lattice,
                                                        SearchGridConfig grid_config,
                                                        bool ignore_periodic_self_interactions, size_t &c_count) {
  int const c_x = atom_bin[atom_idx] / (grid_config.K_y * grid_config.K_z);
  int const c_y = (atom_bin[atom_idx] / grid_config.K_z) % grid_config.K_y;
  int const c_z = atom_bin[atom_idx] % grid_config.K_z;

  for (int dx = -grid_config.max_dx; dx <= grid_config.max_dx; ++dx) {
    int const nx_bin = c_x + dx;
    int const shift_x = (nx_bin >= 0) ? (nx_bin / grid_config.K_x) : ((nx_bin - grid_config.K_x + 1) / grid_config.K_x);
    int const wrap_x = nx_bin - shift_x * grid_config.K_x;

    for (int dy = -grid_config.max_dy; dy <= grid_config.max_dy; ++dy) {
      int const ny_bin = c_y + dy;
      int const shift_y =
          (ny_bin >= 0) ? (ny_bin / grid_config.K_y) : ((ny_bin - grid_config.K_y + 1) / grid_config.K_y);
      int const wrap_y = ny_bin - shift_y * grid_config.K_y;

      for (int dz = -grid_config.max_dz; dz <= grid_config.max_dz; ++dz) {
        int const nz_bin = c_z + dz;
        int const shift_z =
            (nz_bin >= 0) ? (nz_bin / grid_config.K_z) : ((nz_bin - grid_config.K_z + 1) / grid_config.K_z);
        int const wrap_z = nz_bin - shift_z * grid_config.K_z;

        int const n_bin_idx = wrap_x * (grid_config.K_y * grid_config.K_z) + wrap_y * grid_config.K_z + wrap_z;

        real_t const disp_x = std::fma(
            static_cast<real_t>(shift_z), lattice[2].x(),
            std::fma(static_cast<real_t>(shift_y), lattice[1].x(), static_cast<real_t>(shift_x) * lattice[0].x()));
        real_t const disp_y = std::fma(
            static_cast<real_t>(shift_z), lattice[2].y(),
            std::fma(static_cast<real_t>(shift_y), lattice[1].y(), static_cast<real_t>(shift_x) * lattice[0].y()));
        real_t const disp_z = std::fma(
            static_cast<real_t>(shift_z), lattice[2].z(),
            std::fma(static_cast<real_t>(shift_y), lattice[1].z(), static_cast<real_t>(shift_x) * lattice[0].z()));
        correlation::math::Vector3<real_t> const disp = {disp_x, disp_y, disp_z};

        collectCandidatesFromBin(atom_idx, disp, n_bin_idx, cell_list, CandidateGatherMode::Partitioned, coords,
                                 ignore_periodic_self_interactions, c_count);
      }
    }
  }
}

void ThreadLocalDistances::collectCandidates(size_t atom_idx, const SoACoordinates &coords,
                                             const std::vector<int> &atom_bin, const FlatCellList &cell_list,
                                             const correlation::math::Matrix3<real_t> &lattice,
                                             SearchGridConfig grid_config, bool ignore_periodic_self_interactions) {
  size_t c_count = 0;
  if (grid_config.is_small_cell) {
    collectSmallCellCandidates(atom_idx, coords, cell_list, lattice, grid_config, ignore_periodic_self_interactions,
                               c_count);
  } else {
    collectPartitionedCandidates(atom_idx, coords, atom_bin, cell_list, lattice, grid_config,
                                 ignore_periodic_self_interactions, c_count);
  }
  candidate_count = c_count;
}

void ThreadLocalDistances::appendCandidate(const correlation::math::Vector3<real_t> &shifted_pos, size_t j_idx,
                                           size_t &c_count) {
  if (soa_x.size() <= c_count) {
    soa_x.push_back(shifted_pos.x());
    soa_y.push_back(shifted_pos.y());
    soa_z.push_back(shifted_pos.z());
    candidate_j.push_back(j_idx);
    dsq_scratch.push_back(0.0);
  } else {
    soa_x[c_count] = shifted_pos.x();
    soa_y[c_count] = shifted_pos.y();
    soa_z[c_count] = shifted_pos.z();
    candidate_j[c_count] = j_idx;
  }
  c_count++;
}

void ThreadLocalDistances::collectCandidatesFromBin(size_t atom_idx, const correlation::math::Vector3<real_t> &disp,
                                                    int n_bin_idx, const FlatCellList &cell_list,
                                                    CandidateGatherMode mode, const SoACoordinates &coords,
                                                    bool ignore_periodic_self_interactions, size_t &c_count) {

  size_t const start = cell_list.offsets[n_bin_idx];
  size_t const end = cell_list.offsets[n_bin_idx + 1];

  for (size_t offset = start; offset < end; ++offset) {
    size_t const j_idx = cell_list.indices[offset];
    if (shouldSkipCandidate(mode, atom_idx, j_idx, disp, ignore_periodic_self_interactions)) {
      continue;
    }

    correlation::math::Vector3<real_t> const shifted_pos{coords.x[j_idx] + disp.x(), coords.y[j_idx] + disp.y(),
                                                         coords.z[j_idx] + disp.z()};

    appendCandidate(shifted_pos, j_idx, c_count);
  }
}

void ThreadLocalDistances::updateHistogram(const BondCandidate &cand) {
  if (num_bins == 0 || r_bin_width <= 0.0) {
    return;
  }

  auto const dist = static_cast<real_t>(std::sqrt(cand.d_sq));
  if (dist >= r_max) {
    return;
  }

  auto const bin_idx = static_cast<size_t>(dist / r_bin_width);
  if (bin_idx < num_bins) {
    local_histograms[cand.type_A][cand.type_B][bin_idx] += 1.0;
    if (cand.atom_idx != cand.j_idx && cand.type_A != cand.type_B) {
      local_histograms[cand.type_B][cand.type_A][bin_idx] += 1.0;
    }
  }
}

void ThreadLocalDistances::updateBonds(const BondCandidate &cand, const SoACoordinates &coords,
                                       const correlation::analysis::BondCutoffMatrix &bond_cutoffs,
                                       const correlation::math::Vector3<real_t> &atom_pos) {
  real_t min_bond_dist_sq = 0.0;
  real_t max_bond_dist_sq = 0.0;
  if (static_cast<size_t>(cand.type_A) < bond_cutoffs.size() &&
      static_cast<size_t>(cand.type_B) < bond_cutoffs[cand.type_A].size()) {
    const auto &range = bond_cutoffs[cand.type_A][cand.type_B];
    min_bond_dist_sq = range.min_sq;
    max_bond_dist_sq = range.max_sq;
  }

  if (cand.d_sq < min_bond_dist_sq || cand.d_sq > max_bond_dist_sq) {
    return;
  }

  // Recover displacement vector coordinates precisely to minimize rounding/cancellation errors
  real_t const disp_x = soa_x[cand.scratch_idx] - coords.x[cand.j_idx];
  real_t const disp_y = soa_y[cand.scratch_idx] - coords.y[cand.j_idx];
  real_t const disp_z = soa_z[cand.scratch_idx] - coords.z[cand.j_idx];

  real_t const r_x = (coords.x[cand.j_idx] - atom_pos.x()) + disp_x;
  real_t const r_y = (coords.y[cand.j_idx] - atom_pos.y()) + disp_y;
  real_t const r_z = (coords.z[cand.j_idx] - atom_pos.z()) + disp_z;

  // Recompute consistent distance directly from the exact relative vector
  real_t const precise_dist = std::sqrt(r_x * r_x + r_y * r_y + r_z * r_z);
  correlation::math::Vector3<real_t> const r_ij = {r_x, r_y, r_z};

  bonds.push_back(ThreadLocalBond{
      .from = cand.atom_idx,
      .to = cand.j_idx,
      .distance = precise_dist,
      .r_ij = r_ij,
  });

  bonds.push_back(ThreadLocalBond{
      .from = cand.j_idx,
      .to = cand.atom_idx,
      .distance = precise_dist,
      .r_ij = -1.0 * r_ij,
  });
}

void ThreadLocalDistances::computeDistances(size_t atom_idx, const std::vector<correlation::core::Atom> &atoms,
                                            const SoACoordinates &coords, real_t cutoff_sq,
                                            const correlation::analysis::BondCutoffMatrix &bond_cutoffs) {

  size_t const c_count = candidate_count;
  if (c_count == 0) {
    return;
  }

  const auto &atom_A = atoms[atom_idx];
  const int type_A = atom_A.element_id();
  const real_t a_x = coords.x[atom_idx];
  const real_t a_y = coords.y[atom_idx];
  const real_t a_z = coords.z[atom_idx];
  correlation::math::Vector3<real_t> const atom_pos{a_x, a_y, a_z};

  const correlation::math::PositionBlock block{
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
    int const type_B = atoms[j_idx].element_id();

    const BondCandidate cand{
        .atom_idx = atom_idx,
        .j_idx = j_idx,
        .scratch_idx = k,
        .type_A = type_A,
        .type_B = type_B,
        .d_sq = d_sq,
    };

    updateHistogram(cand);
    updateBonds(cand, coords, bond_cutoffs, atom_pos);
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
                   RawHistogramTensor *out_histograms, DistanceCalculationConfig hist_config,
                   correlation::core::NeighborGraph &out_graph) {
  if (gpu::has_gpu_device()) {
    std::vector<std::vector<real_t>> max_cutoffs_sq(bond_cutoffs.size());
    for (size_t i = 0; i < bond_cutoffs.size(); ++i) {
      max_cutoffs_sq[i].resize(bond_cutoffs[i].size());
      for (size_t j = 0; j < bond_cutoffs[i].size(); ++j) {
        max_cutoffs_sq[i][j] = bond_cutoffs[i][j].max_sq;
      }
    }
    gpu::compute_distances_gpu(cell, cutoff_sq, max_cutoffs_sq, ignore_periodic_self_interactions, out_histograms,
                               hist_config, out_graph);
    return true;
  }
  return false;
}
#endif

/**
 * @brief Computes search grid dimensions and cutoff parameters.
 */
SearchGridConfig buildSearchGridConfig(const correlation::core::Cell &cell, real_t cutoff_sq) {
  const auto &lattice = cell.latticeVectors();
  const math::Vector3<real_t> lattice_a = lattice[0];
  const math::Vector3<real_t> lattice_b = lattice[1];
  const math::Vector3<real_t> lattice_c = lattice[2];

  const math::Vector3<real_t> cross_bc = math::cross(lattice_b, lattice_c);
  const math::Vector3<real_t> cross_ca = math::cross(lattice_c, lattice_a);
  const math::Vector3<real_t> cross_ab = math::cross(lattice_a, lattice_b);

  const real_t volume = cell.volume();
  const real_t width_x = volume / math::norm(cross_bc);
  const real_t width_y = volume / math::norm(cross_ca);
  const real_t width_z = volume / math::norm(cross_ab);

  const real_t cutoff = std::sqrt(cutoff_sq);

  const int K_x = std::max(1, static_cast<int>(std::floor(width_x / cutoff)));
  const int K_y = std::max(1, static_cast<int>(std::floor(width_y / cutoff)));
  const int K_z = std::max(1, static_cast<int>(std::floor(width_z / cutoff)));

  const bool is_small_cell = (K_x == 1 && K_y == 1 && K_z == 1);

  const int max_dx = (K_x == 1) ? std::max(1, static_cast<int>(std::ceil(cutoff / width_x))) : 1;
  const int max_dy = (K_y == 1) ? std::max(1, static_cast<int>(std::ceil(cutoff / width_y))) : 1;
  const int max_dz = (K_z == 1) ? std::max(1, static_cast<int>(std::ceil(cutoff / width_z))) : 1;

  return SearchGridConfig{
      .K_x = K_x,
      .K_y = K_y,
      .K_z = K_z,
      .max_dx = max_dx,
      .max_dy = max_dy,
      .max_dz = max_dz,
      .is_small_cell = is_small_cell,
  };
}

/**
 * @struct WrappedPositions
 * @brief Wrapped atom coordinates and spatial bin partitioning data.
 */
struct WrappedPositions {
  std::vector<real_t> wrapped_x;
  std::vector<real_t> wrapped_y;
  std::vector<real_t> wrapped_z;
  std::vector<int> atom_bin;
  std::vector<size_t> bin_counts;
};

/**
 * @brief Computes wrapped positions and bins for all atoms in the cell.
 */
WrappedPositions buildWrappedPositionsAndBins(const correlation::core::Cell &cell, SearchGridConfig grid_config) {
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
                        RawHistogramTensor *out_histograms, size_t num_bins,
                        correlation::core::NeighborGraph &out_graph) {
  for (const auto &local_results : ets) {
    if (out_histograms != nullptr && num_bins > 0) {
      for (size_t i = 0; i < num_elements; ++i) {
        for (size_t j_idx = 0; j_idx < num_elements; ++j_idx) {
          for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
            (*out_histograms)[i][j_idx][bin_idx] += local_results.local_histograms[i][j_idx][bin_idx];
          }
        }
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
                                 bool ignore_periodic_self_interactions, correlation::core::NeighborGraph &out_graph,
                                 RawHistogramTensor *out_histograms, DistanceCalculationConfig hist_config) {
  validateInputs(cutoff_sq, bond_cutoffs);

  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const size_t num_elements = cell.elements().size();
  const auto &lattice = cell.latticeVectors();

  if (out_histograms != nullptr && hist_config.num_bins > 0) {
    out_histograms->assign(
        num_elements, std::vector<std::vector<real_t>>(num_elements, std::vector<real_t>(hist_config.num_bins, 0.0)));
  }

#if defined(CORRELATION_USE_CUDA) || defined(CORRELATION_USE_HIP)
  if (tryComputeGpu(cell, cutoff_sq, bond_cutoffs, ignore_periodic_self_interactions, out_histograms, hist_config,
                    out_graph)) {
    return;
  }
#endif

  const SearchGridConfig grid_config = buildSearchGridConfig(cell, cutoff_sq);
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
      .r_max = hist_config.r_max,
      .r_bin_width = hist_config.r_bin_width,
      .num_bins = hist_config.num_bins,
  });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, atom_count), [&](const tbb::blocked_range<size_t> &range) {
    ThreadLocalDistances &local_results = ets.local();
    for (size_t i = range.begin(); i != range.end(); ++i) {
      local_results.collectCandidates(i, coords, wrapped.atom_bin, cell_list, lattice, grid_config,
                                      ignore_periodic_self_interactions);
      local_results.computeDistances(i, atoms, coords, cutoff_sq, bond_cutoffs);
    }
  });

  mergeThreadResults(ets, num_elements, out_histograms, hist_config.num_bins, out_graph);
}

} // namespace correlation::calculators
