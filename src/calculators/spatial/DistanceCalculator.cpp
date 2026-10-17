/**
 * @file DistanceCalculator.cpp
 * @brief Implementation of pairwise distance calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/spatial/DistanceCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"
#include "math/SIMDUtils.hpp"

#if defined(CORRELATION_USE_CUDA) || defined(CORRELATION_USE_HIP)
#include "calculators/gpu/GPUDistanceCalculator.hpp"
#endif

#include <atomic>
#include <cmath>
#include <span>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for_each.h>

namespace correlation::calculators {

void DistanceCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions & /*dists*/,
    const correlation::analysis::AnalysisSettings & /*settings*/) const {}

namespace {
enum class CandidateGatherMode : uint8_t {
  PARTITIONED,
  PrimaryCellOnly,
  PeriodicImageShift,
};

struct SearchGridConfig {
  int k_x;
  int k_y;
  int k_z;
  int max_dx;
  int max_dy;
  int max_dz;
  bool is_small_cell;
  std::vector<correlation::math::Vector3<real_t>> small_cell_shifts;
  std::array<int, 27> interior_bin_deltas;
};

/**
 * @struct SoACoordinates
 * @brief Structure of Arrays (SoA) view of Cartesian positions.
 */
struct SoACoordinates {
  const real_t *x;
  const real_t *y;
  const real_t *z;
};

/**
 * @struct BondCandidate
 * @brief Represents a candidate pair evaluated for bonding.
 */
struct BondCandidate {
  size_t atom_idx;
  size_t j_idx;
  size_t scratch_idx;
  int type_a;
  int type_b;
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
 * @struct FlatBondCutoff
 * @brief Flat lookup entry for species-pair bond distance cutoffs.
 */
struct FlatBondCutoff {
  real_t min_sq = 0.0;
  real_t max_sq = 0.0;
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
[[nodiscard]] constexpr bool shouldSkipCandidate(CandidateGatherMode mode, size_t atom_idx,
                                                 size_t j_idx,
                                                 const correlation::math::Vector3<real_t> &disp,
                                                 bool ignore_periodic_self_interactions) noexcept {
  switch (mode) {
  case CandidateGatherMode::PrimaryCellOnly:
    return j_idx <= atom_idx;
  case CandidateGatherMode::PARTITIONED: {
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
 */
struct ThreadLocalDistances {
  std::vector<real_t>
      flat_histograms; ///< Flat raw histogram bins [(e1 * num_elements + e2) * num_bins + bin].
  std::vector<ThreadLocalBond> bonds; ///< Local flat bond/neighbor lists.
  real_t r_max = 0.0;
  real_t r_max_sq = 0.0;
  real_t r_bin_width = 0.0;
  size_t num_bins = 0;
  size_t num_elements = 0;
  uint64_t last_epoch = 0;

  std::vector<real_t> soa_x;       ///< Scratch for candidate x-coordinates.
  std::vector<real_t> soa_y;       ///< Scratch for candidate y-coordinates.
  std::vector<real_t> soa_z;       ///< Scratch for candidate z-coordinates.
  std::vector<real_t> dsq_scratch; ///< Scratch for squared distance results.
  std::vector<size_t> candidate_j; ///< Scratch for candidate original atom indices.
  size_t candidate_count = 0;      ///< Count of candidates populated in scratch buffers.

  void growScratchBuffers(size_t new_cap) {
    soa_x.resize(new_cap);
    soa_y.resize(new_cap);
    soa_z.resize(new_cap);
    candidate_j.resize(new_cap);
    dsq_scratch.resize(new_cap);
  }

  void reset(ThreadLocalConfig config) {
    r_max = config.r_max;
    r_max_sq = config.r_max * config.r_max;
    r_bin_width = config.r_bin_width;
    num_bins = config.num_bins;
    num_elements = config.num_elements;
    candidate_count = 0;
    bonds.clear();

    const size_t total_bins = num_elements * num_elements * num_bins;
    if (flat_histograms.size() != total_bins) {
      flat_histograms.assign(total_bins, 0.0);
    } else {
      std::ranges::fill(flat_histograms, 0.0);
    }
  }

  explicit ThreadLocalDistances(ThreadLocalConfig config = {}) {
    reset(config);
    growScratchBuffers(1024);
    bonds.reserve(2048);
  }

  void collectSmallCellCandidates(size_t orig_atom_idx, const SoACoordinates &sorted_coords,
                                  const size_t *orig_indices, const size_t *bin_offsets,
                                  const SearchGridConfig &grid_config,
                                  bool ignore_periodic_self_interactions, size_t &c_count);

  void collectPartitionedInterior(size_t orig_atom_idx, int current_bin,
                                  const SoACoordinates &sorted_coords, const size_t *orig_indices,
                                  const size_t *bin_offsets, const SearchGridConfig &grid_config,
                                  bool ignore_periodic_self_interactions, size_t &c_count);

  void collectPartitionedBoundary(size_t orig_atom_idx, int current_bin,
                                  const SoACoordinates &sorted_coords, const size_t *orig_indices,
                                  const size_t *bin_offsets,
                                  const correlation::math::Matrix3<real_t> &lattice,
                                  const SearchGridConfig &grid_config,
                                  bool ignore_periodic_self_interactions, size_t &c_count);

  void collectPartitionedCandidates(size_t orig_atom_idx, int current_bin,
                                    const SoACoordinates &sorted_coords, const size_t *orig_indices,
                                    const size_t *bin_offsets,
                                    const correlation::math::Matrix3<real_t> &lattice,
                                    const SearchGridConfig &grid_config,
                                    bool ignore_periodic_self_interactions, size_t &c_count);

  void collectCandidates(size_t orig_atom_idx, int current_bin, const SoACoordinates &sorted_coords,
                         const size_t *orig_indices, const size_t *bin_offsets,
                         const correlation::math::Matrix3<real_t> &lattice,
                         const SearchGridConfig &grid_config,
                         bool ignore_periodic_self_interactions);

  void collectCandidatesFromBin(size_t orig_atom_idx,
                                const correlation::math::Vector3<real_t> &disp, size_t bin_start,
                                size_t bin_end, CandidateGatherMode mode,
                                const SoACoordinates &sorted_coords, const size_t *orig_indices,
                                bool ignore_periodic_self_interactions, size_t &c_count);

  void updateHistogram(const BondCandidate &cand);

  void updateBonds(const BondCandidate &cand, std::span<const FlatBondCutoff> flat_cutoffs,
                   const correlation::math::Vector3<real_t> &atom_pos);

  void computeDistances(size_t orig_atom_idx, const correlation::math::Vector3<real_t> &atom_pos,
                        const std::vector<correlation::core::Atom> &atoms, real_t cutoff_sq,
                        std::span<const FlatBondCutoff> flat_cutoffs);
};

void ThreadLocalDistances::collectSmallCellCandidates(
    size_t orig_atom_idx, const SoACoordinates &sorted_coords, const size_t *orig_indices,
    const size_t *bin_offsets, const SearchGridConfig &grid_config,
    bool ignore_periodic_self_interactions, size_t &c_count) {
  size_t const bin_start = bin_offsets[0];
  size_t const bin_end = bin_offsets[1];

  // 1. Primary unit cell (zero displacement)
  collectCandidatesFromBin(orig_atom_idx, {0.0, 0.0, 0.0}, bin_start, bin_end,
                           CandidateGatherMode::PrimaryCellOnly, sorted_coords, orig_indices,
                           ignore_periodic_self_interactions, c_count);

  // 2. Precomputed lexicographically positive multi-image shifts (n > 0)
  for (const auto &disp : grid_config.small_cell_shifts) {
    collectCandidatesFromBin(orig_atom_idx, disp, bin_start, bin_end,
                             CandidateGatherMode::PeriodicImageShift, sorted_coords, orig_indices,
                             ignore_periodic_self_interactions, c_count);
  }
}

constexpr void computePeriodicShiftAndBin(int coord, int delta, int dim_k, int &shift,
                                          int &wrap) noexcept {
  int const shifted_coord = coord + delta;
  shift = (shifted_coord >= 0) ? (shifted_coord / dim_k) : ((shifted_coord - dim_k + 1) / dim_k);
  wrap = shifted_coord - shift * dim_k;
}

void ThreadLocalDistances::collectPartitionedInterior(
    size_t orig_atom_idx, int current_bin, const SoACoordinates &sorted_coords,
    const size_t *orig_indices, const size_t *bin_offsets, const SearchGridConfig &grid_config,
    bool ignore_periodic_self_interactions, size_t &c_count) {
  const correlation::math::Vector3<real_t> zero_disp{0.0, 0.0, 0.0};
  for (int const delta : grid_config.interior_bin_deltas) {
    int const neighbor_bin = current_bin + delta;
    size_t const bin_start = bin_offsets[neighbor_bin];
    size_t const bin_end = bin_offsets[neighbor_bin + 1];
    collectCandidatesFromBin(orig_atom_idx, zero_disp, bin_start, bin_end,
                             CandidateGatherMode::PARTITIONED, sorted_coords, orig_indices,
                             ignore_periodic_self_interactions, c_count);
  }
}

void ThreadLocalDistances::collectPartitionedBoundary(
    size_t orig_atom_idx, int current_bin, const SoACoordinates &sorted_coords,
    const size_t *orig_indices, const size_t *bin_offsets,
    const correlation::math::Matrix3<real_t> &lattice, const SearchGridConfig &grid_config,
    bool ignore_periodic_self_interactions, size_t &c_count) {
  int const c_x = current_bin / (grid_config.k_y * grid_config.k_z);
  int const c_y = (current_bin / grid_config.k_z) % grid_config.k_y;
  int const c_z = current_bin % grid_config.k_z;

  for (int dx = -grid_config.max_dx; dx <= grid_config.max_dx; ++dx) {
    int shift_x = 0;
    int wrap_x = 0;
    computePeriodicShiftAndBin(c_x, dx, grid_config.k_x, shift_x, wrap_x);

    for (int dy = -grid_config.max_dy; dy <= grid_config.max_dy; ++dy) {
      int shift_y = 0;
      int wrap_y = 0;
      computePeriodicShiftAndBin(c_y, dy, grid_config.k_y, shift_y, wrap_y);

      for (int dz = -grid_config.max_dz; dz <= grid_config.max_dz; ++dz) {
        int shift_z = 0;
        int wrap_z = 0;
        computePeriodicShiftAndBin(c_z, dz, grid_config.k_z, shift_z, wrap_z);

        int const n_bin_idx =
            wrap_x * (grid_config.k_y * grid_config.k_z) + wrap_y * grid_config.k_z + wrap_z;

        real_t const disp_x = std::fma(static_cast<real_t>(shift_z), lattice[2].x(),
                                       std::fma(static_cast<real_t>(shift_y), lattice[1].x(),
                                                static_cast<real_t>(shift_x) * lattice[0].x()));
        real_t const disp_y = std::fma(static_cast<real_t>(shift_z), lattice[2].y(),
                                       std::fma(static_cast<real_t>(shift_y), lattice[1].y(),
                                                static_cast<real_t>(shift_x) * lattice[0].y()));
        real_t const disp_z = std::fma(static_cast<real_t>(shift_z), lattice[2].z(),
                                       std::fma(static_cast<real_t>(shift_y), lattice[1].z(),
                                                static_cast<real_t>(shift_x) * lattice[0].z()));
        correlation::math::Vector3<real_t> const disp = {disp_x, disp_y, disp_z};

        size_t const bin_start = bin_offsets[n_bin_idx];
        size_t const bin_end = bin_offsets[n_bin_idx + 1];

        collectCandidatesFromBin(orig_atom_idx, disp, bin_start, bin_end,
                                 CandidateGatherMode::PARTITIONED, sorted_coords, orig_indices,
                                 ignore_periodic_self_interactions, c_count);
      }
    }
  }
}

void ThreadLocalDistances::collectPartitionedCandidates(
    size_t orig_atom_idx, int current_bin, const SoACoordinates &sorted_coords,
    const size_t *orig_indices, const size_t *bin_offsets,
    const correlation::math::Matrix3<real_t> &lattice, const SearchGridConfig &grid_config,
    bool ignore_periodic_self_interactions, size_t &c_count) {
  int const c_x = current_bin / (grid_config.k_y * grid_config.k_z);
  int const c_y = (current_bin / grid_config.k_z) % grid_config.k_y;
  int const c_z = current_bin % grid_config.k_z;

  bool const is_interior = (c_x > 0 && c_x < grid_config.k_x - 1) &&
                           (c_y > 0 && c_y < grid_config.k_y - 1) &&
                           (c_z > 0 && c_z < grid_config.k_z - 1);

  if (is_interior) {
    collectPartitionedInterior(orig_atom_idx, current_bin, sorted_coords, orig_indices, bin_offsets,
                               grid_config, ignore_periodic_self_interactions, c_count);
  } else {
    collectPartitionedBoundary(orig_atom_idx, current_bin, sorted_coords, orig_indices, bin_offsets,
                               lattice, grid_config, ignore_periodic_self_interactions, c_count);
  }
}

void ThreadLocalDistances::collectCandidates(size_t orig_atom_idx, int current_bin,
                                             const SoACoordinates &sorted_coords,
                                             const size_t *orig_indices, const size_t *bin_offsets,
                                             const correlation::math::Matrix3<real_t> &lattice,
                                             const SearchGridConfig &grid_config,
                                             bool ignore_periodic_self_interactions) {
  size_t c_count = 0;
  if (grid_config.is_small_cell) {
    collectSmallCellCandidates(orig_atom_idx, sorted_coords, orig_indices, bin_offsets, grid_config,
                               ignore_periodic_self_interactions, c_count);
  } else {
    collectPartitionedCandidates(orig_atom_idx, current_bin, sorted_coords, orig_indices,
                                 bin_offsets, lattice, grid_config,
                                 ignore_periodic_self_interactions, c_count);
  }
  candidate_count = c_count;
}

void ThreadLocalDistances::collectCandidatesFromBin(
    size_t orig_atom_idx, const correlation::math::Vector3<real_t> &disp, size_t bin_start,
    size_t bin_end, CandidateGatherMode mode, const SoACoordinates &sorted_coords,
    const size_t *orig_indices, bool ignore_periodic_self_interactions, size_t &c_count) {

  for (size_t slot = bin_start; slot < bin_end; ++slot) {
    size_t const j_idx = orig_indices[slot];
    if (shouldSkipCandidate(mode, orig_atom_idx, j_idx, disp, ignore_periodic_self_interactions)) {
      continue;
    }

    if (c_count >= soa_x.size()) {
      growScratchBuffers(soa_x.size() * 2);
    }

    soa_x[c_count] = sorted_coords.x[slot] + disp.x();
    soa_y[c_count] = sorted_coords.y[slot] + disp.y();
    soa_z[c_count] = sorted_coords.z[slot] + disp.z();
    candidate_j[c_count] = j_idx;
    c_count++;
  }
}

void ThreadLocalDistances::updateHistogram(const BondCandidate &cand) {
  if (num_bins == 0 || r_bin_width <= 0.0 || cand.d_sq >= r_max_sq) {
    return;
  }

  auto const dist = static_cast<real_t>(std::sqrt(cand.d_sq));
  auto const bin_idx = static_cast<size_t>(dist / r_bin_width);
  if (bin_idx < num_bins) {
    const size_t offset_ab =
        (static_cast<size_t>(cand.type_a) * num_elements + static_cast<size_t>(cand.type_b)) *
            num_bins +
        bin_idx;
    flat_histograms[offset_ab] += 1.0;
    if (cand.atom_idx != cand.j_idx && cand.type_a != cand.type_b) {
      const size_t offset_ba =
          (static_cast<size_t>(cand.type_b) * num_elements + static_cast<size_t>(cand.type_a)) *
              num_bins +
          bin_idx;
      flat_histograms[offset_ba] += 1.0;
    }
  }
}

void ThreadLocalDistances::updateBonds(const BondCandidate &cand,
                                       std::span<const FlatBondCutoff> flat_cutoffs,
                                       const correlation::math::Vector3<real_t> &atom_pos) {
  const size_t cutoff_idx =
      static_cast<size_t>(cand.type_a) * num_elements + static_cast<size_t>(cand.type_b);
  if (cutoff_idx >= flat_cutoffs.size()) {
    return;
  }
  const auto &range = flat_cutoffs[cutoff_idx];
  if (cand.d_sq < range.min_sq || cand.d_sq > range.max_sq) {
    return;
  }

  real_t const r_x = soa_x[cand.scratch_idx] - atom_pos.x();
  real_t const r_y = soa_y[cand.scratch_idx] - atom_pos.y();
  real_t const r_z = soa_z[cand.scratch_idx] - atom_pos.z();

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

void ThreadLocalDistances::computeDistances(size_t orig_atom_idx,
                                            const correlation::math::Vector3<real_t> &atom_pos,
                                            const std::vector<correlation::core::Atom> &atoms,
                                            real_t cutoff_sq,
                                            std::span<const FlatBondCutoff> flat_cutoffs) {

  size_t const c_count = candidate_count;
  if (c_count == 0) {
    return;
  }

  const int type_a = atoms[orig_atom_idx].elementId();
  const real_t a_x = atom_pos.x();
  const real_t a_y = atom_pos.y();
  const real_t a_z = atom_pos.z();

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
    int const type_b = atoms[j_idx].elementId();

    const BondCandidate cand{
        .atom_idx = orig_atom_idx,
        .j_idx = j_idx,
        .scratch_idx = k,
        .type_a = type_a,
        .type_b = type_b,
        .d_sq = d_sq,
    };

    updateHistogram(cand);
    updateBonds(cand, flat_cutoffs, atom_pos);
  }
}

[[nodiscard]] std::vector<FlatBondCutoff>
buildFlatBondCutoffs(const correlation::analysis::BondCutoffMatrix &bond_cutoffs,
                     size_t num_elements) {
  std::vector<FlatBondCutoff> flat_cutoffs(num_elements * num_elements, FlatBondCutoff{});
  for (size_t i = 0; i < num_elements && i < bond_cutoffs.size(); ++i) {
    for (size_t j = 0; j < num_elements && j < bond_cutoffs[i].size(); ++j) {
      const auto &range = bond_cutoffs(i, j);
      flat_cutoffs[i * num_elements + j] = FlatBondCutoff{
          .min_sq = range.min_sq,
          .max_sq = range.max_sq,
      };
    }
  }
  return flat_cutoffs;
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
        throw std::invalid_argument(
            "Minimum bond cutoff cannot be greater than maximum bond cutoff.");
      }
      constexpr auto K_EPSILON = static_cast<real_t>(1e-4);
      if (range.max_sq > cutoff_sq + K_EPSILON) {
        throw std::invalid_argument(
            "Maximum bond cutoff squared cannot exceed search cutoff squared.");
      }
    }
  }
}
struct ShiftBounds {
  int max_dx;
  int max_dy;
  int max_dz;
  real_t width_x;
  real_t width_y;
  real_t width_z;
  real_t cutoff;
};

std::vector<correlation::math::Vector3<real_t>>
computeSmallCellShifts(const correlation::math::Matrix3<real_t> &lattice,
                       const ShiftBounds &bounds) {
  std::vector<correlation::math::Vector3<real_t>> shifts;
  for (int dx = -bounds.max_dx; dx <= bounds.max_dx; ++dx) {
    for (int dy = -bounds.max_dy; dy <= bounds.max_dy; ++dy) {
      for (int dz = -bounds.max_dz; dz <= bounds.max_dz; ++dz) {
        const bool is_positive_shift =
            (dx > 0) || (dx == 0 && dy > 0) || (dx == 0 && dy == 0 && dz > 0);
        if (!is_positive_shift) {
          continue;
        }

        const real_t gap_x = static_cast<real_t>(std::max(0, std::abs(dx) - 1)) * bounds.width_x;
        const real_t gap_y = static_cast<real_t>(std::max(0, std::abs(dy) - 1)) * bounds.width_y;
        const real_t gap_z = static_cast<real_t>(std::max(0, std::abs(dz) - 1)) * bounds.width_z;
        if (std::max({gap_x, gap_y, gap_z}) > bounds.cutoff) {
          continue;
        }

        const real_t disp_x = std::fma(static_cast<real_t>(dz), lattice[2].x(),
                                       std::fma(static_cast<real_t>(dy), lattice[1].x(),
                                                static_cast<real_t>(dx) * lattice[0].x()));
        const real_t disp_y = std::fma(static_cast<real_t>(dz), lattice[2].y(),
                                       std::fma(static_cast<real_t>(dy), lattice[1].y(),
                                                static_cast<real_t>(dx) * lattice[0].y()));
        const real_t disp_z = std::fma(static_cast<real_t>(dz), lattice[2].z(),
                                       std::fma(static_cast<real_t>(dy), lattice[1].z(),
                                                static_cast<real_t>(dx) * lattice[0].z()));
        shifts.emplace_back(disp_x, disp_y, disp_z);
      }
    }
  }
  return shifts;
}

/**
 * @brief Computes search grid dimensions and cutoff parameters with numerical degeneracy guards.
 */
SearchGridConfig buildSearchGridConfig(const correlation::core::Cell &cell, real_t cutoff_sq) {
  constexpr auto K_ZERO = static_cast<real_t>(0.0);
  const real_t volume = std::abs(cell.volume());

  if (cutoff_sq <= K_ZERO || volume <= K_ZERO) {
    return SearchGridConfig{
        .k_x = 1,
        .k_y = 1,
        .k_z = 1,
        .max_dx = 1,
        .max_dy = 1,
        .max_dz = 1,
        .is_small_cell = true,
        .small_cell_shifts = {},
        .interior_bin_deltas = {},
    };
  }

  const auto &lattice = cell.latticeVectors();
  const math::Vector3<real_t> lattice_a = lattice[0];
  const math::Vector3<real_t> lattice_b = lattice[1];
  const math::Vector3<real_t> lattice_c = lattice[2];

  const math::Vector3<real_t> cross_bc = math::cross(lattice_b, lattice_c);
  const math::Vector3<real_t> cross_ca = math::cross(lattice_c, lattice_a);
  const math::Vector3<real_t> cross_ab = math::cross(lattice_a, lattice_b);

  const real_t norm_bc = math::norm(cross_bc);
  const real_t norm_ca = math::norm(cross_ca);
  const real_t norm_ab = math::norm(cross_ab);

  if (norm_bc <= K_ZERO || norm_ca <= K_ZERO || norm_ab <= K_ZERO) {
    return SearchGridConfig{
        .k_x = 1,
        .k_y = 1,
        .k_z = 1,
        .max_dx = 1,
        .max_dy = 1,
        .max_dz = 1,
        .is_small_cell = true,
        .small_cell_shifts = {},
        .interior_bin_deltas = {},
    };
  }

  const real_t width_x = volume / norm_bc;
  const real_t width_y = volume / norm_ca;
  const real_t width_z = volume / norm_ab;

  const real_t cutoff = std::sqrt(cutoff_sq);

  const int k_x = std::max(1, static_cast<int>(std::floor(width_x / cutoff)));
  const int k_y = std::max(1, static_cast<int>(std::floor(width_y / cutoff)));
  const int k_z = std::max(1, static_cast<int>(std::floor(width_z / cutoff)));

  const bool is_small_cell = (k_x == 1 && k_y == 1 && k_z == 1);

  const int max_dx = (k_x == 1) ? std::max(1, static_cast<int>(std::ceil(cutoff / width_x))) : 1;
  const int max_dy = (k_y == 1) ? std::max(1, static_cast<int>(std::ceil(cutoff / width_y))) : 1;
  const int max_dz = (k_z == 1) ? std::max(1, static_cast<int>(std::ceil(cutoff / width_z))) : 1;

  std::vector<correlation::math::Vector3<real_t>> small_cell_shifts;
  if (is_small_cell) {
    small_cell_shifts = computeSmallCellShifts(lattice, ShiftBounds{
                                                            .max_dx = max_dx,
                                                            .max_dy = max_dy,
                                                            .max_dz = max_dz,
                                                            .width_x = width_x,
                                                            .width_y = width_y,
                                                            .width_z = width_z,
                                                            .cutoff = cutoff,
                                                        });
  }

  std::array<int, 27> interior_bin_deltas{};
  if (!is_small_cell) {
    size_t d_idx = 0;
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          interior_bin_deltas.at(d_idx++) = dx * (k_y * k_z) + dy * k_z + dz;
        }
      }
    }
  }

  return SearchGridConfig{
      .k_x = k_x,
      .k_y = k_y,
      .k_z = k_z,
      .max_dx = max_dx,
      .max_dy = max_dy,
      .max_dz = max_dz,
      .is_small_cell = is_small_cell,
      .small_cell_shifts = std::move(small_cell_shifts),
      .interior_bin_deltas = interior_bin_deltas,
  };
}

#if defined(CORRELATION_USE_CUDA) || defined(CORRELATION_USE_HIP)
/**
 * @brief Dispatches calculation to GPU if a compatible device is present and system is eligible.
 */
bool tryComputeGpu(const correlation::core::Cell &cell, real_t cutoff_sq,
                   const correlation::analysis::BondCutoffMatrix &bond_cutoffs,
                   bool ignore_periodic_self_interactions, RawHistogramTensor *out_histograms,
                   DistanceCalculationConfig hist_config,
                   correlation::core::NeighborGraph &out_graph,
                   const SearchGridConfig &grid_config) {
  if (grid_config.is_small_cell) {
    return false;
  }
  for (const auto &row : bond_cutoffs) {
    for (const auto &b : row) {
      if (b.min_sq > 0.0) {
        return false;
      }
    }
  }
  if (gpu::has_gpu_device()) {
    std::vector<std::vector<real_t>> max_cutoffs_sq(bond_cutoffs.size());
    for (size_t i = 0; i < bond_cutoffs.size(); ++i) {
      max_cutoffs_sq[i].resize(bond_cutoffs[i].size());
      for (size_t j = 0; j < bond_cutoffs[i].size(); ++j) {
        max_cutoffs_sq[i][j] = bond_cutoffs[i][j].max_sq;
      }
    }
    gpu::compute_distances_gpu(cell, cutoff_sq, max_cutoffs_sq, ignore_periodic_self_interactions,
                               out_histograms, hist_config, out_graph);
    return true;
  }
  return false;
}
#endif

/**
 * @brief Computes wrapped Cartesian coordinates and 3D grid bin for a single atom.
 */
inline void wrapAtomCoordinates(const correlation::math::Vector3<real_t> &pos,
                                const correlation::math::Matrix3<real_t> &inv_lattice,
                                const correlation::math::Matrix3<real_t> &lattice,
                                const SearchGridConfig &grid_config, real_t &wrapped_x,
                                real_t &wrapped_y, real_t &wrapped_z, int &bin_idx) noexcept {
  const correlation::math::Vector3<real_t> frac = inv_lattice * pos;
  const real_t f_x = frac.x() - std::floor(frac.x());
  const real_t f_y = frac.y() - std::floor(frac.y());
  const real_t f_z = frac.z() - std::floor(frac.z());

  wrapped_x = std::fma(f_z, lattice[2].x(), std::fma(f_y, lattice[1].x(), f_x * lattice[0].x()));
  wrapped_y = std::fma(f_z, lattice[2].y(), std::fma(f_y, lattice[1].y(), f_x * lattice[0].y()));
  wrapped_z = std::fma(f_z, lattice[2].z(), std::fma(f_y, lattice[1].z(), f_x * lattice[0].z()));

  const int c_x =
      std::clamp(static_cast<int>(std::floor(f_x * static_cast<real_t>(grid_config.k_x))), 0,
                 grid_config.k_x - 1);
  const int c_y =
      std::clamp(static_cast<int>(std::floor(f_y * static_cast<real_t>(grid_config.k_y))), 0,
                 grid_config.k_y - 1);
  const int c_z =
      std::clamp(static_cast<int>(std::floor(f_z * static_cast<real_t>(grid_config.k_z))), 0,
                 grid_config.k_z - 1);

  bin_idx = c_x * (grid_config.k_y * grid_config.k_z) + c_y * grid_config.k_z + c_z;
}

/**
 * @struct BinSortedPositions
 * @brief Contiguous SoA coordinates, atom mapping, and bin slices sorted by spatial bin.
 */
struct BinSortedPositions {
  std::vector<real_t> x_sorted;     ///< Contiguous bin-sorted x coordinates.
  std::vector<real_t> y_sorted;     ///< Contiguous bin-sorted y coordinates.
  std::vector<real_t> z_sorted;     ///< Contiguous bin-sorted z coordinates.
  std::vector<size_t> orig_indices; ///< Maps sorted slot to original atom index.
  std::vector<int> sorted_atom_bin; ///< Spatial bin index for each sorted atom slot.
  std::vector<size_t> bin_offsets;  ///< Slices [bin_offsets[b], bin_offsets[b+1]) per bin.
};

/**
 * @brief Constructs bin-sorted SoA positions and bin slice offsets.
 */
BinSortedPositions buildBinSortedPositions(const correlation::core::Cell &cell,
                                           const SearchGridConfig &grid_config) {
  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const auto &lattice = cell.latticeVectors();
  const auto &inv_lattice = cell.inverseLatticeVectors();
  const auto num_bins = static_cast<size_t>(grid_config.k_x) *
                        static_cast<size_t>(grid_config.k_y) * static_cast<size_t>(grid_config.k_z);

  std::vector<real_t> wrapped_x(atom_count);
  std::vector<real_t> wrapped_y(atom_count);
  std::vector<real_t> wrapped_z(atom_count);
  std::vector<int> atom_bin(atom_count);
  std::vector<size_t> bin_counts(num_bins, 0);

  if (atom_count < 2048) {
    for (size_t i = 0; i < atom_count; ++i) {
      wrapAtomCoordinates(atoms[i].position(), inv_lattice, lattice, grid_config, wrapped_x[i],
                          wrapped_y[i], wrapped_z[i], atom_bin[i]);
      bin_counts[atom_bin[i]]++;
    }
  } else {
    tbb::enumerable_thread_specific<std::vector<size_t>> local_counts(
        [num_bins] { return std::vector<size_t>(num_bins, 0); });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, atom_count), [&](const tbb::blocked_range<size_t> &range) {
          auto &counts = local_counts.local();
          for (size_t i = range.begin(); i != range.end(); ++i) {
            wrapAtomCoordinates(atoms[i].position(), inv_lattice, lattice, grid_config,
                                wrapped_x[i], wrapped_y[i], wrapped_z[i], atom_bin[i]);
            counts[atom_bin[i]]++;
          }
        });

    for (const auto &local_thread_counts : local_counts) {
      for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
        bin_counts[bin_idx] += local_thread_counts[bin_idx];
      }
    }
  }

  std::vector<size_t> bin_offsets(num_bins + 1, 0);
  for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
    bin_offsets[bin_idx + 1] = bin_offsets[bin_idx] + bin_counts[bin_idx];
  }

  std::vector<real_t> x_sorted(atom_count);
  std::vector<real_t> y_sorted(atom_count);
  std::vector<real_t> z_sorted(atom_count);
  std::vector<size_t> orig_indices(atom_count);
  std::vector<int> sorted_atom_bin(atom_count);

  std::vector<size_t> cursors = bin_offsets;
  for (size_t i = 0; i < atom_count; ++i) {
    const int bin_idx = atom_bin[i];
    const size_t slot = cursors[bin_idx]++;
    x_sorted[slot] = wrapped_x[i];
    y_sorted[slot] = wrapped_y[i];
    z_sorted[slot] = wrapped_z[i];
    orig_indices[slot] = i;
    sorted_atom_bin[slot] = bin_idx;
  }

  return BinSortedPositions{
      .x_sorted = std::move(x_sorted),
      .y_sorted = std::move(y_sorted),
      .z_sorted = std::move(z_sorted),
      .orig_indices = std::move(orig_indices),
      .sorted_atom_bin = std::move(sorted_atom_bin),
      .bin_offsets = std::move(bin_offsets),
  };
}

/**
 * @brief Merges thread-local distance histograms into the output tensor.
 */
void mergeHistograms(const tbb::enumerable_thread_specific<ThreadLocalDistances> &ets,
                     uint64_t current_epoch, size_t num_elements, size_t num_bins,
                     RawHistogramTensor &out_histograms) {
  for (const auto &local_results : ets) {
    if (local_results.last_epoch != current_epoch || local_results.flat_histograms.empty()) {
      continue;
    }
    for (size_t i = 0; i < num_elements; ++i) {
      for (size_t j_idx = 0; j_idx < num_elements; ++j_idx) {
        const size_t base_offset = (i * num_elements + j_idx) * num_bins;
        for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
          out_histograms[i][j_idx][bin_idx] += local_results.flat_histograms[base_offset + bin_idx];
        }
      }
    }
  }
}

void countRowBonds(const tbb::enumerable_thread_specific<ThreadLocalDistances> &ets,
                   uint64_t current_epoch, size_t atom_count, std::span<size_t> row_counts,
                   size_t &total_bonds) {
  for (const auto &local_results : ets) {
    if (local_results.last_epoch != current_epoch) {
      continue;
    }
    total_bonds += local_results.bonds.size();
    for (const auto &bond : local_results.bonds) {
      if (bond.from < atom_count) {
        row_counts[bond.from]++;
      }
    }
  }
}

void populateCSRBonds(const tbb::enumerable_thread_specific<ThreadLocalDistances> &ets,
                      uint64_t current_epoch, size_t atom_count, std::span<size_t> write_pos,
                      std::span<correlation::core::Neighbor> neighbors) {
  for (const auto &local_results : ets) {
    if (local_results.last_epoch != current_epoch) {
      continue;
    }
    for (const auto &bond : local_results.bonds) {
      if (bond.from < atom_count) {
        const size_t slot = write_pos[bond.from]++;
        neighbors[slot] = correlation::core::Neighbor{
            .index = static_cast<correlation::core::AtomID>(bond.to),
            .distance = bond.distance,
            .r_ij = bond.r_ij,
        };
      }
    }
  }
}

void sortCSRRows(size_t atom_count, std::span<const size_t> offsets,
                 std::span<correlation::core::Neighbor> neighbors,
                 std::span<correlation::core::AtomID> targets) {
  tbb::parallel_for(tbb::blocked_range<size_t>(0, atom_count),
                    [&](const tbb::blocked_range<size_t> &range) {
                      for (size_t node_idx = range.begin(); node_idx != range.end(); ++node_idx) {
                        const size_t start = offsets[node_idx];
                        const size_t end = offsets[node_idx + 1];
                        std::sort(neighbors.begin() + static_cast<std::ptrdiff_t>(start),
                                  neighbors.begin() + static_cast<std::ptrdiff_t>(end),
                                  [](const correlation::core::Neighbor &nbr_a,
                                     const correlation::core::Neighbor &nbr_b) {
                                    if (nbr_a.index != nbr_b.index) {
                                      return nbr_a.index < nbr_b.index;
                                    }
                                    if (nbr_a.distance != nbr_b.distance) {
                                      return nbr_a.distance < nbr_b.distance;
                                    }
                                    if (nbr_a.r_ij.x() != nbr_b.r_ij.x()) {
                                      return nbr_a.r_ij.x() < nbr_b.r_ij.x();
                                    }
                                    if (nbr_a.r_ij.y() != nbr_b.r_ij.y()) {
                                      return nbr_a.r_ij.y() < nbr_b.r_ij.y();
                                    }
                                    return nbr_a.r_ij.z() < nbr_b.r_ij.z();
                                  });
                        for (size_t k = start; k < end; ++k) {
                          targets[k] = neighbors[k].index;
                        }
                      }
                    });
}

/**
 * @brief Builds a deterministic CSR NeighborGraph from thread-local bonds.
 */
void buildDeterministicCSR(const tbb::enumerable_thread_specific<ThreadLocalDistances> &ets,
                           uint64_t current_epoch, size_t atom_count,
                           correlation::core::NeighborGraph &out_graph) {
  std::vector<size_t> row_counts(atom_count, 0);
  size_t total_bonds = 0;
  countRowBonds(ets, current_epoch, atom_count, row_counts, total_bonds);

  std::vector<size_t> offsets(atom_count + 1, 0);
  for (size_t i = 0; i < atom_count; ++i) {
    offsets[i + 1] = offsets[i] + row_counts[i];
  }

  std::vector<correlation::core::Neighbor> neighbors(total_bonds);
  std::vector<correlation::core::AtomID> targets(total_bonds);
  std::vector<size_t> write_pos = offsets;
  populateCSRBonds(ets, current_epoch, atom_count, write_pos, neighbors);
  sortCSRRows(atom_count, offsets, neighbors, targets);

  out_graph.assignCSR(atom_count, std::move(offsets), std::move(neighbors), std::move(targets));
}

/**
 * @brief Merges thread-local distance calculations into output containers.
 */
void mergeThreadResults(const tbb::enumerable_thread_specific<ThreadLocalDistances> &ets,
                        uint64_t current_epoch, size_t num_elements,
                        RawHistogramTensor *out_histograms, size_t num_bins,
                        correlation::core::NeighborGraph &out_graph, size_t atom_count) {
  if (out_histograms != nullptr && num_bins > 0) {
    mergeHistograms(ets, current_epoch, num_elements, num_bins, *out_histograms);
  }
  buildDeterministicCSR(ets, current_epoch, atom_count, out_graph);
}

uint64_t getNextQueryEpoch() {
  static std::atomic<uint64_t> s_distance_query_epoch{1};
  return s_distance_query_epoch.fetch_add(1, std::memory_order_relaxed);
}

tbb::enumerable_thread_specific<ThreadLocalDistances> &getDistanceCalculatorWorkspace() {
  static tbb::enumerable_thread_specific<ThreadLocalDistances> workspace(
      [] { return ThreadLocalDistances(); });
  return workspace;
}

} // namespace

void DistanceCalculator::compute(const correlation::core::Cell &cell, real_t cutoff_sq,
                                 const correlation::analysis::BondCutoffMatrix &bond_cutoffs,
                                 bool ignore_periodic_self_interactions,
                                 correlation::core::NeighborGraph &out_graph,
                                 RawHistogramTensor *out_histograms,
                                 DistanceCalculationConfig hist_config) {
  validateInputs(cutoff_sq, bond_cutoffs);

  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const size_t num_elements = cell.elements().size();
  const auto &lattice = cell.latticeVectors();

  if (out_graph.nodeCount() != atom_count) {
    out_graph = correlation::core::NeighborGraph(atom_count);
  }

  if (out_histograms != nullptr && hist_config.num_bins > 0) {
    out_histograms->assign(num_elements,
                           std::vector<std::vector<real_t>>(
                               num_elements, std::vector<real_t>(hist_config.num_bins, 0.0)));
  }

  const SearchGridConfig grid_config = buildSearchGridConfig(cell, cutoff_sq);

#if defined(CORRELATION_USE_CUDA) || defined(CORRELATION_USE_HIP)
  if (tryComputeGpu(cell, cutoff_sq, bond_cutoffs, ignore_periodic_self_interactions,
                    out_histograms, hist_config, out_graph, grid_config)) {
    return;
  }
#endif
  const BinSortedPositions bin_sorted = buildBinSortedPositions(cell, grid_config);

  const SoACoordinates sorted_coords{
      .x = bin_sorted.x_sorted.data(),
      .y = bin_sorted.y_sorted.data(),
      .z = bin_sorted.z_sorted.data(),
  };

  const auto flat_cutoffs = buildFlatBondCutoffs(bond_cutoffs, num_elements);
  const uint64_t epoch = getNextQueryEpoch();
  auto &ets = getDistanceCalculatorWorkspace();

  const ThreadLocalConfig tl_config{
      .num_elements = num_elements,
      .atom_count = atom_count,
      .r_max = hist_config.r_max,
      .r_bin_width = hist_config.r_bin_width,
      .num_bins = hist_config.num_bins,
  };

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, atom_count), [&](const tbb::blocked_range<size_t> &range) {
        ThreadLocalDistances &local_results = ets.local();
        if (local_results.last_epoch != epoch) {
          local_results.reset(tl_config);
          local_results.last_epoch = epoch;
        }

        for (size_t slot = range.begin(); slot != range.end(); ++slot) {
          size_t const orig_idx = bin_sorted.orig_indices[slot];
          int const current_bin = bin_sorted.sorted_atom_bin[slot];
          correlation::math::Vector3<real_t> const atom_pos{
              bin_sorted.x_sorted[slot], bin_sorted.y_sorted[slot], bin_sorted.z_sorted[slot]};

          local_results.collectCandidates(orig_idx, current_bin, sorted_coords,
                                          bin_sorted.orig_indices.data(),
                                          bin_sorted.bin_offsets.data(), lattice, grid_config,
                                          ignore_periodic_self_interactions);
          local_results.computeDistances(orig_idx, atom_pos, atoms, cutoff_sq, flat_cutoffs);
        }
      });

  mergeThreadResults(ets, epoch, num_elements, out_histograms, hist_config.num_bins, out_graph,
                     atom_count);
}

} // namespace correlation::calculators
