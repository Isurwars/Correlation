/**
 * @file CellListTraversal.hpp
 * @brief High-performance bin-sorted cell-list traversal and pairwise visitor engine.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Cell.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"
#include "math/SIMDUtils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace correlation::calculators::detail {

enum class CandidateGatherMode : uint8_t {
  PARTITIONED,
  PrimaryCellOnly,
  PeriodicImageShift,
};

struct SearchGridConfig {
  int k_x{1};
  int k_y{1};
  int k_z{1};
  int max_dx{1};
  int max_dy{1};
  int max_dz{1};
  bool is_small_cell{false};
  std::vector<correlation::math::Vector3<real_t>> small_cell_shifts;
  std::array<int, 27> interior_bin_deltas{};
};

struct SoACoordinates {
  const real_t *x{nullptr};
  const real_t *y{nullptr};
  const real_t *z{nullptr};
};

struct BinSortedPositions {
  std::vector<real_t> x_sorted;
  std::vector<real_t> y_sorted;
  std::vector<real_t> z_sorted;
  std::vector<size_t> orig_indices;
  std::vector<int> sorted_atom_bin;
  std::vector<size_t> bin_offsets;
};

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
      const bool is_zero_disp = (disp.x() == 0.0 && disp.y() == 0.0 && disp.z() == 0.0);
      return is_zero_disp || ignore_periodic_self_interactions;
    }
    return false;
  }
  case CandidateGatherMode::PeriodicImageShift:
    return atom_idx == j_idx && ignore_periodic_self_interactions;
  }
  return false;
}

constexpr void computePeriodicShiftAndBin(int coord, int delta, int dim_k, int &shift,
                                          int &wrap) noexcept {
  const int shifted_coord = coord + delta;
  shift = (shifted_coord >= 0) ? (shifted_coord / dim_k) : ((shifted_coord - dim_k + 1) / dim_k);
  wrap = shifted_coord - shift * dim_k;
}

struct CandidateGatherer {
  std::vector<real_t> soa_x;
  std::vector<real_t> soa_y;
  std::vector<real_t> soa_z;
  std::vector<real_t> dsq_scratch;
  std::vector<size_t> candidate_j;
  size_t candidate_count{0};

  void growScratchBuffers(size_t new_cap) {
    soa_x.resize(new_cap);
    soa_y.resize(new_cap);
    soa_z.resize(new_cap);
    candidate_j.resize(new_cap);
    dsq_scratch.resize(new_cap);
  }

  void collectCandidatesFromBin(size_t orig_atom_idx,
                                const correlation::math::Vector3<real_t> &disp, size_t bin_start,
                                size_t bin_end, CandidateGatherMode mode,
                                const SoACoordinates &sorted_coords, const size_t *orig_indices,
                                bool ignore_periodic_self_interactions) {
    for (size_t slot = bin_start; slot < bin_end; ++slot) {
      const size_t j_idx = orig_indices[slot];
      if (shouldSkipCandidate(mode, orig_atom_idx, j_idx, disp,
                              ignore_periodic_self_interactions)) {
        continue;
      }
      if (candidate_count >= soa_x.size()) {
        growScratchBuffers(soa_x.size() * 2);
      }
      soa_x[candidate_count] = sorted_coords.x[slot] + disp.x();
      soa_y[candidate_count] = sorted_coords.y[slot] + disp.y();
      soa_z[candidate_count] = sorted_coords.z[slot] + disp.z();
      candidate_j[candidate_count] = j_idx;
      candidate_count++;
    }
  }

  void collectSmallCellCandidates(size_t orig_atom_idx, const SoACoordinates &sorted_coords,
                                  const size_t *orig_indices, const size_t *bin_offsets,
                                  const SearchGridConfig &grid_config,
                                  bool ignore_periodic_self_interactions) {
    const size_t bin_start = bin_offsets[0];
    const size_t bin_end = bin_offsets[1];
    collectCandidatesFromBin(orig_atom_idx, {0.0, 0.0, 0.0}, bin_start, bin_end,
                             CandidateGatherMode::PrimaryCellOnly, sorted_coords, orig_indices,
                             ignore_periodic_self_interactions);
    for (const auto &disp : grid_config.small_cell_shifts) {
      collectCandidatesFromBin(orig_atom_idx, disp, bin_start, bin_end,
                               CandidateGatherMode::PeriodicImageShift, sorted_coords, orig_indices,
                               ignore_periodic_self_interactions);
    }
  }

  void collectPartitionedInterior(size_t orig_atom_idx, int current_bin,
                                  const SoACoordinates &sorted_coords, const size_t *orig_indices,
                                  const size_t *bin_offsets, const SearchGridConfig &grid_config,
                                  bool ignore_periodic_self_interactions) {
    const correlation::math::Vector3<real_t> zero_disp{0.0, 0.0, 0.0};
    for (const int delta : grid_config.interior_bin_deltas) {
      const int neighbor_bin = current_bin + delta;
      const size_t bin_start = bin_offsets[neighbor_bin];
      const size_t bin_end = bin_offsets[neighbor_bin + 1];
      collectCandidatesFromBin(orig_atom_idx, zero_disp, bin_start, bin_end,
                               CandidateGatherMode::PARTITIONED, sorted_coords, orig_indices,
                               ignore_periodic_self_interactions);
    }
  }

  void collectBoundaryDisplacementSlice(size_t orig_atom_idx, int c_x, int dx, int c_y, int dy,
                                        int c_z, const SoACoordinates &sorted_coords,
                                        const size_t *orig_indices, const size_t *bin_offsets,
                                        const correlation::math::Matrix3<real_t> &lattice,
                                        const SearchGridConfig &grid_config,
                                        bool ignore_periodic_self_interactions) {
    int shift_x = 0;
    int wrap_x = 0;
    computePeriodicShiftAndBin(c_x, dx, grid_config.k_x, shift_x, wrap_x);

    int shift_y = 0;
    int wrap_y = 0;
    computePeriodicShiftAndBin(c_y, dy, grid_config.k_y, shift_y, wrap_y);

    for (int dz = -grid_config.max_dz; dz <= grid_config.max_dz; ++dz) {
      int shift_z = 0;
      int wrap_z = 0;
      computePeriodicShiftAndBin(c_z, dz, grid_config.k_z, shift_z, wrap_z);

      const int n_bin_idx =
          wrap_x * (grid_config.k_y * grid_config.k_z) + wrap_y * grid_config.k_z + wrap_z;

      const real_t disp_x = std::fma(static_cast<real_t>(shift_z), lattice[2].x(),
                                     std::fma(static_cast<real_t>(shift_y), lattice[1].x(),
                                              static_cast<real_t>(shift_x) * lattice[0].x()));
      const real_t disp_y = std::fma(static_cast<real_t>(shift_z), lattice[2].y(),
                                     std::fma(static_cast<real_t>(shift_y), lattice[1].y(),
                                              static_cast<real_t>(shift_x) * lattice[0].y()));
      const real_t disp_z = std::fma(static_cast<real_t>(shift_z), lattice[2].z(),
                                     std::fma(static_cast<real_t>(shift_y), lattice[1].z(),
                                              static_cast<real_t>(shift_x) * lattice[0].z()));
      const correlation::math::Vector3<real_t> disp{disp_x, disp_y, disp_z};

      const size_t bin_start = bin_offsets[n_bin_idx];
      const size_t bin_end = bin_offsets[n_bin_idx + 1];

      collectCandidatesFromBin(orig_atom_idx, disp, bin_start, bin_end,
                               CandidateGatherMode::PARTITIONED, sorted_coords, orig_indices,
                               ignore_periodic_self_interactions);
    }
  }

  void collectPartitionedBoundary(size_t orig_atom_idx, int current_bin,
                                  const SoACoordinates &sorted_coords, const size_t *orig_indices,
                                  const size_t *bin_offsets,
                                  const correlation::math::Matrix3<real_t> &lattice,
                                  const SearchGridConfig &grid_config,
                                  bool ignore_periodic_self_interactions) {
    const int c_x = current_bin / (grid_config.k_y * grid_config.k_z);
    const int c_y = (current_bin / grid_config.k_z) % grid_config.k_y;
    const int c_z = current_bin % grid_config.k_z;

    for (int dx = -grid_config.max_dx; dx <= grid_config.max_dx; ++dx) {
      for (int dy = -grid_config.max_dy; dy <= grid_config.max_dy; ++dy) {
        collectBoundaryDisplacementSlice(orig_atom_idx, c_x, dx, c_y, dy, c_z, sorted_coords,
                                         orig_indices, bin_offsets, lattice, grid_config,
                                         ignore_periodic_self_interactions);
      }
    }
  }

  void collectCandidates(size_t orig_atom_idx, int current_bin, const SoACoordinates &sorted_coords,
                         const size_t *orig_indices, const size_t *bin_offsets,
                         const correlation::math::Matrix3<real_t> &lattice,
                         const SearchGridConfig &grid_config,
                         bool ignore_periodic_self_interactions) {
    candidate_count = 0;
    if (grid_config.is_small_cell) {
      collectSmallCellCandidates(orig_atom_idx, sorted_coords, orig_indices, bin_offsets,
                                 grid_config, ignore_periodic_self_interactions);
    } else {
      const int c_x = current_bin / (grid_config.k_y * grid_config.k_z);
      const int c_y = (current_bin / grid_config.k_z) % grid_config.k_y;
      const int c_z = current_bin % grid_config.k_z;

      const bool is_interior = (c_x > 0 && c_x < grid_config.k_x - 1) &&
                               (c_y > 0 && c_y < grid_config.k_y - 1) &&
                               (c_z > 0 && c_z < grid_config.k_z - 1);

      if (is_interior) {
        collectPartitionedInterior(orig_atom_idx, current_bin, sorted_coords, orig_indices,
                                   bin_offsets, grid_config, ignore_periodic_self_interactions);
      } else {
        collectPartitionedBoundary(orig_atom_idx, current_bin, sorted_coords, orig_indices,
                                   bin_offsets, lattice, grid_config,
                                   ignore_periodic_self_interactions);
      }
    }
  }
};

struct ShiftBounds {
  int max_dx{0};
  int max_dy{0};
  int max_dz{0};
  real_t width_x{0.0};
  real_t width_y{0.0};
  real_t width_z{0.0};
  real_t cutoff{0.0};
};

inline std::vector<correlation::math::Vector3<real_t>>
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

inline SearchGridConfig buildSearchGridConfig(const correlation::core::Cell &cell,
                                              real_t cutoff_sq) {
  constexpr auto K_ZERO = static_cast<real_t>(0.0);
  const real_t volume = std::abs(cell.volume());

  if (cutoff_sq <= K_ZERO || volume <= K_ZERO) {
    return SearchGridConfig{};
  }

  const auto &lattice = cell.latticeVectors();
  const math::Vector3<real_t> lattice_a = lattice[0];
  const math::Vector3<real_t> lattice_b = lattice[1];
  const math::Vector3<real_t> lattice_c = lattice[2];

  const math::Vector3<real_t> cross_bc = math::cross(lattice_b, lattice_c);
  const math::Vector3<real_t> cross_ca = math::cross(lattice_c, lattice_a);
  const math::Vector3<real_t> cross_ab = math::cross(lattice_a, lattice_b);

  const real_t width_x = volume / std::max(math::norm(cross_bc), static_cast<real_t>(1e-9));
  const real_t width_y = volume / std::max(math::norm(cross_ca), static_cast<real_t>(1e-9));
  const real_t width_z = volume / std::max(math::norm(cross_ab), static_cast<real_t>(1e-9));

  const real_t cutoff = std::sqrt(cutoff_sq);
  const bool is_small_cell = (width_x < static_cast<real_t>(2.0) * cutoff) ||
                             (width_y < static_cast<real_t>(2.0) * cutoff) ||
                             (width_z < static_cast<real_t>(2.0) * cutoff);

  if (is_small_cell) {
    const int max_dx = std::max(1, static_cast<int>(std::ceil(cutoff / width_x)));
    const int max_dy = std::max(1, static_cast<int>(std::ceil(cutoff / width_y)));
    const int max_dz = std::max(1, static_cast<int>(std::ceil(cutoff / width_z)));
    const ShiftBounds bounds{
        .max_dx = max_dx,
        .max_dy = max_dy,
        .max_dz = max_dz,
        .width_x = width_x,
        .width_y = width_y,
        .width_z = width_z,
        .cutoff = cutoff,
    };
    return SearchGridConfig{
        .k_x = 1,
        .k_y = 1,
        .k_z = 1,
        .max_dx = max_dx,
        .max_dy = max_dy,
        .max_dz = max_dz,
        .is_small_cell = true,
        .small_cell_shifts = computeSmallCellShifts(lattice, bounds),
        .interior_bin_deltas = {},
    };
  }

  const int k_x = std::max(3, static_cast<int>(std::floor(width_x / cutoff)));
  const int k_y = std::max(3, static_cast<int>(std::floor(width_y / cutoff)));
  const int k_z = std::max(3, static_cast<int>(std::floor(width_z / cutoff)));

  const real_t d_x = width_x / static_cast<real_t>(k_x);
  const real_t d_y = width_y / static_cast<real_t>(k_y);
  const real_t d_z = width_z / static_cast<real_t>(k_z);

  const int max_dx = std::max(1, static_cast<int>(std::ceil(cutoff / d_x)));
  const int max_dy = std::max(1, static_cast<int>(std::ceil(cutoff / d_y)));
  const int max_dz = std::max(1, static_cast<int>(std::ceil(cutoff / d_z)));

  std::array<int, 27> deltas{};
  size_t idx = 0;
  for (int dx = -1; dx <= 1; ++dx) {
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dz = -1; dz <= 1; ++dz) {
        deltas[idx++] = dx * (k_y * k_z) + dy * k_z + dz;
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
      .is_small_cell = false,
      .small_cell_shifts = {},
      .interior_bin_deltas = deltas,
  };
}

inline BinSortedPositions buildBinSortedPositions(const correlation::core::Cell &cell,
                                                  const SearchGridConfig &grid_config) {
  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const auto &inv_lattice = cell.inverseLatticeVectors();

  std::vector<real_t> wrapped_x(atom_count);
  std::vector<real_t> wrapped_y(atom_count);
  std::vector<real_t> wrapped_z(atom_count);
  std::vector<int> atom_bin(atom_count);

  const size_t num_bins = static_cast<size_t>(grid_config.k_x * grid_config.k_y * grid_config.k_z);
  std::vector<size_t> bin_counts(num_bins, 0);

  const real_t k_x_f = static_cast<real_t>(grid_config.k_x);
  const real_t k_y_f = static_cast<real_t>(grid_config.k_y);
  const real_t k_z_f = static_cast<real_t>(grid_config.k_z);

  for (size_t i = 0; i < atom_count; ++i) {
    const auto &pos = atoms[i].position();
    real_t u =
        inv_lattice[0].x() * pos.x() + inv_lattice[1].x() * pos.y() + inv_lattice[2].x() * pos.z();
    real_t v =
        inv_lattice[0].y() * pos.x() + inv_lattice[1].y() * pos.y() + inv_lattice[2].y() * pos.z();
    real_t w =
        inv_lattice[0].z() * pos.x() + inv_lattice[1].z() * pos.y() + inv_lattice[2].z() * pos.z();

    u -= std::floor(u);
    v -= std::floor(v);
    w -= std::floor(w);

    const auto &lattice = cell.latticeVectors();
    wrapped_x[i] = u * lattice[0].x() + v * lattice[1].x() + w * lattice[2].x();
    wrapped_y[i] = u * lattice[0].y() + v * lattice[1].y() + w * lattice[2].y();
    wrapped_z[i] = u * lattice[0].z() + v * lattice[1].z() + w * lattice[2].z();

    int b_x = static_cast<int>(u * k_x_f);
    int b_y = static_cast<int>(v * k_y_f);
    int b_z = static_cast<int>(w * k_z_f);

    b_x = std::clamp(b_x, 0, grid_config.k_x - 1);
    b_y = std::clamp(b_y, 0, grid_config.k_y - 1);
    b_z = std::clamp(b_z, 0, grid_config.k_z - 1);

    const int bin_idx = b_x * (grid_config.k_y * grid_config.k_z) + b_y * grid_config.k_z + b_z;
    atom_bin[i] = bin_idx;
    bin_counts[bin_idx]++;
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

template <typename Visitor>
void forEachPairWithin(const correlation::core::Cell &cell, real_t r_max,
                       bool ignore_periodic_self_interactions, Visitor &&visitor) {
  if (r_max <= static_cast<real_t>(0.0)) {
    return;
  }
  const real_t cutoff_sq = r_max * r_max;
  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  if (atom_count == 0) {
    return;
  }

  const auto &lattice = cell.latticeVectors();
  const SearchGridConfig grid_config = buildSearchGridConfig(cell, cutoff_sq);
  const BinSortedPositions bin_sorted = buildBinSortedPositions(cell, grid_config);

  const SoACoordinates sorted_coords{
      .x = bin_sorted.x_sorted.data(),
      .y = bin_sorted.y_sorted.data(),
      .z = bin_sorted.z_sorted.data(),
  };

  tbb::enumerable_thread_specific<CandidateGatherer> ets([] {
    CandidateGatherer gatherer;
    gatherer.growScratchBuffers(1024);
    return gatherer;
  });

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, atom_count), [&](const tbb::blocked_range<size_t> &range) {
        auto &gatherer = ets.local();
        for (size_t slot = range.begin(); slot != range.end(); ++slot) {
          const size_t orig_idx = bin_sorted.orig_indices[slot];
          const int current_bin = bin_sorted.sorted_atom_bin[slot];
          const correlation::math::Vector3<real_t> atom_pos{
              bin_sorted.x_sorted[slot], bin_sorted.y_sorted[slot], bin_sorted.z_sorted[slot]};

          gatherer.collectCandidates(orig_idx, current_bin, sorted_coords,
                                     bin_sorted.orig_indices.data(), bin_sorted.bin_offsets.data(),
                                     lattice, grid_config, ignore_periodic_self_interactions);

          const size_t c_count = gatherer.candidate_count;
          if (c_count == 0) {
            continue;
          }

          const correlation::math::PositionBlock block{
              .x = gatherer.soa_x.data(),
              .y = gatherer.soa_y.data(),
              .z = gatherer.soa_z.data(),
              .count = c_count,
          };
          correlation::math::compute_dsq_block(atom_pos.x(), atom_pos.y(), atom_pos.z(), block,
                                               gatherer.dsq_scratch.data());

          for (size_t k = 0; k < c_count; ++k) {
            const real_t d_sq = gatherer.dsq_scratch[k];
            if (d_sq >= cutoff_sq) {
              continue;
            }
            const size_t j_idx = gatherer.candidate_j[k];
            const real_t r_x = gatherer.soa_x[k] - atom_pos.x();
            const real_t r_y = gatherer.soa_y[k] - atom_pos.y();
            const real_t r_z = gatherer.soa_z[k] - atom_pos.z();
            const correlation::math::Vector3<real_t> r_ij{r_x, r_y, r_z};

            visitor(orig_idx, j_idx, d_sq, r_ij);
          }
        }
      });
}

} // namespace correlation::calculators::detail
