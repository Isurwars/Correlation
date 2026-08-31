/**
 * @file GPUDistanceCalculator.cu
 * @brief GPU implementation of pairwise distance calculation supporting float and double precision.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/GPUDistanceCalculator.hpp"
#include "core/DeviceBuffer.hpp"
#include "core/GPUErrorCheck.hpp"
#include "core/GPUPortability.hpp"
#include "math/LinearAlgebra.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace correlation::calculators::gpu {

namespace {

template <typename T> struct GPUBond {
  int from;
  int to;
  T distance;
  T r_x;
  T r_y;
  T r_z;
};

template <typename T> struct GPULattice {
  T v0_x, v0_y, v0_z;
  T v1_x, v1_y, v1_z;
  T v2_x, v2_y, v2_z;
};

struct GPUSearchGrid {
  int K_x, K_y, K_z;
  int max_dx, max_dy, max_dz;
};

template <typename T> struct GPUPosition {
  T x;
  T y;
  T z;
};

template <typename T> struct GPUAtomData {
  const T *CORRELATION_RESTRICT wrapped_x;
  const T *CORRELATION_RESTRICT wrapped_y;
  const T *CORRELATION_RESTRICT wrapped_z;
  const int *CORRELATION_RESTRICT element_ids;
  const int *CORRELATION_RESTRICT atom_bin;
};

struct GPUBinData {
  const unsigned long long *CORRELATION_RESTRICT offsets;
  const unsigned long long *CORRELATION_RESTRICT indices;
};

CORRELATION_DEVICE CORRELATION_FORCEINLINE void wrap_coordinate(int bin, int k_val, int &wrap, int &shift) {
  shift = (bin >= 0) ? (bin / k_val) : ((bin - k_val + 1) / k_val);
  wrap = bin - shift * k_val;
}

CORRELATION_DEVICE CORRELATION_FORCEINLINE bool should_skip_atom_pair(int i_val, int j_val, bool zero_disp,
                                                                      bool ignore_periodic_self) {
  if (j_val < i_val) {
    return true;
  }
  if (i_val == j_val && (zero_disp || ignore_periodic_self)) {
    return true;
  }
  return false;
}

struct AtomPair {
  int i_val;
  int j_val;
  int type_A;
  int type_B;
};

template <typename T> struct PairPositions {
  GPUPosition<T> atom_pos;
  GPUPosition<T> shifted;
};

template <typename T>
CORRELATION_DEVICE CORRELATION_FORCEINLINE void
accumulate_histograms(AtomPair pair, T dist, int num_elements, T hist_r_max, T hist_r_bin_width, int hist_num_bins,
                      unsigned long long *CORRELATION_RESTRICT d_histograms) {
  if (d_histograms == nullptr || hist_num_bins <= 0 || hist_r_bin_width <= static_cast<T>(0.0) || dist >= hist_r_max) {
    return;
  }

  int const bin_idx = static_cast<int>(dist / hist_r_bin_width);
  if (bin_idx >= 0 && bin_idx < hist_num_bins) {
    atomicAdd(&d_histograms[pair.type_A * (num_elements * hist_num_bins) + pair.type_B * hist_num_bins + bin_idx],
              1ULL);
    if (pair.i_val != pair.j_val && pair.type_A != pair.type_B) {
      atomicAdd(&d_histograms[pair.type_B * (num_elements * hist_num_bins) + pair.type_A * hist_num_bins + bin_idx],
                1ULL);
    }
  }
}

template <bool WriteBonds, typename T>
CORRELATION_DEVICE CORRELATION_FORCEINLINE void
record_bond_interaction(AtomPair pair, T d_sq, PairPositions<T> positions, GPUAtomData<T> atoms, int num_elements,
                        const T *CORRELATION_RESTRICT bond_cutoffs_sq,
                        unsigned long long *CORRELATION_RESTRICT bond_counter, GPUBond<T> *CORRELATION_RESTRICT bonds) {
  if (bond_cutoffs_sq == nullptr || pair.type_A >= num_elements || pair.type_B >= num_elements) {
    return;
  }

  T const max_bond_dist_sq = bond_cutoffs_sq[pair.type_A * num_elements + pair.type_B];
  if (max_bond_dist_sq <= static_cast<T>(0.0) || d_sq > max_bond_dist_sq) {
    return;
  }

  if constexpr (WriteBonds) {
    T const precise_disp_x = positions.shifted.x - atoms.wrapped_x[pair.j_val];
    T const precise_disp_y = positions.shifted.y - atoms.wrapped_y[pair.j_val];
    T const precise_disp_z = positions.shifted.z - atoms.wrapped_z[pair.j_val];

    T const r_x = (atoms.wrapped_x[pair.j_val] - positions.atom_pos.x) + precise_disp_x;
    T const r_y = (atoms.wrapped_y[pair.j_val] - positions.atom_pos.y) + precise_disp_y;
    T const r_z = (atoms.wrapped_z[pair.j_val] - positions.atom_pos.z) + precise_disp_z;

    T const precise_dist = sqrt(r_x * r_x + r_y * r_y + r_z * r_z);

    unsigned long long const idx = atomicAdd(bond_counter, 1ULL);
    bonds[idx] = GPUBond<T>{
        .from = pair.i_val,
        .to = pair.j_val,
        .distance = precise_dist,
        .r_x = r_x,
        .r_y = r_y,
        .r_z = r_z,
    };
  } else if (bond_counter != nullptr) {
    atomicAdd(bond_counter, 1ULL);
  }
}

template <bool WriteBonds, typename T>
CORRELATION_DEVICE CORRELATION_FORCEINLINE void
process_bin(int i_val, GPUPosition<T> atom_pos, int type_A, GPUPosition<T> disp, int n_bin_idx, GPUAtomData<T> atoms,
            bool zero_disp, GPUBinData bins, T cutoff_sq, const T *CORRELATION_RESTRICT bond_cutoffs_sq,
            int num_elements, bool ignore_periodic_self_interactions,
            unsigned long long *CORRELATION_RESTRICT bond_counter, GPUBond<T> *CORRELATION_RESTRICT bonds, T hist_r_max,
            T hist_r_bin_width, int hist_num_bins, unsigned long long *CORRELATION_RESTRICT d_histograms) {

  unsigned long long const start = bins.offsets[n_bin_idx];
  unsigned long long const end = bins.offsets[n_bin_idx + 1];

  for (unsigned long long offset = start; offset < end; ++offset) {
    int const j_val = static_cast<int>(bins.indices[offset]);
    if (should_skip_atom_pair(i_val, j_val, zero_disp, ignore_periodic_self_interactions)) {
      continue;
    }

    T const shifted_x = atoms.wrapped_x[j_val] + disp.x;
    T const shifted_y = atoms.wrapped_y[j_val] + disp.y;
    T const shifted_z = atoms.wrapped_z[j_val] + disp.z;

    T const dx_val = shifted_x - atom_pos.x;
    T const dy_val = shifted_y - atom_pos.y;
    T const dz_val = shifted_z - atom_pos.z;

    T const d_sq = dx_val * dx_val + dy_val * dy_val + dz_val * dz_val;
    if (d_sq >= cutoff_sq) {
      continue;
    }

    T const dist = sqrt(d_sq);
    int const type_B = atoms.element_ids[j_val];
    AtomPair const pair{.i_val = i_val, .j_val = j_val, .type_A = type_A, .type_B = type_B};
    GPUPosition<T> const shifted{.x = shifted_x, .y = shifted_y, .z = shifted_z};

    accumulate_histograms(pair, dist, num_elements, hist_r_max, hist_r_bin_width, hist_num_bins, d_histograms);

    PairPositions<T> const positions{.atom_pos = atom_pos, .shifted = shifted};
    record_bond_interaction<WriteBonds>(pair, d_sq, positions, atoms, num_elements, bond_cutoffs_sq, bond_counter,
                                        bonds);
  }
}

template <bool WriteBonds, typename T>
CORRELATION_GLOBAL void distance_kernel(GPUAtomData<T> atoms, GPUBinData bins, GPULattice<T> lattice,
                                        GPUSearchGrid grid, T cutoff_sq, const T *CORRELATION_RESTRICT bond_cutoffs_sq,
                                        int num_elements, bool ignore_periodic_self_interactions, int num_atoms,
                                        unsigned long long *CORRELATION_RESTRICT bond_counter,
                                        GPUBond<T> *CORRELATION_RESTRICT bonds, T hist_r_max, T hist_r_bin_width,
                                        int hist_num_bins, unsigned long long *CORRELATION_RESTRICT d_histograms) {

  int const i_val = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
  if (i_val >= num_atoms) {
    return;
  }

  T const a_x = atoms.wrapped_x[i_val];
  T const a_y = atoms.wrapped_y[i_val];
  T const a_z = atoms.wrapped_z[i_val];
  int const type_A = atoms.element_ids[i_val];

  int const c_bin = atoms.atom_bin[i_val];
  int const c_x = c_bin / (grid.K_y * grid.K_z);
  int const c_y = (c_bin / grid.K_z) % grid.K_y;
  int const c_z = c_bin % grid.K_z;

  for (int dx = -grid.max_dx; dx <= grid.max_dx; ++dx) {
    int wrap_x = 0;
    int shift_x = 0;
    wrap_coordinate(c_x + dx, grid.K_x, wrap_x, shift_x);

    for (int dy = -grid.max_dy; dy <= grid.max_dy; ++dy) {
      int wrap_y = 0;
      int shift_y = 0;
      wrap_coordinate(c_y + dy, grid.K_y, wrap_y, shift_y);

      for (int dz = -grid.max_dz; dz <= grid.max_dz; ++dz) {
        int wrap_z = 0;
        int shift_z = 0;
        wrap_coordinate(c_z + dz, grid.K_z, wrap_z, shift_z);

        int const n_bin_idx = wrap_x * (grid.K_y * grid.K_z) + wrap_y * grid.K_z + wrap_z;

        T const disp_x = static_cast<T>(shift_z) * lattice.v2_x + static_cast<T>(shift_y) * lattice.v1_x +
                         static_cast<T>(shift_x) * lattice.v0_x;
        T const disp_y = static_cast<T>(shift_z) * lattice.v2_y + static_cast<T>(shift_y) * lattice.v1_y +
                         static_cast<T>(shift_x) * lattice.v0_y;
        T const disp_z = static_cast<T>(shift_z) * lattice.v2_z + static_cast<T>(shift_y) * lattice.v1_z +
                         static_cast<T>(shift_x) * lattice.v0_z;

        bool const zero_disp = (shift_x == 0 && shift_y == 0 && shift_z == 0);

        process_bin<WriteBonds, T>(i_val,
                                   GPUPosition<T>{
                                       .x = a_x,
                                       .y = a_y,
                                       .z = a_z,
                                   },
                                   type_A,
                                   GPUPosition<T>{
                                       .x = disp_x,
                                       .y = disp_y,
                                       .z = disp_z,
                                   },
                                   n_bin_idx, atoms, zero_disp, bins, cutoff_sq, bond_cutoffs_sq, num_elements,
                                   ignore_periodic_self_interactions, bond_counter, bonds, hist_r_max, hist_r_bin_width,
                                   hist_num_bins, d_histograms);
      }
    }
  }
}

template <typename T>
std::vector<T> flatten_bond_cutoffs(const std::vector<std::vector<T>> &bond_cutoffs_sq, size_t num_elements) {
  std::vector<T> flat_cutoffs(num_elements * num_elements, static_cast<T>(0.0));
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
      if (i < bond_cutoffs_sq.size() && j < bond_cutoffs_sq[i].size()) {
        flat_cutoffs[i * num_elements + j] = bond_cutoffs_sq[i][j];
      }
    }
  }
  return flat_cutoffs;
}

template <typename T>
void unpack_gpu_bonds(const std::vector<GPUBond<T>> &host_bonds, correlation::core::NeighborGraph &out_graph) {
  for (const auto &bond : host_bonds) {
    out_graph.addDirectedEdge(
        static_cast<core::AtomID>(bond.from), static_cast<core::AtomID>(bond.to), static_cast<real_t>(bond.distance),
        correlation::math::Vector3<real_t>(static_cast<real_t>(bond.r_x), static_cast<real_t>(bond.r_y),
                                           static_cast<real_t>(bond.r_z)));
    if (bond.from != bond.to) {
      out_graph.addDirectedEdge(
          static_cast<core::AtomID>(bond.to), static_cast<core::AtomID>(bond.from), static_cast<real_t>(bond.distance),
          correlation::math::Vector3<real_t>(-static_cast<real_t>(bond.r_x), -static_cast<real_t>(bond.r_y),
                                             -static_cast<real_t>(bond.r_z)));
    }
  }
}

template <typename T> struct SpatialPartitionData {
  GPUSearchGrid grid{};
  std::vector<T> wrapped_x;
  std::vector<T> wrapped_y;
  std::vector<T> wrapped_z;
  std::vector<int> element_ids;
  std::vector<int> atom_bin;
  std::vector<unsigned long long> bin_offsets;
  std::vector<unsigned long long> bin_indices;
  size_t num_bins{0};
  bool ignore_periodic_self_interactions{false};
};

template <typename T>
SpatialPartitionData<T> build_spatial_partition(const correlation::core::Cell &cell, T cutoff_sq,
                                                bool ignore_periodic_self_interactions) {
  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  const auto &lattice = cell.latticeVectors();

  T const vol = static_cast<T>(cell.volume());
  T const width_x = vol / static_cast<T>(correlation::math::norm(correlation::math::cross(lattice[1], lattice[2])));
  T const width_y = vol / static_cast<T>(correlation::math::norm(correlation::math::cross(lattice[0], lattice[2])));
  T const width_z = vol / static_cast<T>(correlation::math::norm(correlation::math::cross(lattice[0], lattice[1])));

  T const cutoff = std::sqrt(cutoff_sq);

  int const k_x = std::max(1, static_cast<int>(std::floor(width_x / cutoff)));
  int const k_y = std::max(1, static_cast<int>(std::floor(width_y / cutoff)));
  int const k_z = std::max(1, static_cast<int>(std::floor(width_z / cutoff)));

  int const max_dx = (k_x == 1) ? static_cast<int>(std::ceil(cutoff / width_x)) : 1;
  int const max_dy = (k_y == 1) ? static_cast<int>(std::ceil(cutoff / width_y)) : 1;
  int const max_dz = (k_z == 1) ? static_cast<int>(std::ceil(cutoff / width_z)) : 1;

  bool ignore_self = ignore_periodic_self_interactions;
  if (max_dx + max_dy + max_dz > 8) {
    ignore_self = false;
  }

  size_t const num_bins = static_cast<size_t>(k_x) * static_cast<size_t>(k_y) * static_cast<size_t>(k_z);
  std::vector<size_t> bin_counts(num_bins, 0);
  std::vector<int> atom_bin(atom_count);

  std::vector<T> wrapped_x(atom_count);
  std::vector<T> wrapped_y(atom_count);
  std::vector<T> wrapped_z(atom_count);
  std::vector<int> element_ids(atom_count);

  for (size_t i = 0; i < atom_count; ++i) {
    correlation::math::Vector3<real_t> const frac = cell.inverseLatticeVectors() * atoms[i].position();
    T const f_x = static_cast<T>(frac.x() - std::floor(frac.x()));
    T const f_y = static_cast<T>(frac.y() - std::floor(frac.y()));
    T const f_z = static_cast<T>(frac.z() - std::floor(frac.z()));

    wrapped_x[i] = static_cast<T>(std::fma(f_z, lattice[2].x(), std::fma(f_y, lattice[1].x(), f_x * lattice[0].x())));
    wrapped_y[i] = static_cast<T>(std::fma(f_z, lattice[2].y(), std::fma(f_y, lattice[1].y(), f_x * lattice[0].y())));
    wrapped_z[i] = static_cast<T>(std::fma(f_z, lattice[2].z(), std::fma(f_y, lattice[1].z(), f_x * lattice[0].z())));

    int const c_x = std::clamp(static_cast<int>(std::floor(f_x * k_x)), 0, k_x - 1);
    int const c_y = std::clamp(static_cast<int>(std::floor(f_y * k_y)), 0, k_y - 1);
    int const c_z = std::clamp(static_cast<int>(std::floor(f_z * k_z)), 0, k_z - 1);

    int const bin_idx = c_x * (k_y * k_z) + c_y * k_z + c_z;
    atom_bin[i] = bin_idx;
    bin_counts[bin_idx]++;
    element_ids[i] = atoms[i].element_id();
  }

  std::vector<unsigned long long> bin_offsets(num_bins + 1, 0);
  for (size_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
    bin_offsets[bin_idx + 1] = bin_offsets[bin_idx] + bin_counts[bin_idx];
  }

  std::vector<unsigned long long> bin_indices(atom_count);
  std::vector<unsigned long long> insertion_cursors = bin_offsets;
  for (size_t i = 0; i < atom_count; ++i) {
    int const bin_idx = atom_bin[i];
    bin_indices[insertion_cursors[bin_idx]++] = i;
  }

  return SpatialPartitionData<T>{
      .grid = {.K_x = k_x, .K_y = k_y, .K_z = k_z, .max_dx = max_dx, .max_dy = max_dy, .max_dz = max_dz},
      .wrapped_x = std::move(wrapped_x),
      .wrapped_y = std::move(wrapped_y),
      .wrapped_z = std::move(wrapped_z),
      .element_ids = std::move(element_ids),
      .atom_bin = std::move(atom_bin),
      .bin_offsets = std::move(bin_offsets),
      .bin_indices = std::move(bin_indices),
      .num_bins = num_bins,
      .ignore_periodic_self_interactions = ignore_self,
  };
}

template <typename T> bool has_active_bonds(const std::vector<std::vector<T>> &bond_cutoffs_sq) {
  for (const auto &row : bond_cutoffs_sq) {
    for (T val : row) {
      if (val > static_cast<T>(0.0)) {
        return true;
      }
    }
  }
  return false;
}

void copy_and_unpack_histograms(const correlation::core::gpu::DeviceBuffer<unsigned long long> &d_histograms,
                                size_t num_elements, const DistanceCalculationConfig &hist_config,
                                RawHistogramTensor &out_histograms) {
  size_t const hist_tensor_size = num_elements * num_elements * hist_config.num_bins;
  std::vector<unsigned long long> host_histograms(hist_tensor_size);
  d_histograms.copyToHost(host_histograms.data(), hist_tensor_size);

  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
      for (size_t bin_idx = 0; bin_idx < hist_config.num_bins; ++bin_idx) {
        size_t const flat_idx = i * (num_elements * hist_config.num_bins) + j * hist_config.num_bins + bin_idx;
        out_histograms[i][j][bin_idx] = static_cast<real_t>(host_histograms[flat_idx]);
      }
    }
  }
}

template <typename T> struct KernelLaunchConfig {
  int grid_size{0};
  int block_size{256};
};

template <typename T> struct TwoPassExecutionContext {
  GPUAtomData<T> gpu_atoms{};
  GPUBinData gpu_bins{};
  GPULattice<T> gpu_lattice{};
  GPUSearchGrid gpu_grid{};
  T cutoff_sq{};
  size_t num_elements{0};
  size_t atom_count{0};
  bool ignore_periodic_self_interactions{false};
  KernelLaunchConfig<T> launch_config{};
  DistanceCalculationConfig hist_config{};
  unsigned long long *d_histograms_ptr{nullptr};
};

template <typename T>
void execute_two_pass_bonds(const TwoPassExecutionContext<T> &ctx, const std::vector<std::vector<T>> &bond_cutoffs_sq,
                            correlation::core::NeighborGraph &out_graph) {
  std::vector<T> const h_bond_cutoffs_sq = flatten_bond_cutoffs(bond_cutoffs_sq, ctx.num_elements);
  correlation::core::gpu::DeviceBuffer<T> d_bond_cutoffs_sq(ctx.num_elements * ctx.num_elements);
  d_bond_cutoffs_sq.copyFromHost(h_bond_cutoffs_sq.data(), ctx.num_elements * ctx.num_elements);

  correlation::core::gpu::DeviceBuffer<unsigned long long> d_bond_counter(1);
  unsigned long long const zero_val = 0;
  d_bond_counter.setScalar(zero_val);

  // Pass 1: Count bonds (pass d_histograms = nullptr to prevent double-counting)
  hipLaunchKernelGGL((distance_kernel<false, T>), ctx.launch_config.grid_size, ctx.launch_config.block_size, 0, 0,
                     ctx.gpu_atoms, ctx.gpu_bins, ctx.gpu_lattice, ctx.gpu_grid, ctx.cutoff_sq, d_bond_cutoffs_sq.get(),
                     static_cast<int>(ctx.num_elements), ctx.ignore_periodic_self_interactions,
                     static_cast<int>(ctx.atom_count), d_bond_counter.get(), nullptr, static_cast<T>(0),
                     static_cast<T>(0), 0, nullptr);
  correlation::core::gpu::hipCheck(hipDeviceSynchronize());

  unsigned long long h_bond_count = 0;
  d_bond_counter.copyToHost(&h_bond_count, 1);

  correlation::core::gpu::DeviceBuffer<GPUBond<T>> const d_bonds(h_bond_count);
  d_bond_counter.setScalar(zero_val);

  // Pass 2: Write bonds and accumulate device histograms
  hipLaunchKernelGGL((distance_kernel<true, T>), ctx.launch_config.grid_size, ctx.launch_config.block_size, 0, 0,
                     ctx.gpu_atoms, ctx.gpu_bins, ctx.gpu_lattice, ctx.gpu_grid, ctx.cutoff_sq, d_bond_cutoffs_sq.get(),
                     static_cast<int>(ctx.num_elements), ctx.ignore_periodic_self_interactions,
                     static_cast<int>(ctx.atom_count), d_bond_counter.get(), d_bonds.get(),
                     static_cast<T>(ctx.hist_config.r_max), static_cast<T>(ctx.hist_config.r_bin_width),
                     static_cast<int>(ctx.hist_config.num_bins), ctx.d_histograms_ptr);
  correlation::core::gpu::hipCheck(hipDeviceSynchronize());

  std::vector<GPUBond<T>> host_bonds(h_bond_count);
  if (h_bond_count > 0) {
    d_bonds.copyToHost(host_bonds.data(), h_bond_count);
  }
  unpack_gpu_bonds(host_bonds, out_graph);
}

} // namespace

bool has_gpu_device() {
  int device_count = 0;
  hipError_t const err = hipGetDeviceCount(&device_count);
  return (err == hipSuccess && device_count > 0);
}

template <typename T>
void compute_distances_gpu(const correlation::core::Cell &cell, T cutoff_sq,
                           const std::vector<std::vector<T>> &bond_cutoffs_sq, bool ignore_periodic_self_interactions,
                           RawHistogramTensor *out_histograms, DistanceCalculationConfig hist_config,
                           correlation::core::NeighborGraph &out_graph) {
  const size_t atom_count = cell.atoms().size();
  if (atom_count == 0) {
    return;
  }
  const size_t num_elements = cell.elements().size();
  const auto &lattice = cell.latticeVectors();

  auto const partition = build_spatial_partition(cell, cutoff_sq, ignore_periodic_self_interactions);

  bool const has_bonds = has_active_bonds(bond_cutoffs_sq);
  bool const has_histograms = (out_histograms != nullptr && hist_config.num_bins > 0 && hist_config.r_bin_width > 0.0);
  size_t const hist_tensor_size = has_histograms ? (num_elements * num_elements * hist_config.num_bins) : 0;

  // --- RAII device allocations ---
  using correlation::core::gpu::DeviceBuffer;

  DeviceBuffer<T> d_wrapped_x(atom_count);
  DeviceBuffer<T> d_wrapped_y(atom_count);
  DeviceBuffer<T> d_wrapped_z(atom_count);
  DeviceBuffer<int> d_element_ids(atom_count);
  DeviceBuffer<int> d_atom_bin(atom_count);
  DeviceBuffer<unsigned long long> d_bin_offsets(partition.num_bins + 1);
  DeviceBuffer<unsigned long long> d_bin_indices(atom_count);
  DeviceBuffer<unsigned long long> d_histograms(hist_tensor_size > 0 ? hist_tensor_size : 1);

  d_wrapped_x.copyFromHost(partition.wrapped_x.data(), atom_count);
  d_wrapped_y.copyFromHost(partition.wrapped_y.data(), atom_count);
  d_wrapped_z.copyFromHost(partition.wrapped_z.data(), atom_count);
  d_element_ids.copyFromHost(partition.element_ids.data(), atom_count);
  d_atom_bin.copyFromHost(partition.atom_bin.data(), atom_count);
  d_bin_offsets.copyFromHost(partition.bin_offsets.data(), partition.num_bins + 1);
  d_bin_indices.copyFromHost(partition.bin_indices.data(), atom_count);

  if (has_histograms) {
    std::vector<unsigned long long> const zero_hist(hist_tensor_size, 0ULL);
    d_histograms.copyFromHost(zero_hist.data(), hist_tensor_size);
  }

  GPULattice<T> const gpu_lattice{
      static_cast<T>(lattice[0].x()), static_cast<T>(lattice[0].y()), static_cast<T>(lattice[0].z()),
      static_cast<T>(lattice[1].x()), static_cast<T>(lattice[1].y()), static_cast<T>(lattice[1].z()),
      static_cast<T>(lattice[2].x()), static_cast<T>(lattice[2].y()), static_cast<T>(lattice[2].z())};
  GPUAtomData<T> const gpu_atoms{d_wrapped_x.get(), d_wrapped_y.get(), d_wrapped_z.get(), d_element_ids.get(),
                                 d_atom_bin.get()};
  GPUBinData const gpu_bins{.offsets = d_bin_offsets.get(), .indices = d_bin_indices.get()};

  int const block_size = 256;
  int const grid_size = (static_cast<int>(atom_count) + block_size - 1) / block_size;

  if (!has_bonds) {
    // Single-pass direct histogramming kernel (no pair/bond allocation overhead)
    hipLaunchKernelGGL((distance_kernel<false, T>), grid_size, block_size, 0, 0, gpu_atoms, gpu_bins, gpu_lattice,
                       partition.grid, cutoff_sq, nullptr, static_cast<int>(num_elements),
                       partition.ignore_periodic_self_interactions, static_cast<int>(atom_count), nullptr, nullptr,
                       static_cast<T>(hist_config.r_max), static_cast<T>(hist_config.r_bin_width),
                       static_cast<int>(hist_config.num_bins), has_histograms ? d_histograms.get() : nullptr);
    correlation::core::gpu::hipCheck(hipDeviceSynchronize());
  } else {
    TwoPassExecutionContext<T> const ctx{
        .gpu_atoms = gpu_atoms,
        .gpu_bins = gpu_bins,
        .gpu_lattice = gpu_lattice,
        .gpu_grid = partition.grid,
        .cutoff_sq = cutoff_sq,
        .num_elements = num_elements,
        .atom_count = atom_count,
        .ignore_periodic_self_interactions = partition.ignore_periodic_self_interactions,
        .launch_config = {.grid_size = grid_size, .block_size = block_size},
        .hist_config = hist_config,
        .d_histograms_ptr = has_histograms ? d_histograms.get() : nullptr,
    };
    execute_two_pass_bonds(ctx, bond_cutoffs_sq, out_graph);
  }

  // Copy histograms to host
  if (has_histograms) {
    copy_and_unpack_histograms(d_histograms, num_elements, hist_config, *out_histograms);
  }
}

template void compute_distances_gpu<float>(const correlation::core::Cell &cell, float cutoff_sq,
                                           const std::vector<std::vector<float>> &bond_cutoffs_sq,
                                           bool ignore_periodic_self_interactions, RawHistogramTensor *out_histograms,
                                           DistanceCalculationConfig hist_config,
                                           correlation::core::NeighborGraph &out_graph);

template void compute_distances_gpu<double>(const correlation::core::Cell &cell, double cutoff_sq,
                                            const std::vector<std::vector<double>> &bond_cutoffs_sq,
                                            bool ignore_periodic_self_interactions, RawHistogramTensor *out_histograms,
                                            DistanceCalculationConfig hist_config,
                                            correlation::core::NeighborGraph &out_graph);

} // namespace correlation::calculators::gpu
