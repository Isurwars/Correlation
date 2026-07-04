/**
 * @file GPUDistanceCalculator.cu
 * @brief GPU implementation of pairwise distance calculation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/GPUDistanceCalculator.hpp"
#include "core/GPUPortability.hpp"
#include "math/LinearAlgebra.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace correlation::calculators::gpu {

namespace {

struct GPUDistance {
  int from;
  int to;
  double distance;
};

struct GPUBond {
  int from;
  int to;
  double distance;
  double r_x;
  double r_y;
  double r_z;
};

struct GPULattice {
  double v0_x, v0_y, v0_z;
  double v1_x, v1_y, v1_z;
  double v2_x, v2_y, v2_z;
};

struct GPUSearchGrid {
  int K_x, K_y, K_z;
  int max_dx, max_dy, max_dz;
};

struct GPUPosition {
  double x;
  double y;
  double z;
};

struct GPUAtomData {
  const double *__restrict__ wrapped_x;
  const double *__restrict__ wrapped_y;
  const double *__restrict__ wrapped_z;
  const int *__restrict__ element_ids;
  const int *__restrict__ atom_bin;
};

struct GPUBinData {
  const unsigned long long *__restrict__ offsets;
  const unsigned long long *__restrict__ indices;
};

__device__ __forceinline__ void wrap_coordinate(int bin, int k_val, int &wrap, int &shift) {
  shift = (bin >= 0) ? (bin / k_val) : ((bin - k_val + 1) / k_val);
  wrap = bin - shift * k_val;
}

template <bool WriteMode>
__device__ __forceinline__ void
process_bin(int i_val, GPUPosition atom_pos, int type_A, GPUPosition disp, int n_bin_idx, GPUAtomData atoms,
            bool zero_disp, GPUBinData bins, double cutoff_sq, const double *__restrict__ bond_cutoffs_sq,
            int num_elements, unsigned long long *distance_counter, bool ignore_periodic_self_interactions,
            unsigned long long *bond_counter, GPUDistance *__restrict__ distances, GPUBond *__restrict__ bonds) {

  unsigned long long start = bins.offsets[n_bin_idx];
  unsigned long long end = bins.offsets[n_bin_idx + 1];

  for (unsigned long long offset = start; offset < end; ++offset) {
    int j_val = static_cast<int>(bins.indices[offset]);
    if (j_val < i_val) {
      continue;
    }
    if (i_val == j_val && zero_disp) {
      continue;
    }
    if (i_val == j_val && ignore_periodic_self_interactions) {
      continue;
    }

    double shifted_x = atoms.wrapped_x[j_val] + disp.x;
    double shifted_y = atoms.wrapped_y[j_val] + disp.y;
    double shifted_z = atoms.wrapped_z[j_val] + disp.z;

    double dx_val = shifted_x - atom_pos.x;
    double dy_val = shifted_y - atom_pos.y;
    double dz_val = shifted_z - atom_pos.z;

    double d_sq = dx_val * dx_val + dy_val * dy_val + dz_val * dz_val;
    if (d_sq >= cutoff_sq) {
      continue;
    }

    double dist = sqrt(d_sq);
    int type_B = atoms.element_ids[j_val];

    // Increment distance count / write
    if (WriteMode) {
      unsigned long long idx = atomicAdd(distance_counter, 1ULL);
      distances[idx] = GPUDistance{.from = i_val, .to = j_val, .distance = dist};
    } else {
      atomicAdd(distance_counter, 1ULL);
    }

    double max_bond_dist_sq = 0.0;
    if (type_A < num_elements && type_B < num_elements) {
      max_bond_dist_sq = bond_cutoffs_sq[type_A * num_elements + type_B];
    }

    if (d_sq <= max_bond_dist_sq) {
      // Recover displacement vector coordinates precisely to minimize rounding/cancellation errors
      double precise_disp_x = shifted_x - atoms.wrapped_x[j_val];
      double precise_disp_y = shifted_y - atoms.wrapped_y[j_val];
      double precise_disp_z = shifted_z - atoms.wrapped_z[j_val];

      double r_x = (atoms.wrapped_x[j_val] - atom_pos.x) + precise_disp_x;
      double r_y = (atoms.wrapped_y[j_val] - atom_pos.y) + precise_disp_y;
      double r_z = (atoms.wrapped_z[j_val] - atom_pos.z) + precise_disp_z;

      // Recompute consistent distance directly from the exact relative vector
      double precise_dist = sqrt(r_x * r_x + r_y * r_y + r_z * r_z);

      if (WriteMode) {
        unsigned long long idx = atomicAdd(bond_counter, 1ULL);
        bonds[idx] = GPUBond{.from = i_val, .to = j_val, .distance = precise_dist, .r_x = r_x, .r_y = r_y, .r_z = r_z};
      } else {
        atomicAdd(bond_counter, 1ULL);
      }
    }
  }
}

template <bool WriteMode>
__global__ void distance_kernel(GPUAtomData atoms, GPUBinData bins, GPULattice lattice, GPUSearchGrid grid,
                                double cutoff_sq, const double *__restrict__ bond_cutoffs_sq, int num_elements,
                                bool ignore_periodic_self_interactions, int num_atoms,
                                unsigned long long *distance_counter, unsigned long long *bond_counter,
                                GPUDistance *__restrict__ distances, GPUBond *__restrict__ bonds) {

  int i_val = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
  if (i_val >= num_atoms) {
    return;
  }

  double a_x = atoms.wrapped_x[i_val];
  double a_y = atoms.wrapped_y[i_val];
  double a_z = atoms.wrapped_z[i_val];
  int type_A = atoms.element_ids[i_val];

  int c_bin = atoms.atom_bin[i_val];
  int c_x = c_bin / (grid.K_y * grid.K_z);
  int c_y = (c_bin / grid.K_z) % grid.K_y;
  int c_z = c_bin % grid.K_z;

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

        int n_bin_idx = wrap_x * (grid.K_y * grid.K_z) + wrap_y * grid.K_z + wrap_z;

        // Displacement vector coordinates
        double disp_x = shift_z * lattice.v2_x + shift_y * lattice.v1_x + shift_x * lattice.v0_x;
        double disp_y = shift_z * lattice.v2_y + shift_y * lattice.v1_y + shift_x * lattice.v0_y;
        double disp_z = shift_z * lattice.v2_z + shift_y * lattice.v1_z + shift_x * lattice.v0_z;

        bool zero_disp = (shift_x == 0 && shift_y == 0 && shift_z == 0);

        process_bin<WriteMode>(i_val, GPUPosition{.x = a_x, .y = a_y, .z = a_z}, type_A,
                               GPUPosition{.x = disp_x, .y = disp_y, .z = disp_z}, n_bin_idx, atoms, zero_disp, bins,
                               cutoff_sq, bond_cutoffs_sq, num_elements, distance_counter,
                               ignore_periodic_self_interactions, bond_counter, distances, bonds);
      }
    }
  }
}

std::vector<double> flatten_bond_cutoffs(const std::vector<std::vector<double>> &bond_cutoffs_sq, size_t num_elements) {
  std::vector<double> flat_cutoffs(num_elements * num_elements, 0.0);
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
      if (i < bond_cutoffs_sq.size() && j < bond_cutoffs_sq[i].size()) {
        flat_cutoffs[i * num_elements + j] = bond_cutoffs_sq[i][j];
      }
    }
  }
  return flat_cutoffs;
}

void unpack_gpu_results(const std::vector<GPUDistance> &host_distances,
                        const std::vector<GPUBond> &host_bonds,
                        const std::vector<int> &element_ids,
                        DistanceTensor &out_distances,
                        correlation::core::NeighborGraph &out_graph) {
  for (const auto &dist : host_distances) {
    const int type_A = element_ids[dist.from];
    const int type_B = element_ids[dist.to];
    out_distances[type_A][type_B].push_back(dist.distance);
    if (dist.from != dist.to && type_A != type_B) {
      out_distances[type_B][type_A].push_back(dist.distance);
    }
  }

  for (const auto &bond : host_bonds) {
    out_graph.addDirectedEdge(bond.from, bond.to, bond.distance, correlation::math::Vector3<double>(bond.r_x, bond.r_y, bond.r_z));
    if (bond.from != bond.to) {
      out_graph.addDirectedEdge(bond.to, bond.from, bond.distance, correlation::math::Vector3<double>(-bond.r_x, -bond.r_y, -bond.r_z));
    }
  }
}

} // namespace

bool has_gpu_device() {
  int device_count = 0;
  hipError_t err = hipGetDeviceCount(&device_count);
  return (err == hipSuccess && device_count > 0);
}

void compute_distances_gpu(const correlation::core::Cell &cell, double cutoff_sq,
                           const std::vector<std::vector<double>> &bond_cutoffs_sq,
                           bool ignore_periodic_self_interactions, DistanceTensor &out_distances,
                           correlation::core::NeighborGraph &out_graph) {
  const auto &atoms = cell.atoms();
  const size_t atom_count = atoms.size();
  if (atom_count == 0) {
    return;
  }
  const size_t num_elements = cell.elements().size();
  const auto &lattice = cell.latticeVectors();

  // 1. Recreate CPU Cell-List configuration
  double vol = cell.volume();
  double width_x = vol / correlation::math::norm(correlation::math::cross(lattice[1], lattice[2]));
  double width_y = vol / correlation::math::norm(correlation::math::cross(lattice[0], lattice[2]));
  double width_z = vol / correlation::math::norm(correlation::math::cross(lattice[0], lattice[1]));

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
  std::vector<int> element_ids(atom_count);

  for (size_t i = 0; i < atom_count; ++i) {
    correlation::math::Vector3<double> frac = cell.inverseLatticeVectors() * atoms[i].position();
    double f_x = frac.x() - std::floor(frac.x());
    double f_y = frac.y() - std::floor(frac.y());
    double f_z = frac.z() - std::floor(frac.z());

    wrapped_x[i] = std::fma(f_z, lattice[2].x(), std::fma(f_y, lattice[1].x(), f_x * lattice[0].x()));
    wrapped_y[i] = std::fma(f_z, lattice[2].y(), std::fma(f_y, lattice[1].y(), f_x * lattice[0].y()));
    wrapped_z[i] = std::fma(f_z, lattice[2].z(), std::fma(f_y, lattice[1].z(), f_x * lattice[0].z()));

    int c_x = std::clamp(static_cast<int>(std::floor(f_x * K_x)), 0, K_x - 1);
    int c_y = std::clamp(static_cast<int>(std::floor(f_y * K_y)), 0, K_y - 1);
    int c_z = std::clamp(static_cast<int>(std::floor(f_z * K_z)), 0, K_z - 1);

    int bin_idx = c_x * (K_y * K_z) + c_y * K_z + c_z;
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
    int bin_idx = atom_bin[i];
    bin_indices[insertion_cursors[bin_idx]++] = i;
  }

  // Flat 2D vector for bond cutoffs
  std::vector<double> h_bond_cutoffs_sq = flatten_bond_cutoffs(bond_cutoffs_sq, num_elements);

  // 2. Allocate and copy to GPU device memory
  double *d_wrapped_x = nullptr;
  double *d_wrapped_y = nullptr;
  double *d_wrapped_z = nullptr;
  int *d_element_ids = nullptr;
  int *d_atom_bin = nullptr;
  unsigned long long *d_bin_offsets = nullptr;
  unsigned long long *d_bin_indices = nullptr;
  double *d_bond_cutoffs_sq = nullptr;
  unsigned long long *d_distance_counter = nullptr;
  unsigned long long *d_bond_counter = nullptr;

  hipMalloc(&d_wrapped_x, atom_count * sizeof(double));
  hipMalloc(&d_wrapped_y, atom_count * sizeof(double));
  hipMalloc(&d_wrapped_z, atom_count * sizeof(double));
  hipMalloc(&d_element_ids, atom_count * sizeof(int));
  hipMalloc(&d_atom_bin, atom_count * sizeof(int));
  hipMalloc(&d_bin_offsets, (num_bins + 1) * sizeof(unsigned long long));
  hipMalloc(&d_bin_indices, atom_count * sizeof(unsigned long long));
  hipMalloc(&d_bond_cutoffs_sq, num_elements * num_elements * sizeof(double));
  hipMalloc(&d_distance_counter, sizeof(unsigned long long));
  hipMalloc(&d_bond_counter, sizeof(unsigned long long));

  hipMemcpy(d_wrapped_x, wrapped_x.data(), atom_count * sizeof(double), hipMemcpyHostToDevice);
  hipMemcpy(d_wrapped_y, wrapped_y.data(), atom_count * sizeof(double), hipMemcpyHostToDevice);
  hipMemcpy(d_wrapped_z, wrapped_z.data(), atom_count * sizeof(double), hipMemcpyHostToDevice);
  hipMemcpy(d_element_ids, element_ids.data(), atom_count * sizeof(int), hipMemcpyHostToDevice);
  hipMemcpy(d_atom_bin, atom_bin.data(), atom_count * sizeof(int), hipMemcpyHostToDevice);
  hipMemcpy(d_bin_offsets, bin_offsets.data(), (num_bins + 1) * sizeof(unsigned long long), hipMemcpyHostToDevice);
  hipMemcpy(d_bin_indices, bin_indices.data(), atom_count * sizeof(unsigned long long), hipMemcpyHostToDevice);
  hipMemcpy(d_bond_cutoffs_sq, h_bond_cutoffs_sq.data(), num_elements * num_elements * sizeof(double),
            hipMemcpyHostToDevice);

  // Reset counters to 0
  unsigned long long zero_val = 0;
  hipMemcpy(d_distance_counter, &zero_val, sizeof(unsigned long long), hipMemcpyHostToDevice);
  hipMemcpy(d_bond_counter, &zero_val, sizeof(unsigned long long), hipMemcpyHostToDevice);

  // Setup struct parameters
  GPULattice gpu_lattice{lattice[0].x(), lattice[0].y(), lattice[0].z(), lattice[1].x(), lattice[1].y(),
                         lattice[1].z(), lattice[2].x(), lattice[2].y(), lattice[2].z()};
  GPUSearchGrid gpu_grid{K_x, K_y, K_z, max_dx, max_dy, max_dz};
  GPUAtomData gpu_atoms{d_wrapped_x, d_wrapped_y, d_wrapped_z, d_element_ids, d_atom_bin};
  GPUBinData gpu_bins{d_bin_offsets, d_bin_indices};

  // 3. Launch Pass 1 to count distances and bonds
  int block_size = 256;
  int grid_size = (static_cast<int>(atom_count) + block_size - 1) / block_size;

  hipLaunchKernelGGL(distance_kernel<false>, grid_size, block_size, 0, 0, gpu_atoms, gpu_bins, gpu_lattice, gpu_grid,
                     cutoff_sq, d_bond_cutoffs_sq, static_cast<int>(num_elements), ignore_periodic_self_interactions,
                     static_cast<int>(atom_count), d_distance_counter, d_bond_counter, nullptr, nullptr);
  hipDeviceSynchronize();

  // Copy counts back to host
  unsigned long long h_distance_count = 0;
  unsigned long long h_bond_count = 0;
  hipMemcpy(&h_distance_count, d_distance_counter, sizeof(unsigned long long), hipMemcpyDeviceToHost);
  hipMemcpy(&h_bond_count, d_bond_counter, sizeof(unsigned long long), hipMemcpyDeviceToHost);

  // 4. Allocate outputs on GPU and launch Pass 2 (WriteMode = true)
  GPUDistance *d_distances = nullptr;
  GPUBond *d_bonds = nullptr;

  if (h_distance_count > 0) {
    hipMalloc(&d_distances, h_distance_count * sizeof(GPUDistance));
  }
  if (h_bond_count > 0) {
    hipMalloc(&d_bonds, h_bond_count * sizeof(GPUBond));
  }

  // Reset counters
  hipMemcpy(d_distance_counter, &zero_val, sizeof(unsigned long long), hipMemcpyHostToDevice);
  hipMemcpy(d_bond_counter, &zero_val, sizeof(unsigned long long), hipMemcpyHostToDevice);

  hipLaunchKernelGGL(distance_kernel<true>, grid_size, block_size, 0, 0, gpu_atoms, gpu_bins, gpu_lattice, gpu_grid,
                     cutoff_sq, d_bond_cutoffs_sq, static_cast<int>(num_elements), ignore_periodic_self_interactions,
                     static_cast<int>(atom_count), d_distance_counter, d_bond_counter, d_distances, d_bonds);
  hipDeviceSynchronize();

  // 5. Copy outputs back to host
  std::vector<GPUDistance> host_distances(h_distance_count);
  std::vector<GPUBond> host_bonds(h_bond_count);

  if (h_distance_count > 0) {
    hipMemcpy(host_distances.data(), d_distances, h_distance_count * sizeof(GPUDistance), hipMemcpyDeviceToHost);
  }
  if (h_bond_count > 0) {
    hipMemcpy(host_bonds.data(), d_bonds, h_bond_count * sizeof(GPUBond), hipMemcpyDeviceToHost);
  }

  // Free device memory
  hipFree(d_wrapped_x);
  hipFree(d_wrapped_y);
  hipFree(d_wrapped_z);
  hipFree(d_element_ids);
  hipFree(d_atom_bin);
  hipFree(d_bin_offsets);
  hipFree(d_bin_indices);
  hipFree(d_bond_cutoffs_sq);
  hipFree(d_distance_counter);
  hipFree(d_bond_counter);
  if (d_distances != nullptr) {
    hipFree(d_distances);
  }
  if (d_bonds != nullptr) {
    hipFree(d_bonds);
  }

  // 6. Unpack results on CPU
  unpack_gpu_results(host_distances, host_bonds, element_ids, out_distances, out_graph);
}

} // namespace correlation::calculators::gpu
