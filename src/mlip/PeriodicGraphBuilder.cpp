/**
 * @file PeriodicGraphBuilder.cpp
 * @brief Implementation of PeriodicGraphBuilder for atomic GNN inference.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "mlip/PeriodicGraphBuilder.hpp"
#include "math/LinearAlgebra.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string_view>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <unordered_map>

namespace correlation::mlip {

namespace {

struct EdgeTuple {
  int64_t src;
  int64_t dst;
  real_t shift_x;
  real_t shift_y;
  real_t shift_z;
  real_t vec_x;
  real_t vec_y;
  real_t vec_z;
  real_t distance;
};

/// @brief Candidate neighbor pair with periodic image shift and squared distance.
struct CandidateEdge {
  size_t atom_i;
  size_t atom_j;
  int nx;
  int ny;
  int nz;
  real_t vec_x;
  real_t vec_y;
  real_t vec_z;
  real_t dist_sq;
  real_t distance;
};

// Periodic search grid boundaries
struct SearchBoundaries {
  int max_nx{0};
  int max_ny{0};
  int max_nz{0};
};

[[nodiscard]] SearchBoundaries computeSearchBoundaries(const correlation::core::Cell &cell, real_t cutoff) {
  SearchBoundaries bounds;
  if (cell.volume() <= 0.0) {
    return bounds;
  }

  const auto &matrix = cell.latticeVectors();
  const auto lat_a = correlation::math::Vector3<real_t>{matrix[0][0], matrix[0][1], matrix[0][2]};
  const auto lat_b = correlation::math::Vector3<real_t>{matrix[1][0], matrix[1][1], matrix[1][2]};
  const auto lat_c = correlation::math::Vector3<real_t>{matrix[2][0], matrix[2][1], matrix[2][2]};

  const auto b_cross_c = correlation::math::cross(lat_b, lat_c);
  const auto c_cross_a = correlation::math::cross(lat_c, lat_a);
  const auto a_cross_b = correlation::math::cross(lat_a, lat_b);

  const real_t norm_bc = correlation::math::norm(b_cross_c);
  const real_t norm_ca = correlation::math::norm(c_cross_a);
  const real_t norm_ab = correlation::math::norm(a_cross_b);

  if (norm_bc > 1e-6) {
    const real_t h_a = cell.volume() / norm_bc;
    bounds.max_nx = static_cast<int>(std::ceil(cutoff / h_a));
  }
  if (norm_ca > 1e-6) {
    const real_t h_b = cell.volume() / norm_ca;
    bounds.max_ny = static_cast<int>(std::ceil(cutoff / h_b));
  }
  if (norm_ab > 1e-6) {
    const real_t h_c = cell.volume() / norm_ab;
    bounds.max_nz = static_cast<int>(std::ceil(cutoff / h_c));
  }

  // Safety cap on search boundaries to prevent runaway memory on pathological thin cells
  bounds.max_nx = std::clamp(bounds.max_nx, 0, 8);
  bounds.max_ny = std::clamp(bounds.max_ny, 0, 8);
  bounds.max_nz = std::clamp(bounds.max_nz, 0, 8);

  return bounds;
}

/// @brief Tests whether a candidate edge is a self-interaction that should be skipped.
[[nodiscard]] bool isSelfInteraction(const CandidateEdge &edge, bool include_self_loops) noexcept {
  if (include_self_loops) {
    return false;
  }
  const bool is_same_atom = (edge.atom_i == edge.atom_j);
  const bool is_zero_shift = (edge.nx == 0 && edge.ny == 0 && edge.nz == 0);
  return is_same_atom && is_zero_shift && edge.dist_sq < 1e-12;
}

/// @brief Collects all periodic neighbor edges for atom @p atom_i within @p cutoff_sq.
void collectEdgesForAtom(size_t atom_i, const std::vector<correlation::core::Atom> &atoms,
                         const correlation::math::Vector3<real_t> &lat_a,
                         const correlation::math::Vector3<real_t> &lat_b,
                         const correlation::math::Vector3<real_t> &lat_c, const SearchBoundaries &bounds,
                         real_t cutoff_sq, bool include_self_loops, std::vector<EdgeTuple> &edges_out) {
  const auto &pos_i = atoms[atom_i].position();
  const size_t number_atoms = atoms.size();

  for (size_t j = 0; j < number_atoms; ++j) {
    const auto &pos_j = atoms[j].position();

    for (int nx = -bounds.max_nx; nx <= bounds.max_nx; ++nx) {
      for (int ny = -bounds.max_ny; ny <= bounds.max_ny; ++ny) {
        for (int nz = -bounds.max_nz; nz <= bounds.max_nz; ++nz) {
          const auto shift_vec =
              lat_a * static_cast<real_t>(nx) + lat_b * static_cast<real_t>(ny) + lat_c * static_cast<real_t>(nz);
          const auto r_ij = pos_j + shift_vec - pos_i;
          const real_t dist_sq = correlation::math::dot(r_ij, r_ij);

          if (dist_sq > cutoff_sq) {
            continue;
          }

          const real_t dist = std::sqrt(dist_sq);
          const CandidateEdge candidate{.atom_i = atom_i,
                                        .atom_j = j,
                                        .nx = nx,
                                        .ny = ny,
                                        .nz = nz,
                                        .vec_x = r_ij.x(),
                                        .vec_y = r_ij.y(),
                                        .vec_z = r_ij.z(),
                                        .dist_sq = dist_sq,
                                        .distance = dist};
          if (isSelfInteraction(candidate, include_self_loops)) {
            continue;
          }

          edges_out.push_back(EdgeTuple{
              .src = static_cast<int64_t>(candidate.atom_i),
              .dst = static_cast<int64_t>(candidate.atom_j),
              .shift_x = static_cast<real_t>(candidate.nx),
              .shift_y = static_cast<real_t>(candidate.ny),
              .shift_z = static_cast<real_t>(candidate.nz),
              .vec_x = candidate.vec_x,
              .vec_y = candidate.vec_y,
              .vec_z = candidate.vec_z,
              .distance = candidate.distance,
          });
        }
      }
    }
  }
}

// Precomputed constants for real orthonormal spherical harmonics Y_lm up to l=3
constexpr real_t k_pi_val = std::numbers::pi_v<real_t>;

// l=0: sqrt(1 / (4 * pi))
constexpr real_t k_y00 = static_cast<real_t>(0.5) * std::numbers::inv_sqrtpi_v<real_t>;

[[nodiscard]] constexpr real_t computeY0() noexcept { return k_y00; }

[[nodiscard]] std::array<real_t, 3> computeY1(real_t dir_x, real_t dir_y, real_t dir_z, real_t inv_r) noexcept {
  const real_t factor_y1 = std::sqrt(static_cast<real_t>(3.0) / (static_cast<real_t>(4.0) * k_pi_val));
  return {
      factor_y1 * dir_y * inv_r,
      factor_y1 * dir_z * inv_r,
      factor_y1 * dir_x * inv_r,
  };
}

[[nodiscard]] std::array<real_t, 5> computeY2(real_t dir_x, real_t dir_y, real_t dir_z, real_t inv_r2) noexcept {
  const real_t factor_m21 = static_cast<real_t>(0.5) * std::sqrt(static_cast<real_t>(15.0) / k_pi_val);
  const real_t factor_0 = static_cast<real_t>(0.25) * std::sqrt(static_cast<real_t>(5.0) / k_pi_val);
  const real_t factor_2 = static_cast<real_t>(0.25) * std::sqrt(static_cast<real_t>(15.0) / k_pi_val);
  const real_t sq_x = dir_x * dir_x;
  const real_t sq_y = dir_y * dir_y;
  const real_t sq_z = dir_z * dir_z;
  return {
      factor_m21 * dir_x * dir_y * inv_r2,
      factor_m21 * dir_y * dir_z * inv_r2,
      factor_0 * (static_cast<real_t>(2.0) * sq_z - sq_x - sq_y) * inv_r2,
      factor_m21 * dir_x * dir_z * inv_r2,
      factor_2 * (sq_x - sq_y) * inv_r2,
  };
}

[[nodiscard]] std::array<real_t, 7> computeY3(real_t dir_x, real_t dir_y, real_t dir_z, real_t inv_r3) noexcept {
  const real_t factor_3 =
      static_cast<real_t>(0.25) * std::sqrt(static_cast<real_t>(35.0) / (static_cast<real_t>(2.0) * k_pi_val));
  const real_t factor_m2 = static_cast<real_t>(0.5) * std::sqrt(static_cast<real_t>(105.0) / k_pi_val);
  const real_t factor_1 =
      static_cast<real_t>(0.25) * std::sqrt(static_cast<real_t>(21.0) / (static_cast<real_t>(2.0) * k_pi_val));
  const real_t factor_0 = static_cast<real_t>(0.25) * std::sqrt(static_cast<real_t>(7.0) / k_pi_val);
  const real_t factor_2 = static_cast<real_t>(0.25) * std::sqrt(static_cast<real_t>(105.0) / k_pi_val);

  const real_t sq_x = dir_x * dir_x;
  const real_t sq_y = dir_y * dir_y;
  const real_t sq_z = dir_z * dir_z;
  return {
      factor_3 * dir_y * (static_cast<real_t>(3.0) * sq_x - sq_y) * inv_r3,
      factor_m2 * dir_x * dir_y * dir_z * inv_r3,
      factor_1 * dir_y * (static_cast<real_t>(4.0) * sq_z - sq_x - sq_y) * inv_r3,
      factor_0 * dir_z *
          (static_cast<real_t>(2.0) * sq_z - static_cast<real_t>(3.0) * sq_x - static_cast<real_t>(3.0) * sq_y) *
          inv_r3,
      factor_1 * dir_x * (static_cast<real_t>(4.0) * sq_z - sq_x - sq_y) * inv_r3,
      factor_2 * dir_z * (sq_x - sq_y) * inv_r3,
      factor_3 * dir_x * (sq_x - static_cast<real_t>(3.0) * sq_y) * inv_r3,
  };
}

[[nodiscard]] constexpr real_t computeOrbY0() noexcept { return static_cast<real_t>(1.0); }

[[nodiscard]] std::array<real_t, 3> computeOrbY1(real_t dir_x, real_t dir_y, real_t dir_z) noexcept {
  const real_t factor_y1 = std::numbers::sqrt3_v<real_t>;
  return {factor_y1 * dir_x, factor_y1 * dir_y, factor_y1 * dir_z};
}

[[nodiscard]] std::array<real_t, 5> computeOrbY2(real_t dir_x, real_t dir_y, real_t dir_z) noexcept {
  const real_t sqrt3 = std::numbers::sqrt3_v<real_t>;
  const real_t sqrt5 = std::sqrt(static_cast<real_t>(5.0));

  const real_t sh_2_0 = sqrt3 * dir_x * dir_z;
  const real_t sh_2_1 = sqrt3 * dir_x * dir_y;
  const real_t sq_y = dir_y * dir_y;
  const real_t sq_xz = dir_x * dir_x + dir_z * dir_z;
  const real_t sh_2_2 = sq_y - static_cast<real_t>(0.5) * sq_xz;
  const real_t sh_2_3 = sqrt3 * dir_y * dir_z;
  const real_t sh_2_4 = sqrt3 * static_cast<real_t>(0.5) * (dir_z * dir_z - dir_x * dir_x);

  return {
      sqrt5 * sh_2_0, sqrt5 * sh_2_1, sqrt5 * sh_2_2, sqrt5 * sh_2_3, sqrt5 * sh_2_4,
  };
}

[[nodiscard]] std::array<real_t, 7> computeOrbY3(real_t dir_x, real_t dir_y, real_t dir_z) noexcept {
  const real_t sqrt3 = std::numbers::sqrt3_v<real_t>;
  const real_t sqrt7 = std::sqrt(static_cast<real_t>(7.0));

  const real_t sh_2_0 = sqrt3 * dir_x * dir_z;
  const real_t sq_y = dir_y * dir_y;
  const real_t sq_xz = dir_x * dir_x + dir_z * dir_z;
  const real_t sh_2_4 = sqrt3 * static_cast<real_t>(0.5) * (dir_z * dir_z - dir_x * dir_x);

  const real_t factor_30_36 = std::sqrt(static_cast<real_t>(5.0) / static_cast<real_t>(6.0));
  const real_t factor_31_35 = std::sqrt(static_cast<real_t>(5.0));
  const real_t factor_32_34 = std::sqrt(static_cast<real_t>(3.0) / static_cast<real_t>(8.0));

  const real_t sh_3_0 = factor_30_36 * (sh_2_0 * dir_z + sh_2_4 * dir_x);
  const real_t sh_3_1 = factor_31_35 * sh_2_0 * dir_y;
  const real_t sh_3_2 = factor_32_34 * (static_cast<real_t>(4.0) * sq_y - sq_xz) * dir_x;
  const real_t sh_3_3 =
      static_cast<real_t>(0.5) * dir_y * (static_cast<real_t>(2.0) * sq_y - static_cast<real_t>(3.0) * sq_xz);
  const real_t sh_3_4 = factor_32_34 * (static_cast<real_t>(4.0) * sq_y - sq_xz) * dir_z;
  const real_t sh_3_5 = factor_31_35 * sh_2_4 * dir_y;
  const real_t sh_3_6 = factor_30_36 * (sh_2_4 * dir_z - sh_2_0 * dir_x);

  return {
      sqrt7 * sh_3_0, sqrt7 * sh_3_1, sqrt7 * sh_3_2, sqrt7 * sh_3_3, sqrt7 * sh_3_4, sqrt7 * sh_3_5, sqrt7 * sh_3_6,
  };
}

struct OrbUnitCutoffResult {
  real_t unit_x{0.0};
  real_t unit_y{0.0};
  real_t unit_z{0.0};
  real_t cutoff_val{0.0};
};

[[nodiscard]] inline OrbUnitCutoffResult computeOrbUnitAndCutoff(const correlation::math::Vector3<real_t> &vec,
                                                                 real_t edge_dist, real_t cutoff_radius) noexcept {
  OrbUnitCutoffResult res;
  if (edge_dist > static_cast<real_t>(1e-10)) {
    const real_t inv_dist = static_cast<real_t>(1.0) / edge_dist;
    res.unit_x = vec.x() * inv_dist;
    res.unit_y = vec.y() * inv_dist;
    res.unit_z = vec.z() * inv_dist;
  }
  res.cutoff_val = PeriodicGraphBuilder::computeOrbCutoffEnvelope(edge_dist, cutoff_radius);
  return res;
}

struct OrbEdgeContext {
  real_t cutoff_radius{6.0};
  real_t bessel_prefactor{0.0};
  size_t num_rbf{8};
  size_t num_sh{16};
  size_t effective_l_max{3};
  bool compute_orb_features{true};
  std::span<const real_t> bessel_weights;
};

struct OrbLocalBuffers {
  std::span<real_t> sh_buffer;
  std::span<real_t> rbf_buffer;
};

inline void evaluateOrbEdge(size_t idx_edge, const PeriodicGraphData &data_in, const OrbEdgeContext &ctx,
                            PeriodicGraphData &data_out, const OrbLocalBuffers &buffers) noexcept {
  const real_t vec_x = data_in.edge_vectors_flat[idx_edge * 3 + 0];
  const real_t vec_y = data_in.edge_vectors_flat[idx_edge * 3 + 1];
  const real_t vec_z = data_in.edge_vectors_flat[idx_edge * 3 + 2];
  const real_t edge_dist = data_in.edge_distances[idx_edge];

  const auto geom =
      computeOrbUnitAndCutoff(correlation::math::Vector3<real_t>{vec_x, vec_y, vec_z}, edge_dist, ctx.cutoff_radius);
  data_out.edge_unit_vectors_flat[idx_edge * 3 + 0] = geom.unit_x;
  data_out.edge_unit_vectors_flat[idx_edge * 3 + 1] = geom.unit_y;
  data_out.edge_unit_vectors_flat[idx_edge * 3 + 2] = geom.unit_z;
  data_out.edge_cutoff_envelope_flat[idx_edge] = geom.cutoff_val;

  for (size_t idx_rbf = 0; idx_rbf < ctx.num_rbf; ++idx_rbf) {
    const real_t weight_n = ctx.bessel_weights[idx_rbf];
    real_t rbf_val = 0.0;
    if (edge_dist < static_cast<real_t>(1e-8)) {
      rbf_val = ctx.bessel_prefactor * weight_n;
    } else {
      rbf_val = ctx.bessel_prefactor * std::sin(weight_n * edge_dist) / edge_dist;
    }
    buffers.rbf_buffer[idx_rbf] = rbf_val;
    data_out.edge_radial_basis_flat[idx_edge * ctx.num_rbf + idx_rbf] = rbf_val;
  }

  PeriodicGraphBuilder::computeOrbSphericalHarmonics(
      correlation::math::Vector3<real_t>{geom.unit_x, geom.unit_y, geom.unit_z}, ctx.effective_l_max,
      buffers.sh_buffer);
  for (size_t idx_sh = 0; idx_sh < ctx.num_sh; ++idx_sh) {
    data_out.edge_spherical_harmonics_flat[idx_edge * ctx.num_sh + idx_sh] = buffers.sh_buffer[idx_sh];
  }

  if (ctx.compute_orb_features) {
    const size_t orb_feat_dim = ctx.num_rbf * ctx.num_sh;
    for (size_t idx_rbf = 0; idx_rbf < ctx.num_rbf; ++idx_rbf) {
      const real_t scaled_rbf = geom.cutoff_val * buffers.rbf_buffer[idx_rbf];
      for (size_t idx_sh = 0; idx_sh < ctx.num_sh; ++idx_sh) {
        data_out.edge_orb_features_flat[idx_edge * orb_feat_dim + idx_rbf * ctx.num_sh + idx_sh] =
            scaled_rbf * buffers.sh_buffer[idx_sh];
      }
    }
  }
}

/// @brief Flattens thread-local edge buffers into COO-format output arrays.
void flattenEdges(const tbb::enumerable_thread_specific<std::vector<EdgeTuple>> &local_edges, PeriodicGraphData &data,
                  size_t l_max) {
  size_t total_edges = 0;
  for (const auto &vec : local_edges) {
    total_edges += vec.size();
  }
  data.edge_count = total_edges;

  data.edge_index_flat.resize(2 * total_edges);
  data.edge_shifts_flat.resize(total_edges * 3);
  data.edge_vectors_flat.resize(total_edges * 3);
  data.edge_distances.resize(total_edges);

  const size_t effective_l_max = std::min(l_max, size_t{3});
  const size_t num_sh = (effective_l_max + 1) * (effective_l_max + 1);
  if (l_max > 0) {
    data.edge_spherical_harmonics_flat.resize(total_edges * num_sh);
  } else {
    data.edge_spherical_harmonics_flat.clear();
  }

  size_t edge_idx = 0;
  for (const auto &vec : local_edges) {
    for (const auto &edge : vec) {
      data.edge_index_flat[0 * total_edges + edge_idx] = edge.src;
      data.edge_index_flat[1 * total_edges + edge_idx] = edge.dst;

      data.edge_shifts_flat[edge_idx * 3 + 0] = edge.shift_x;
      data.edge_shifts_flat[edge_idx * 3 + 1] = edge.shift_y;
      data.edge_shifts_flat[edge_idx * 3 + 2] = edge.shift_z;

      data.edge_vectors_flat[edge_idx * 3 + 0] = edge.vec_x;
      data.edge_vectors_flat[edge_idx * 3 + 1] = edge.vec_y;
      data.edge_vectors_flat[edge_idx * 3 + 2] = edge.vec_z;

      data.edge_distances[edge_idx] = edge.distance;

      if (l_max > 0) {
        const std::span<real_t> sh_span(data.edge_spherical_harmonics_flat.data() + edge_idx * num_sh, num_sh);
        PeriodicGraphBuilder::computeSphericalHarmonics(
            correlation::math::Vector3<real_t>{edge.vec_x, edge.vec_y, edge.vec_z}, effective_l_max, sh_span);
      }

      ++edge_idx;
    }
  }
}

} // namespace

int64_t PeriodicGraphBuilder::getAtomicNumber(std::string_view symbol) noexcept {
  static const std::unordered_map<std::string_view, int64_t> k_symbol_to_z = {
      {"H", 1},    {"He", 2},   {"Li", 3},   {"Be", 4},   {"B", 5},    {"C", 6},    {"N", 7},    {"O", 8},
      {"F", 9},    {"Ne", 10},  {"Na", 11},  {"Mg", 12},  {"Al", 13},  {"Si", 14},  {"P", 15},   {"S", 16},
      {"Cl", 17},  {"Ar", 18},  {"K", 19},   {"Ca", 20},  {"Sc", 21},  {"Ti", 22},  {"V", 23},   {"Cr", 24},
      {"Mn", 25},  {"Fe", 26},  {"Co", 27},  {"Ni", 28},  {"Cu", 29},  {"Zn", 30},  {"Ga", 31},  {"Ge", 32},
      {"As", 33},  {"Se", 34},  {"Br", 35},  {"Kr", 36},  {"Rb", 37},  {"Sr", 38},  {"Y", 39},   {"Zr", 40},
      {"Nb", 41},  {"Mo", 42},  {"Tc", 43},  {"Ru", 44},  {"Rh", 45},  {"Pd", 46},  {"Ag", 47},  {"Cd", 48},
      {"In", 49},  {"Sn", 50},  {"Sb", 51},  {"Te", 52},  {"I", 53},   {"Xe", 54},  {"Cs", 55},  {"Ba", 56},
      {"La", 57},  {"Ce", 58},  {"Pr", 59},  {"Nd", 60},  {"Pm", 61},  {"Sm", 62},  {"Eu", 63},  {"Gd", 64},
      {"Tb", 65},  {"Dy", 66},  {"Ho", 67},  {"Er", 68},  {"Tm", 69},  {"Yb", 70},  {"Lu", 71},  {"Hf", 72},
      {"Ta", 73},  {"W", 74},   {"Re", 75},  {"Os", 76},  {"Ir", 77},  {"Pt", 78},  {"Au", 79},  {"Hg", 80},
      {"Tl", 81},  {"Pb", 82},  {"Bi", 83},  {"Po", 84},  {"At", 85},  {"Rn", 86},  {"Fr", 87},  {"Ra", 88},
      {"Ac", 89},  {"Th", 90},  {"Pa", 91},  {"U", 92},   {"Np", 93},  {"Pu", 94},  {"Am", 95},  {"Cm", 96},
      {"Bk", 97},  {"Cf", 98},  {"Es", 99},  {"Fm", 100}, {"Md", 101}, {"No", 102}, {"Lr", 103}, {"Rf", 104},
      {"Db", 105}, {"Sg", 106}, {"Bh", 107}, {"Hs", 108}, {"Mt", 109}, {"Ds", 110}, {"Rg", 111}, {"Cn", 112},
      {"Nh", 113}, {"Fl", 114}, {"Mc", 115}, {"Lv", 116}, {"Ts", 117}, {"Og", 118}};

  auto itx = k_symbol_to_z.find(symbol);
  if (itx != k_symbol_to_z.end()) {
    return itx->second;
  }
  return 0;
}

PeriodicGraphData PeriodicGraphBuilder::buildGraph(const correlation::core::Cell &cell, real_t cutoff_radius,
                                                   bool include_self_loops, size_t l_max) {
  PeriodicGraphData data;
  const auto &atoms = cell.atoms();
  const size_t number_atoms = atoms.size();
  data.atom_count = number_atoms;

  if (number_atoms == 0) {
    return data;
  }

  // Populate positions and atomic numbers
  data.positions_flat.resize(number_atoms * 3);
  data.atomic_numbers.resize(number_atoms);

  for (size_t i = 0; i < number_atoms; ++i) {
    const auto &pos = atoms[i].position();
    data.positions_flat[i * 3 + 0] = pos.x();
    data.positions_flat[i * 3 + 1] = pos.y();
    data.positions_flat[i * 3 + 2] = pos.z();
    data.atomic_numbers[i] = getAtomicNumber(atoms[i].element().symbol);
  }

  // Populate cell lattice matrix
  const auto &lat = cell.latticeVectors();
  data.cell_flat = {lat[0][0], lat[0][1], lat[0][2], lat[1][0], lat[1][1], lat[1][2], lat[2][0], lat[2][1], lat[2][2]};

  const auto lat_a = correlation::math::Vector3<real_t>{lat[0][0], lat[0][1], lat[0][2]};
  const auto lat_b = correlation::math::Vector3<real_t>{lat[1][0], lat[1][1], lat[1][2]};
  const auto lat_c = correlation::math::Vector3<real_t>{lat[2][0], lat[2][1], lat[2][2]};

  const real_t cutoff_sq = cutoff_radius * cutoff_radius;
  const auto bounds = computeSearchBoundaries(cell, cutoff_radius);

  // Thread-local edge accumulator
  tbb::enumerable_thread_specific<std::vector<EdgeTuple>> local_edges;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, number_atoms), [&](const tbb::blocked_range<size_t> &range) {
    auto &my_edges = local_edges.local();
    for (size_t i = range.begin(); i < range.end(); ++i) {
      collectEdgesForAtom(i, atoms, lat_a, lat_b, lat_c, bounds, cutoff_sq, include_self_loops, my_edges);
    }
  });

  flattenEdges(local_edges, data, l_max);

  return data;
}

real_t PeriodicGraphBuilder::computeCutoffEnvelope(real_t distance, real_t cutoff_radius) noexcept {
  if (cutoff_radius <= static_cast<real_t>(0.0)) {
    return static_cast<real_t>(0.0);
  }
  const real_t scaled_distance = distance / cutoff_radius;
  if (scaled_distance >= static_cast<real_t>(1.0)) {
    return static_cast<real_t>(0.0);
  }
  if (scaled_distance <= static_cast<real_t>(0.0)) {
    return static_cast<real_t>(1.0);
  }
  // Smooth 5th-order polynomial cutoff envelope: 1 - 6u^5 + 15u^4 - 10u^3
  const real_t dist_cubed = scaled_distance * scaled_distance * scaled_distance;
  const real_t dist_fourth = dist_cubed * scaled_distance;
  const real_t dist_fifth = dist_fourth * scaled_distance;
  return static_cast<real_t>(1.0) - static_cast<real_t>(6.0) * dist_fifth + static_cast<real_t>(15.0) * dist_fourth -
         static_cast<real_t>(10.0) * dist_cubed;
}

std::vector<real_t> PeriodicGraphBuilder::computeBesselBasis(real_t distance, real_t cutoff_radius, size_t num_basis) {
  std::vector<real_t> basis(num_basis, static_cast<real_t>(0.0));
  if (num_basis == 0 || cutoff_radius <= static_cast<real_t>(0.0)) {
    return basis;
  }

  const real_t env = computeCutoffEnvelope(distance, cutoff_radius);
  if (env <= static_cast<real_t>(0.0)) {
    return basis;
  }

  constexpr real_t k_pi = std::numbers::pi_v<real_t>;
  const real_t norm_factor = std::sqrt(static_cast<real_t>(2.0) / cutoff_radius);

  for (size_t idx_basis = 0; idx_basis < num_basis; ++idx_basis) {
    const auto n_val = static_cast<real_t>(idx_basis + 1);
    const real_t k_n = n_val * k_pi / cutoff_radius;
    auto bessel_val = static_cast<real_t>(0.0);
    if (distance < static_cast<real_t>(1e-8)) {
      bessel_val = k_n;
    } else {
      bessel_val = std::sin(k_n * distance) / distance;
    }
    basis[idx_basis] = norm_factor * bessel_val * env;
  }

  return basis;
}

std::vector<real_t> PeriodicGraphBuilder::computeGaussianRBF(real_t distance, const GaussianRBFConfig &config) {
  std::vector<real_t> basis(config.num_basis, static_cast<real_t>(0.0));
  if (config.num_basis == 0) {
    return basis;
  }
  if (config.num_basis == 1) {
    const real_t center = (config.start + config.stop) * static_cast<real_t>(0.5);
    const real_t diff = distance - center;
    basis[0] = std::exp(-diff * diff);
    return basis;
  }

  const real_t step = (config.stop - config.start) / static_cast<real_t>(config.num_basis - 1);
  const real_t gamma = static_cast<real_t>(0.5) / (step * step);

  for (size_t idx_basis = 0; idx_basis < config.num_basis; ++idx_basis) {
    const real_t center = config.start + static_cast<real_t>(idx_basis) * step;
    const real_t diff = distance - center;
    basis[idx_basis] = std::exp(-gamma * diff * diff);
  }

  return basis;
}

void PeriodicGraphBuilder::computeSphericalHarmonics(const correlation::math::Vector3<real_t> &vec, size_t l_max,
                                                     std::span<real_t> out) noexcept {
  const size_t effective_l_max = std::min(l_max, size_t{3});
  const size_t required_size = (effective_l_max + 1) * (effective_l_max + 1);
  if (out.size() < required_size) {
    return;
  }

  const real_t pos_x = vec.x();
  const real_t pos_y = vec.y();
  const real_t pos_z = vec.z();
  const real_t rad_sq = pos_x * pos_x + pos_y * pos_y + pos_z * pos_z;
  const real_t rad_dist = std::sqrt(rad_sq);

  out[0] = computeY0();

  if (rad_dist < static_cast<real_t>(1e-12)) {
    for (size_t out_idx = 1; out_idx < required_size; ++out_idx) {
      out[out_idx] = static_cast<real_t>(0.0);
    }
    return;
  }

  const real_t inv_r = static_cast<real_t>(1.0) / rad_dist;
  const real_t inv_r2 = inv_r * inv_r;

  if (effective_l_max >= 1) {
    const auto y1_vals = computeY1(pos_x, pos_y, pos_z, inv_r);
    out[1] = y1_vals[0];
    out[2] = y1_vals[1];
    out[3] = y1_vals[2];
  }

  if (effective_l_max >= 2) {
    const auto y2_vals = computeY2(pos_x, pos_y, pos_z, inv_r2);
    out[4] = y2_vals[0];
    out[5] = y2_vals[1];
    out[6] = y2_vals[2];
    out[7] = y2_vals[3];
    out[8] = y2_vals[4];
  }

  if (effective_l_max >= 3) {
    const real_t inv_r3 = inv_r2 * inv_r;
    const auto y3_vals = computeY3(pos_x, pos_y, pos_z, inv_r3);
    out[9] = y3_vals[0];
    out[10] = y3_vals[1];
    out[11] = y3_vals[2];
    out[12] = y3_vals[3];
    out[13] = y3_vals[4];
    out[14] = y3_vals[5];
    out[15] = y3_vals[6];
  }
}

std::vector<real_t> PeriodicGraphBuilder::computeSphericalHarmonics(const correlation::math::Vector3<real_t> &vec,
                                                                    size_t l_max) {
  const size_t effective_l_max = std::min(l_max, size_t{3});
  const size_t required_size = (effective_l_max + 1) * (effective_l_max + 1);
  std::vector<real_t> out(required_size, static_cast<real_t>(0.0));
  computeSphericalHarmonics(vec, effective_l_max, out);
  return out;
}

void PeriodicGraphBuilder::computeOrbSphericalHarmonics(const correlation::math::Vector3<real_t> &unit_vec,
                                                        size_t l_max, std::span<real_t> out) noexcept {
  const size_t effective_l_max = std::min(l_max, size_t{3});
  const size_t required_size = (effective_l_max + 1) * (effective_l_max + 1);
  if (out.size() < required_size) {
    return;
  }

  const real_t pos_x = unit_vec.x();
  const real_t pos_y = unit_vec.y();
  const real_t pos_z = unit_vec.z();
  const real_t rad_sq = pos_x * pos_x + pos_y * pos_y + pos_z * pos_z;

  out[0] = computeOrbY0();

  if (rad_sq < static_cast<real_t>(1e-12)) {
    for (size_t out_idx = 1; out_idx < required_size; ++out_idx) {
      out[out_idx] = static_cast<real_t>(0.0);
    }
    return;
  }

  const real_t inv_r = static_cast<real_t>(1.0) / std::sqrt(rad_sq);
  const real_t norm_x = pos_x * inv_r;
  const real_t norm_y = pos_y * inv_r;
  const real_t norm_z = pos_z * inv_r;

  if (effective_l_max >= 1) {
    const auto y1_vals = computeOrbY1(norm_x, norm_y, norm_z);
    out[1] = y1_vals[0];
    out[2] = y1_vals[1];
    out[3] = y1_vals[2];
  }

  if (effective_l_max >= 2) {
    const auto y2_vals = computeOrbY2(norm_x, norm_y, norm_z);
    out[4] = y2_vals[0];
    out[5] = y2_vals[1];
    out[6] = y2_vals[2];
    out[7] = y2_vals[3];
    out[8] = y2_vals[4];
  }

  if (effective_l_max >= 3) {
    const auto y3_vals = computeOrbY3(norm_x, norm_y, norm_z);
    out[9] = y3_vals[0];
    out[10] = y3_vals[1];
    out[11] = y3_vals[2];
    out[12] = y3_vals[3];
    out[13] = y3_vals[4];
    out[14] = y3_vals[5];
    out[15] = y3_vals[6];
  }
}

std::vector<real_t>
PeriodicGraphBuilder::computeOrbSphericalHarmonics(const correlation::math::Vector3<real_t> &unit_vec, size_t l_max) {
  const size_t effective_l_max = std::min(l_max, size_t{3});
  const size_t required_size = (effective_l_max + 1) * (effective_l_max + 1);
  std::vector<real_t> out(required_size, static_cast<real_t>(0.0));
  computeOrbSphericalHarmonics(unit_vec, effective_l_max, out);
  return out;
}

real_t PeriodicGraphBuilder::computeOrbCutoffEnvelope(real_t distance, real_t cutoff_radius) noexcept {
  if (cutoff_radius <= static_cast<real_t>(0.0) || distance >= cutoff_radius) {
    return static_cast<real_t>(0.0);
  }
  if (distance <= static_cast<real_t>(0.0)) {
    return static_cast<real_t>(1.0);
  }
  const real_t ratio_u = distance / cutoff_radius;
  const real_t ratio_u2 = ratio_u * ratio_u;
  const real_t ratio_u4 = ratio_u2 * ratio_u2;
  const real_t ratio_u5 = ratio_u4 * ratio_u;
  const real_t ratio_u6 = ratio_u5 * ratio_u;
  return static_cast<real_t>(1.0) - static_cast<real_t>(15.0) * ratio_u4 + static_cast<real_t>(24.0) * ratio_u5 -
         static_cast<real_t>(10.0) * ratio_u6;
}

std::vector<real_t> PeriodicGraphBuilder::computeOrbBesselBasis(real_t distance, const OrbDescriptorConfig &config) {
  const size_t num_basis = config.num_rbf;
  const real_t cutoff_radius = config.r_max;
  std::vector<real_t> basis(num_basis, static_cast<real_t>(0.0));
  if (num_basis == 0 || cutoff_radius <= static_cast<real_t>(0.0)) {
    return basis;
  }

  constexpr real_t k_pi = std::numbers::pi_v<real_t>;
  const real_t prefactor = std::sqrt(static_cast<real_t>(2.0) / cutoff_radius);

  for (size_t idx_basis = 0; idx_basis < num_basis; ++idx_basis) {
    const auto n_val = static_cast<real_t>(idx_basis + 1);
    const real_t weight_n = n_val * k_pi / cutoff_radius;
    if (distance < static_cast<real_t>(1e-8)) {
      basis[idx_basis] = prefactor * weight_n;
    } else {
      basis[idx_basis] = prefactor * std::sin(weight_n * distance) / distance;
    }
  }

  return basis;
}

PeriodicGraphData PeriodicGraphBuilder::buildOrbGraph(const correlation::core::Cell &cell,
                                                      const OrbDescriptorConfig &config) {
  const real_t cutoff = config.r_max;
  const bool include_self_loops = config.include_self_loops;
  const size_t effective_l_max = std::min(config.l_max, size_t{3});

  PeriodicGraphData data = buildGraph(cell, cutoff, include_self_loops, 0);
  const size_t total_edges = data.edge_count;

  const size_t num_rbf = config.num_rbf;
  const size_t num_sh = (effective_l_max + 1) * (effective_l_max + 1);

  data.edge_unit_vectors_flat.resize(total_edges * 3);
  data.edge_radial_basis_flat.resize(total_edges * num_rbf);
  data.edge_spherical_harmonics_flat.resize(total_edges * num_sh);
  data.edge_cutoff_envelope_flat.resize(total_edges);

  const size_t orb_feat_dim = num_rbf * num_sh;
  if (config.compute_orb_features) {
    data.edge_orb_features_flat.resize(total_edges * orb_feat_dim);
  } else {
    data.edge_orb_features_flat.clear();
  }

  constexpr real_t k_pi = std::numbers::pi_v<real_t>;
  const real_t cutoff_radius = config.r_max;
  const real_t bessel_prefactor = (cutoff_radius > static_cast<real_t>(0.0))
                                      ? std::sqrt(static_cast<real_t>(2.0) / cutoff_radius)
                                      : static_cast<real_t>(0.0);

  std::vector<real_t> bessel_weights(num_rbf);
  for (size_t idx_weight = 0; idx_weight < num_rbf; ++idx_weight) {
    bessel_weights[idx_weight] = static_cast<real_t>(idx_weight + 1) * k_pi / cutoff_radius;
  }

  const OrbEdgeContext ctx{
      .cutoff_radius = cutoff_radius,
      .bessel_prefactor = bessel_prefactor,
      .num_rbf = num_rbf,
      .num_sh = num_sh,
      .effective_l_max = effective_l_max,
      .compute_orb_features = config.compute_orb_features,
      .bessel_weights = bessel_weights,
  };

  tbb::parallel_for(tbb::blocked_range<size_t>(0, total_edges), [&](const tbb::blocked_range<size_t> &range) {
    std::vector<real_t> local_sh_buf(num_sh);
    std::vector<real_t> local_rbf_buf(num_rbf);
    const std::span<real_t> local_sh(local_sh_buf.data(), num_sh);
    const std::span<real_t> local_rbf(local_rbf_buf.data(), num_rbf);

    for (size_t idx_edge = range.begin(); idx_edge < range.end(); ++idx_edge) {
      evaluateOrbEdge(idx_edge, data, ctx, data, OrbLocalBuffers{.sh_buffer = local_sh, .rbf_buffer = local_rbf});
    }
  });

  return data;
}

} // namespace correlation::mlip
