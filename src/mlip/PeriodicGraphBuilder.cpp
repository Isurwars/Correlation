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

/// @brief Flattens thread-local edge buffers into COO-format output arrays.
void flattenEdges(const tbb::enumerable_thread_specific<std::vector<EdgeTuple>> &local_edges, PeriodicGraphData &data) {
  size_t total_edges = 0;
  for (const auto &vec : local_edges) {
    total_edges += vec.size();
  }
  data.edge_count = total_edges;

  data.edge_index_flat.resize(2 * total_edges);
  data.edge_shifts_flat.resize(total_edges * 3);
  data.edge_vectors_flat.resize(total_edges * 3);
  data.edge_distances.resize(total_edges);

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

      ++edge_idx;
    }
  }
}

} // namespace

int64_t PeriodicGraphBuilder::getAtomicNumber(std::string_view symbol) noexcept {
  static const std::unordered_map<std::string_view, int64_t> kSymbolToZ = {
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

  auto itx = kSymbolToZ.find(symbol);
  if (itx != kSymbolToZ.end()) {
    return itx->second;
  }
  return 0;
}

PeriodicGraphData PeriodicGraphBuilder::buildGraph(const correlation::core::Cell &cell, real_t cutoff_radius,
                                                   bool include_self_loops) {
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

  flattenEdges(local_edges, data);

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
    const real_t n_val = static_cast<real_t>(idx_basis + 1);
    const real_t k_n = n_val * k_pi / cutoff_radius;
    real_t bessel_val = static_cast<real_t>(0.0);
    if (distance < static_cast<real_t>(1e-8)) {
      bessel_val = k_n;
    } else {
      bessel_val = std::sin(k_n * distance) / distance;
    }
    basis[idx_basis] = norm_factor * bessel_val * env;
  }

  return basis;
}

std::vector<real_t> PeriodicGraphBuilder::computeGaussianRBF(real_t distance,
                                                             const GaussianRBFConfig &config) {
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

} // namespace correlation::mlip
