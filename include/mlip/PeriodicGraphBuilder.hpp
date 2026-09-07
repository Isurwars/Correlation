/**
 * @file PeriodicGraphBuilder.hpp
 * @brief High-performance periodic atomic neighbor graph builder for GNN inference.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Cell.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"

#include <array>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

namespace correlation::mlip {

/**
 * @struct PeriodicGraphData
 * @brief Extracted graph tensor buffers for GNN model inference.
 */
struct PeriodicGraphData {
  std::vector<real_t> positions_flat;    /**< [N * 3] Cartesian atomic coordinates in Angstroms. */
  std::vector<int64_t> atomic_numbers;   /**< [N] Atomic numbers (Z). */
  std::vector<int64_t> edge_index_flat;  /**< [2 * E] Directed edge indices (row 0: src, row 1: dst). */
  std::vector<real_t> edge_shifts_flat;  /**< [E * 3] Periodic cell displacement integer shift vectors. */
  std::vector<real_t> edge_vectors_flat; /**< [E * 3] Cartesian displacement vectors r_ij = r_j + R*s - r_i. */
  std::vector<real_t> edge_distances;    /**< [E] Euclidean edge distances ||r_ij||. */
  std::vector<real_t> edge_spherical_harmonics_flat; /**< [E * (l_max + 1)^2] Equivariant spherical harmonics features
                                                        (if l_max > 0). */
  std::array<real_t, 9> cell_flat{};                 /**< [3 * 3] Lattice vectors matrix. */
  size_t atom_count{0};                              /**< Total atom count N. */
  size_t edge_count{0};                              /**< Total directed edge count E. */
};

/**
 * @struct GaussianRBFConfig
 * @brief Configuration parameters for Gaussian radial basis function expansion.
 */
struct GaussianRBFConfig {
  real_t start{static_cast<real_t>(0.0)}; /**< Start center distance in Angstroms. */
  real_t stop{static_cast<real_t>(5.0)};  /**< Stop center distance in Angstroms. */
  size_t num_basis{8};                    /**< Number of Gaussian basis centers. */
};

/**
 * @class PeriodicGraphBuilder
 * @brief Constructs periodic neighbor graphs and extracts tensors for GNN evaluation.
 */
class PeriodicGraphBuilder {
public:
  /**
   * @brief Constructs a periodic atomic graph with periodic boundary condition shift vectors.
   * @param[in] cell Simulation cell containing lattice vectors and atomic positions.
   * @param[in] cutoff_radius Cutoff sphere radius in Angstroms (default: 5.0).
   * @param[in] include_self_loops Whether to include zero-displacement self-loops (default: false).
   * @param[in] l_max Maximum spherical harmonics degree (0..3) to embed (default: 0, no harmonics computed).
   * @return Extracted flat PeriodicGraphData buffers.
   */
  [[nodiscard]] static PeriodicGraphData buildGraph(const correlation::core::Cell &cell,
                                                    real_t cutoff_radius = static_cast<real_t>(5.0),
                                                    bool include_self_loops = false, size_t l_max = 0);

  /**
   * @brief Resolves atomic number Z for a chemical element symbol.
   * @param[in] symbol Element symbol string (e.g. "Si", "Fe", "O").
   * @return Atomic number Z (1..118) or 0 if unknown.
   */
  [[nodiscard]] static int64_t getAtomicNumber(std::string_view symbol) noexcept;

  /**
   * @brief Computes a smooth polynomial cutoff envelope f_c(d) decaying to zero at cutoff.
   * @param[in] distance Interatomic distance d.
   * @param[in] cutoff_radius Cutoff radius r_c.
   * @return Envelope value in range [0, 1].
   */
  [[nodiscard]] static real_t computeCutoffEnvelope(real_t distance, real_t cutoff_radius) noexcept;

  /**
   * @brief Computes spherical Bessel radial basis functions with polynomial envelope.
   * @param[in] distance Interatomic distance d.
   * @param[in] cutoff_radius Cutoff radius r_c.
   * @param[in] num_basis Number of radial Bessel basis functions.
   * @return Vector of Bessel basis values of length @p num_basis.
   */
  [[nodiscard]] static std::vector<real_t> computeBesselBasis(real_t distance, real_t cutoff_radius, size_t num_basis);

  /**
   * @brief Computes Gaussian radial basis functions (RBF) for an interatomic distance.
   * @param[in] distance Interatomic distance d.
   * @param[in] config Gaussian RBF configuration parameters.
   * @return Vector of Gaussian RBF values of length @p config.num_basis.
   */
  [[nodiscard]] static std::vector<real_t> computeGaussianRBF(real_t distance, const GaussianRBFConfig &config);

  /**
   * @brief Evaluates real orthonormal spherical harmonics Y_lm(r) up to l_max <= 3 into a pre-allocated span.
   * @param[in] vec 3D displacement vector r.
   * @param[in] l_max Maximum degree (clamped to 3). Total components evaluated: (l_max + 1)^2.
   * @param[out] out Destination span of size at least (l_max + 1)^2.
   */
  static void computeSphericalHarmonics(const correlation::math::Vector3<real_t> &vec, size_t l_max,
                                        std::span<real_t> out) noexcept;

  /**
   * @brief Evaluates real orthonormal spherical harmonics Y_lm(r) up to l_max <= 3 returning a newly allocated vector.
   * @param[in] vec 3D displacement vector r.
   * @param[in] l_max Maximum degree (clamped to 3). Total components evaluated: (l_max + 1)^2.
   * @return Vector of length (l_max + 1)^2 containing real spherical harmonics.
   */
  [[nodiscard]] static std::vector<real_t> computeSphericalHarmonics(const correlation::math::Vector3<real_t> &vec,
                                                                     size_t l_max);
};

} // namespace correlation::mlip
