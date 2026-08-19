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
#include <string_view>
#include <vector>

namespace correlation::mlip {

/**
 * @struct PeriodicGraphData
 * @brief Extracted graph tensor buffers for GNN model inference.
 */
struct PeriodicGraphData {
  std::vector<real_t> positions_flat;   /**< [N * 3] Cartesian atomic coordinates in Angstroms. */
  std::vector<int64_t> atomic_numbers;  /**< [N] Atomic numbers (Z). */
  std::vector<int64_t> edge_index_flat; /**< [2 * E] Directed edge indices (row 0: src, row 1: dst). */
  std::vector<real_t> edge_shifts_flat; /**< [E * 3] Periodic cell displacement integer shift vectors. */
  std::array<real_t, 9> cell_flat{};    /**< [3 * 3] Lattice vectors matrix. */
  size_t atom_count{0};                 /**< Total atom count N. */
  size_t edge_count{0};                 /**< Total directed edge count E. */
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
   * @return Extracted flat PeriodicGraphData buffers.
   */
  [[nodiscard]] static PeriodicGraphData buildGraph(const correlation::core::Cell &cell,
                                                    real_t cutoff_radius = static_cast<real_t>(5.0),
                                                    bool include_self_loops = false);

  /**
   * @brief Resolves atomic number Z for a chemical element symbol.
   * @param[in] symbol Element symbol string (e.g. "Si", "Fe", "O").
   * @return Atomic number Z (1..118) or 0 if unknown.
   */
  [[nodiscard]] static int64_t getAtomicNumber(std::string_view symbol) noexcept;
};

} // namespace correlation::mlip
