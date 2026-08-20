/**
 * @file StructureAnalyzer.hpp
 * @brief Structural analysis utilities for neighbor detection and bond
 * topology.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include <atomic>
#include <mutex>
#include <vector>

namespace correlation::analysis {

/**
 * @class StructureAnalyzer
 * @brief Computes and stores pairwise distances and bond angles in tensors
 * indexed by element type.
 *
 * This class efficiently finds all atom pairs within a cutoff radius and all
 * unique bond angles, storing the results in multi-dimensional vectors
 * (tensors) suitable for calculating partial distribution functions.
 * Angle and Dihedral tensors are computed lazily on first access.
 */
class StructureAnalyzer {
public:
  /** @brief Tensor for storing raw distance histogram bins [element1][element2][bin_index]. */
  using RawHistogramTensor = std::vector<std::vector<std::vector<real_t>>>;
  /** @brief Tensor for storing bond angles [center][outer1][outer2][angle_index]. */
  using AngleTensor = std::vector<std::vector<std::vector<std::vector<real_t>>>>;
  /** @brief Tensor for storing dihedrals [e1][e2][e3][e4][dihedral_index]. */
  using DihedralTensor = std::vector<std::vector<std::vector<std::vector<std::vector<real_t>>>>>;

  /** @name Constructors & Destructors */
  ///@{
  explicit StructureAnalyzer(const correlation::core::Cell &cell, real_t cutoff, BondCutoffMatrix bond_cutoffs,
                             bool ignore_periodic_self_interactions = true, real_t r_bin_width = 0.02);

  ~StructureAnalyzer() = default;
  StructureAnalyzer(const StructureAnalyzer &) = delete;
  StructureAnalyzer &operator=(const StructureAnalyzer &) = delete;
  StructureAnalyzer(StructureAnalyzer &&other) noexcept;
  StructureAnalyzer &operator=(StructureAnalyzer &&other) noexcept;

  ///@}

  /** @name Accessors */
  ///@{
  /**
   * @brief Gets a multi-dimensional tensor containing raw distance histogram bins.
   * @return The raw histogram tensor `[element1][element2][bin_index]`.
   */
  [[nodiscard]] const RawHistogramTensor &rawHistograms() const { return raw_histograms_; }

  /**
   * @brief Gets a multi-dimensional tensor containing bond angles.
   *        Computed lazily on first access in a thread-safe manner.
   * @return The angle tensor
   * `[center_element][outer_element1][outer_element2][angle_index]`.
   */
  [[nodiscard]] const AngleTensor &angles() const;

  /**
   * @brief Gets a multi-dimensional tensor containing dihedral angles.
   *        Computed lazily on first access in a thread-safe manner.
   * @return The dihedral tensor
   * `[element1][element2][element3][element4][dihedral_index]`.
   */
  [[nodiscard]] const DihedralTensor &dihedrals() const;

  /**
   * @brief Gets the corresponding neighbor graph capturing topological
   * connections.
   * @return Constant reference to the correlation::core::NeighborGraph object.
   */
  [[nodiscard]] const correlation::core::NeighborGraph &neighborGraph() const { return neighbor_graph_; }

  /**
   * @brief Returns true if periodic self-interactions are ignored.
   * @return True if periodic self-interactions are ignored.
   */
  [[nodiscard]] bool getIgnorePeriodicSelfInteractions() const { return ignore_periodic_self_interactions_; }

  ///@}

private:
  void ensureAnglesComputed() const;
  void ensureDihedralsComputed() const;

  correlation::core::Cell cell_;  ///< Reference to the current periodic cell.
  real_t cutoff_sq_;              ///< Squared cutoff for efficiency.
  BondCutoffMatrix bond_cutoffs_; ///< Internal bond cutoffs matrix.

  bool ignore_periodic_self_interactions_; ///< Interaction guard.

  correlation::core::NeighborGraph neighbor_graph_; ///< Graph of topological bonds.
  RawHistogramTensor raw_histograms_;               ///< Raw distance histogram storage.
  mutable AngleTensor angle_tensor_;                ///< Bond angle storage (lazy).
  mutable DihedralTensor dihedral_tensor_;          ///< Dihedral angle storage (lazy).

  mutable std::mutex compute_mutex_;                    ///< Mutex guarding lazy computation.
  mutable std::atomic<bool> angles_computed_{false};    ///< Flag indicating if angles are computed.
  mutable std::atomic<bool> dihedrals_computed_{false}; ///< Flag indicating if dihedrals are computed.
};

} // namespace correlation::analysis
