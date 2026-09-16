/**
 * @file GraphDescriptors.hpp
 * @brief Topological, structural, and spectral graph descriptors for atomic graphs.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/Precision.hpp"
#include "mlip/PeriodicGraphBuilder.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace correlation::mlip {

/**
 * @enum CNALabel
 * @brief Canonical Common Neighbor Analysis structural classification motifs.
 */
enum class CNALabel : std::uint8_t {
  Other = 0, /**< Unclassified or amorphous local environment. */
  FCC = 1,   /**< Face-Centered Cubic coordination motif. */
  HCP = 2,   /**< Hexagonal Close-Packed coordination motif. */
  BCC = 3,   /**< Body-Centered Cubic coordination motif. */
  ICO = 4    /**< Icosahedral coordination motif. */
};

/**
 * @class GraphDescriptors
 * @brief Utility for extracting topological, structural, and spectral descriptors from
 * PeriodicGraphData.
 */
class GraphDescriptors {
public:
  /**
   * @brief Computes per-atom ring statistics embedding using cycle basis detection.
   *
   * For each atom i in [0, N-1] and ring size s in [1, max_size], stores the count
   * of chordless rings of size s that contain atom i.
   *
   * @param[in] graph The periodic neighbor graph tensor buffers.
   * @param[in] max_size Maximum ring size to search for (default: 6).
   * @return Flattened array of size [N * max_size] with per-atom ring counts.
   */
  [[nodiscard]] static std::vector<real_t>
  computeRingStatisticsDescriptor(const PeriodicGraphData &graph, size_t max_size = 6);

  /**
   * @brief Computes per-atom Common Neighbor Analysis (CNA) classification labels.
   *
   * @param[in] graph The periodic neighbor graph tensor buffers.
   * @return Array of size [N] with CNALabel integer values (0=Other, 1=FCC, 2=HCP, 3=BCC, 4=ICO).
   */
  [[nodiscard]] static std::vector<int> computeCNADescriptor(const PeriodicGraphData &graph);

  /**
   * @brief Computes per-atom coordination number embedding.
   *
   * @param[in] graph The periodic neighbor graph tensor buffers.
   * @return Array of size [N] containing the degree of each node.
   */
  [[nodiscard]] static std::vector<real_t>
  computeCoordinationEmbedding(const PeriodicGraphData &graph);

  /**
   * @brief Computes top-k eigenvalues of the graph adjacency matrix.
   *
   * @param[in] graph The periodic neighbor graph tensor buffers.
   * @param[in] k_eigenvalues Number of leading eigenvalues to compute.
   * @return Array of top-k eigenvalues sorted in descending order.
   */
  [[nodiscard]] static std::vector<real_t> computeGraphSpectrum(const PeriodicGraphData &graph,
                                                                size_t k_eigenvalues);

  /**
   * @brief Populates all descriptor fields in PeriodicGraphData in-place.
   *
   * @param[in,out] graph The periodic neighbor graph tensor buffers to enrich.
   * @param[in] max_ring_size Maximum ring size for ring embeddings (default: 6).
   */
  static void populateDescriptors(PeriodicGraphData &graph, size_t max_ring_size = 6);
};

} // namespace correlation::mlip
