/**
 * @file NeighborGraph.hpp
 * @brief Compressed Sparse Row (CSR) neighbor graph of atomic bonds and adjacency.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Atom.hpp"
#include "math/LinearAlgebra.hpp"

#include <atomic>
#include <mutex>
#include <span>
#include <vector>

namespace correlation::core {

/**
 * @brief Strong type wrapper representing an atom index.
 */
struct AtomIndex {
  size_t id; ///< Internal 0-based atom index.
};

/**
 * @brief Represents a neighbor atom in a neighbor list.
 */
struct Neighbor {
  AtomID index{0};            ///< Index of the neighbor in the atom list
  real_t distance{0.0};       ///< Distance to the neighbor
  math::Vector3<real_t> r_ij; ///< Vector from central atom to neighbor
};

/**
 * @brief Represents a directed edge connecting two atoms in the graph.
 */
struct Edge {
  AtomID source{};            ///< Index of the source atom
  AtomID target{};            ///< Index of the target atom
  real_t distance{};          ///< Separation distance (Angstrom)
  math::Vector3<real_t> r_ij; ///< Relative position vector from source to target
};

/**
 * @class NeighborGraph
 * @brief Compressed Sparse Row (CSR) based graph representing atomic connectivity.
 *
 * Stores bonds and neighbor interactions as directed edges with deterministic
 * row sorting by (target_index, distance). Eliminates per-row vector allocations
 * and provides lightweight index-only views for high-throughput traversals.
 */
class NeighborGraph {
public:
  /** @brief Default constructor. */
  explicit NeighborGraph() = default;

  /**
   * @brief Constructs a graph with a fixed number of nodes.
   * @param node_count The total number of atoms/nodes in the graph.
   */
  explicit NeighborGraph(size_t node_count);

  /**
   * @brief Copy constructor.
   * @param other The graph instance to copy from.
   */
  NeighborGraph(const NeighborGraph &other);

  /**
   * @brief Move constructor.
   * @param other The graph instance to move from.
   */
  NeighborGraph(NeighborGraph &&other) noexcept;

  /**
   * @brief Copy assignment operator.
   * @param other The graph instance to copy from.
   * @return Reference to this graph.
   */
  NeighborGraph &operator=(const NeighborGraph &other);

  /**
   * @brief Move assignment operator.
   * @param other The graph instance to move from.
   * @return Reference to this graph.
   */
  NeighborGraph &operator=(NeighborGraph &&other) noexcept;

  /** @brief Destructor. */
  ~NeighborGraph() = default;

  /**
   * @brief Adds a directed edge between two atoms.
   * @param edge The directed edge details containing source, target, distance, and vector.
   */
  void addDirectedEdge(const Edge &edge);

  /**
   * @brief Adds a directed edge between two atoms using individual parameters.
   * @param source Index of the source atom.
   * @param target Index of the target atom.
   * @param distance Separation distance (Angstrom).
   * @param r_ij Relative position vector from source to target.
   */
  void addDirectedEdge(AtomID source, AtomID target, real_t distance,
                       const math::Vector3<real_t> &r_ij) {
    addDirectedEdge(Edge{
        .source = source,
        .target = target,
        .distance = distance,
        .r_ij = r_ij,
    });
  }

  /**
   * @brief Finalizes the graph into Compressed Sparse Row (CSR) storage.
   *
   * Sorts all neighbor rows by (index, distance) deterministically.
   * Automatically invoked before queries if edges have been modified.
   */
  void finalize() const;

  /**
   * @brief Directly assigns Compressed Sparse Row (CSR) storage.
   * @param node_count Number of nodes.
   * @param offsets Row offsets array of size node_count + 1.
   * @param neighbors Flat neighbors array of size offsets.back().
   * @param targets Flat target indices array of size offsets.back().
   */
  void assignCSR(size_t node_count, std::vector<size_t> offsets, std::vector<Neighbor> neighbors,
                 std::vector<AtomID> targets);

  /**
   * @brief Retrieves the neighbor list for a specific atom.
   * @param atom_index Index of the atom.
   * @return Contiguous span of Neighbor objects.
   */
  [[nodiscard]] std::span<const Neighbor> getNeighbors(size_t atom_index) const;

  /**
   * @brief Retrieves a lightweight index-only view of neighbors for a specific atom.
   * @param atom_index Index of the atom.
   * @return Contiguous span of AtomID indices.
   */
  [[nodiscard]] std::span<const AtomID> getNeighborIndices(size_t atom_index) const;

  /**
   * @brief Checks if two atoms are connected by an edge using binary search.
   * @param first_atom First atom index.
   * @param second_atom Second atom index.
   * @return True if second_atom is a neighbor of first_atom.
   */
  [[nodiscard]] bool areConnected(AtomIndex first_atom, AtomIndex second_atom) const;

  /**
   * @brief Generates a flattened dense adjacency matrix.
   * @return A vector of booleans of size N*N.
   */
  [[nodiscard]] std::vector<bool> getDenseAdjacencyMatrix() const;

  /** @return Total number of nodes in the graph. */
  [[nodiscard]] size_t nodeCount() const noexcept;

  /** @return Total number of directed edges in the graph. */
  [[nodiscard]] size_t edgeCount() const noexcept;

  /** @return CSR row offsets span of size nodeCount() + 1. */
  [[nodiscard]] std::span<const size_t> offsets() const noexcept;

  /** @return CSR flat target indices span of size edgeCount(). */
  [[nodiscard]] std::span<const AtomID> targets() const noexcept;

  /** @return CSR flat neighbors span of size edgeCount(). */
  [[nodiscard]] std::span<const Neighbor> neighbors() const noexcept;

private:
  NeighborGraph(const std::scoped_lock<std::mutex> & /*lock*/, const NeighborGraph &other);
  NeighborGraph(const std::scoped_lock<std::mutex> & /*lock*/, NeighborGraph &other) noexcept;

  void ensureFinalized() const;

  size_t node_count_{0};
  mutable std::atomic<bool> is_finalized_{true};
  mutable std::mutex finalize_mutex_;
  mutable std::vector<Edge> unfinalized_edges_;
  mutable std::vector<size_t> offsets_{0};
  mutable std::vector<Neighbor> neighbors_;
  mutable std::vector<AtomID> targets_;
};

} // namespace correlation::core
