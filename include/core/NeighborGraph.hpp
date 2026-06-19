/**
 * @file NeighborGraph.hpp
 * @brief Neighbor graph representation of atomic bonds and adjacency.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Atom.hpp"
#include "math/LinearAlgebra.hpp"

#include <vector>

namespace correlation::core {

struct AtomIndex {
  size_t id;
};

/**
 * @brief Represents a neighbor atom in a neighbor list.
 */
struct Neighbor {
  AtomID index{0};            ///< Index of the neighbor in the atom list
  double distance{0.0};       ///< Distance to the neighbor
  math::Vector3<double> r_ij; ///< Vector from central atom to neighbor
};

/**
 * @class NeighborGraph
 * @brief Adjacency-list based graph representing atomic connectivity and proximity.
 *
 * This class stores bonds and neighbor interactions as directed edges.
 * It is primarily used by calculators that require topological information,
 * such as Ring Statistics (RD) or Steinhardt Order Parameters.
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
   * @brief Adds a directed edge between two atoms.
   *
   * @param source Index of the source atom.
   * @param target Index of the target atom.
   * @param distance Separation distance (Angstrom).
   * @param r_ij Relative position vector (target - source).
   */
  void addDirectedEdge(size_t source, size_t target, double distance, const math::Vector3<double> &r_ij); // NOLINT(bugprone-easily-swappable-parameters)

  /**
   * @brief Retrieves the neighbor list for a specific atom.
   * @param atom_index Index of the atom.
   * @return Constant reference to the vector of Neighbor objects.
   */
  [[nodiscard]] const std::vector<Neighbor> &getNeighbors(size_t atom_index) const;

  /**
   * @brief Checks if two atoms are connected by an edge.
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
  [[nodiscard]] size_t nodeCount() const;

private:
  std::vector<std::vector<Neighbor>> adj_list_; ///< Internal adjacency list.
};

} // namespace correlation::core
