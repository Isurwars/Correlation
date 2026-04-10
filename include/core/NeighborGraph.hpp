/**
 * @file NeighborGraph.hpp
 * @brief Neighbor graph representation of atomic bonds and adjacency.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "core/Atom.hpp"
#include "math/LinearAlgebra.hpp"

#include <vector>

namespace correlation::core {

/**
 * @brief Represents a neighbor atom in a neighbor list.
 */
struct Neighbor {
  AtomID index;               ///< Index of the neighbor in the atom list
  double distance;            ///< Distance to the neighbor
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
   * @param from Index of the source atom.
   * @param to Index of the target atom.
   * @param distance Separation distance (Angstrom).
   * @param r_ij Relative position vector (to - from).
   */
  void addDirectedEdge(size_t from, size_t to, double distance,
                       const math::Vector3<double> &r_ij);

  /**
   * @brief Retrieves the neighbor list for a specific atom.
   * @param atom_index Index of the atom.
   * @return Constant reference to the vector of Neighbor objects.
   */
  [[nodiscard]] const std::vector<Neighbor> &
  getNeighbors(size_t atom_index) const;

  /**
   * @brief Checks if two atoms are connected by an edge.
   * @param i First atom index.
   * @param j Second atom index.
   * @return True if j is a neighbor of i.
   */
  [[nodiscard]] bool areConnected(size_t i, size_t j) const;

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
