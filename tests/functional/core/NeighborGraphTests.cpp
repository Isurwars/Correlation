// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/NeighborGraph.hpp"
#include "math/LinearAlgebra.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {

using namespace correlation::core;
using namespace correlation::math;

namespace {
class NeighborGraphFunctionalTests : public ::testing::Test {};

TEST_F(NeighborGraphFunctionalTests, VerifySimpleCubicTopology) {
  // Construct a neighbor graph representing a Simple Cubic unit cell (8 corners)
  // Indices:
  // 0: (0,0,0)  1: (1,0,0)  2: (0,1,0)  3: (0,0,1)
  // 4: (1,1,0)  5: (1,0,1)  6: (0,1,1)  7: (1,1,1)
  const size_t node_count = 8;
  NeighborGraph graph(node_count);

  const double lattice_a = 3.0; // 3.0 Angstroms
  
  // Define coordinate mapping for convenience
  std::vector<Vector3<double>> coords = {
    {0.0, 0.0, 0.0}, {lattice_a, 0.0, 0.0}, {0.0, lattice_a, 0.0}, {0.0, 0.0, lattice_a},
    {lattice_a, lattice_a, 0.0}, {lattice_a, 0.0, lattice_a}, {0.0, lattice_a, lattice_a}, {lattice_a, lattice_a, lattice_a}
  };

  // Add directed edges for every adjacent pair along grid lines (Simple Cubic nearest neighbors)
  // For SC, every node has neighbors at distance lattice_a
  auto addSymmetricEdge = [&](size_t i, size_t j) {
    Vector3<double> r_ij = coords[j] - coords[i];
    Vector3<double> r_ji = coords[i] - coords[j];
    graph.addDirectedEdge(i, j, lattice_a, r_ij);
    graph.addDirectedEdge(j, i, lattice_a, r_ji);
  };

  // Edges along X
  addSymmetricEdge(0, 1);
  addSymmetricEdge(2, 4);
  addSymmetricEdge(3, 5);
  addSymmetricEdge(6, 7);

  // Edges along Y
  addSymmetricEdge(0, 2);
  addSymmetricEdge(1, 4);
  addSymmetricEdge(3, 6);
  addSymmetricEdge(5, 7);

  // Edges along Z
  addSymmetricEdge(0, 3);
  addSymmetricEdge(1, 5);
  addSymmetricEdge(2, 6);
  addSymmetricEdge(4, 7);

  EXPECT_EQ(graph.nodeCount(), 8);

  // Verify that corner 0 is connected to 1 (X), 2 (Y), and 3 (Z)
  EXPECT_TRUE(graph.areConnected(AtomIndex{0}, AtomIndex{1}));
  EXPECT_TRUE(graph.areConnected(AtomIndex{0}, AtomIndex{2}));
  EXPECT_TRUE(graph.areConnected(AtomIndex{0}, AtomIndex{3}));
  EXPECT_FALSE(graph.areConnected(AtomIndex{0}, AtomIndex{7})); // opposite corner not directly connected

  // Get neighbors of node 0
  const auto &neighbors_0 = graph.getNeighbors(0);
  ASSERT_EQ(neighbors_0.size(), 3);
  
  for (const auto &neigh : neighbors_0) {
    EXPECT_DOUBLE_EQ(neigh.distance, lattice_a);
    // Relative vector norm should match distance
    EXPECT_NEAR(norm(neigh.r_ij), lattice_a, 1e-9);
  }

  // Get dense adjacency matrix
  auto adj_matrix = graph.getDenseAdjacencyMatrix();
  ASSERT_EQ(adj_matrix.size(), node_count * node_count);

  // Assert row 0 entries
  EXPECT_FALSE(adj_matrix[0 * 8 + 0]); // self
  EXPECT_TRUE(adj_matrix[0 * 8 + 1]);  // 0-1
  EXPECT_TRUE(adj_matrix[0 * 8 + 2]);  // 0-2
  EXPECT_TRUE(adj_matrix[0 * 8 + 3]);  // 0-3
  EXPECT_FALSE(adj_matrix[0 * 8 + 4]);
  EXPECT_FALSE(adj_matrix[0 * 8 + 5]);
  EXPECT_FALSE(adj_matrix[0 * 8 + 6]);
  EXPECT_FALSE(adj_matrix[0 * 8 + 7]);
}

TEST_F(NeighborGraphFunctionalTests, VerifySelfInteractionsAndParallelEdges) {
  NeighborGraph graph(3);
  Vector3<double> zero{0.0, 0.0, 0.0};
  Vector3<double> vec{1.5, 0.0, 0.0};

  // Add self-loop to atom 0
  graph.addDirectedEdge(0, 0, 0.0, zero);
  
  // Add multiple parallel directed bonds between atom 1 and 2 (e.g. multi-path PBC image bonds)
  graph.addDirectedEdge(1, 2, 1.5, vec);
  graph.addDirectedEdge(1, 2, 4.5, vec * 3.0);

  EXPECT_TRUE(graph.areConnected(AtomIndex{0}, AtomIndex{0}));
  EXPECT_TRUE(graph.areConnected(AtomIndex{1}, AtomIndex{2}));

  const auto &neighbors_0 = graph.getNeighbors(0);
  ASSERT_EQ(neighbors_0.size(), 1);
  EXPECT_EQ(neighbors_0[0].index, 0);
  EXPECT_DOUBLE_EQ(neighbors_0[0].distance, 0.0);

  const auto &neighbors_1 = graph.getNeighbors(1);
  ASSERT_EQ(neighbors_1.size(), 2);
  EXPECT_EQ(neighbors_1[0].index, 2);
  EXPECT_DOUBLE_EQ(neighbors_1[0].distance, 1.5);
  EXPECT_EQ(neighbors_1[1].index, 2);
  EXPECT_DOUBLE_EQ(neighbors_1[1].distance, 4.5);
}

} // namespace
} // namespace correlation::testing
