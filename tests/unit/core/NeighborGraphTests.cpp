// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/NeighborGraph.hpp"
#include "math/LinearAlgebra.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {
namespace {

using namespace correlation::core;
using namespace correlation::math;

class NeighborGraphTests : public ::testing::Test {};

TEST_F(NeighborGraphTests, DefaultConstructorInitializesEmpty) {
  const NeighborGraph graph{};
  EXPECT_EQ(graph.nodeCount(), 0);
  EXPECT_TRUE(graph.getDenseAdjacencyMatrix().empty());
}

TEST_F(NeighborGraphTests, ParameterConstructorInitializesNodeCount) {
  const NeighborGraph graph(10);
  EXPECT_EQ(graph.nodeCount(), 10);

  auto matrix = graph.getDenseAdjacencyMatrix();
  EXPECT_EQ(matrix.size(), 100);
  for (bool val : matrix) {
    EXPECT_FALSE(val);
  }
}

TEST_F(NeighborGraphTests, AddDirectedEdgeAndGetNeighbors) {
  NeighborGraph graph(5);
  Vector3<real_t> r_ij{1.0, 2.0, 3.0};
  graph.addDirectedEdge(0, 1, 3.74, r_ij);

  EXPECT_TRUE(graph.areConnected(AtomIndex{0}, AtomIndex{1}));
  EXPECT_FALSE(graph.areConnected(AtomIndex{1}, AtomIndex{0}));

  const auto &neighbors = graph.getNeighbors(0);
  ASSERT_EQ(neighbors.size(), 1);
  EXPECT_EQ(neighbors[0].index, 1);
  EXPECT_THAT(neighbors[0].distance, ::correlation::testing::IsRealEq(3.74));
  EXPECT_THAT(neighbors[0].r_ij.x(), ::correlation::testing::IsRealEq(r_ij.x()));
  EXPECT_THAT(neighbors[0].r_ij.y(), ::correlation::testing::IsRealEq(r_ij.y()));
  EXPECT_THAT(neighbors[0].r_ij.z(), ::correlation::testing::IsRealEq(r_ij.z()));
}

TEST_F(NeighborGraphTests, OutOfBoundsEdgeAdditionDoesNotCrash) {
  NeighborGraph graph(3);
  Vector3<real_t> r_ij{0.0, 0.0, 0.0};

  // from index >= size is safely ignored (returns immediately)
  EXPECT_NO_THROW(graph.addDirectedEdge(3, 1, 1.0, r_ij));
  EXPECT_FALSE(graph.areConnected(AtomIndex{3}, AtomIndex{1}));

  // to index >= size is added to the adjacency list but filtered out in the dense matrix
  EXPECT_NO_THROW(graph.addDirectedEdge(0, 5, 1.0, r_ij));
  EXPECT_TRUE(graph.areConnected(AtomIndex{0}, AtomIndex{5}));

  auto matrix = graph.getDenseAdjacencyMatrix();
  ASSERT_EQ(matrix.size(), 9);
  for (bool val : matrix) {
    EXPECT_FALSE(val);
  }
}

TEST_F(NeighborGraphTests, OutOfBoundsGetNeighborsReturnsEmpty) {
  const NeighborGraph graph(3);

  // Index out of bounds should return an empty neighbor vector
  const auto &neighbors = graph.getNeighbors(3);
  EXPECT_TRUE(neighbors.empty());

  const auto &neighbors_large = graph.getNeighbors(100);
  EXPECT_TRUE(neighbors_large.empty());
}

TEST_F(NeighborGraphTests, OutOfBoundsAreConnectedReturnsFalse) {
  NeighborGraph graph(3);
  Vector3<real_t> r_ij{0.0, 0.0, 0.0};
  graph.addDirectedEdge(0, 1, 1.0, r_ij);

  EXPECT_FALSE(graph.areConnected(AtomIndex{3}, AtomIndex{1}));
  EXPECT_FALSE(graph.areConnected(AtomIndex{0}, AtomIndex{3}));
  EXPECT_FALSE(graph.areConnected(AtomIndex{5}, AtomIndex{5}));
}

TEST_F(NeighborGraphTests, DenseAdjacencyMatrixMapping) {
  NeighborGraph graph(3);
  Vector3<real_t> r_ij{0.0, 0.0, 0.0};
  graph.addDirectedEdge(0, 1, 1.0, r_ij);
  graph.addDirectedEdge(1, 2, 1.0, r_ij);
  graph.addDirectedEdge(2, 0, 1.0, r_ij);

  auto matrix = graph.getDenseAdjacencyMatrix();
  ASSERT_EQ(matrix.size(), 9);

  // Row 0: connected to 1
  EXPECT_FALSE(matrix[0]); // 0 -> 0
  EXPECT_TRUE(matrix[1]);  // 0 -> 1
  EXPECT_FALSE(matrix[2]); // 0 -> 2

  // Row 1: connected to 2
  EXPECT_FALSE(matrix[3]); // 1 -> 0
  EXPECT_FALSE(matrix[4]); // 1 -> 1
  EXPECT_TRUE(matrix[5]);  // 1 -> 2

  // Row 2: connected to 0
  EXPECT_TRUE(matrix[6]);  // 2 -> 0
  EXPECT_FALSE(matrix[7]); // 2 -> 1
  EXPECT_FALSE(matrix[8]); // 2 -> 2
}

// --- Extreme / Edge-Case Tests ---

TEST_F(NeighborGraphTests, SelfLoopIsAllowed) {
  NeighborGraph graph(3);
  Vector3<real_t> r_ij{0.0, 0.0, 0.0};

  // Self-loop: atom 0 -> atom 0
  graph.addDirectedEdge(0, 0, 0.0, r_ij);

  EXPECT_TRUE(graph.areConnected(AtomIndex{0}, AtomIndex{0}));
  const auto &neighbors = graph.getNeighbors(0);
  ASSERT_EQ(neighbors.size(), 1);
  EXPECT_EQ(neighbors[0].index, 0);

  // Should appear in the dense matrix diagonal
  auto matrix = graph.getDenseAdjacencyMatrix();
  EXPECT_TRUE(matrix[0]); // (0,0) = true
}

TEST_F(NeighborGraphTests, DuplicateEdgesAreBothStored) {
  NeighborGraph graph(3);
  Vector3<real_t> r_ij1{1.0, 0.0, 0.0};
  Vector3<real_t> r_ij2{2.0, 0.0, 0.0};

  // Add same edge (0 -> 1) twice with different distances
  graph.addDirectedEdge(0, 1, 1.5, r_ij1);
  graph.addDirectedEdge(0, 1, 3.0, r_ij2);

  // Both edges should be stored
  const auto &neighbors = graph.getNeighbors(0);
  ASSERT_EQ(neighbors.size(), 2);
  EXPECT_EQ(neighbors[0].index, 1);
  EXPECT_EQ(neighbors[1].index, 1);
  EXPECT_DOUBLE_EQ(neighbors[0].distance, 1.5);
  EXPECT_DOUBLE_EQ(neighbors[1].distance, 3.0);

  // Dense matrix still shows them as connected
  EXPECT_TRUE(graph.areConnected(AtomIndex{0}, AtomIndex{1}));
}

TEST_F(NeighborGraphTests, LargeGraphPerformance) {
  // Stress test: large fully-connected graph should not crash
  const size_t num_nodes = 500;
  NeighborGraph graph(num_nodes);
  Vector3<real_t> r_ij{1.0, 0.0, 0.0};

  // Connect every pair (i,j) where i < j as a directed edge i -> j
  for (size_t i = 0; i < num_nodes; ++i) {
    for (size_t j = i + 1; j < num_nodes && j < i + 5; ++j) { // Limit to 5 neighbors each
      graph.addDirectedEdge(i, j, 1.0, r_ij);
    }
  }

  EXPECT_EQ(graph.nodeCount(), num_nodes);

  // Verify first node's neighbors (connects to j=1,2,3,4 = 4 neighbors)
  const auto &neighbors = graph.getNeighbors(0);
  EXPECT_EQ(neighbors.size(), 4);

  // Dense matrix should be NxN
  auto matrix = graph.getDenseAdjacencyMatrix();
  EXPECT_EQ(matrix.size(), num_nodes * num_nodes);
}

} // namespace
} // namespace correlation::testing
