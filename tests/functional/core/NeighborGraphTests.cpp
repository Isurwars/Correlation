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
using correlation::real_t;

namespace {
class NeighborGraphFunctionalTests : public ::testing::Test {};

auto buildSimpleCubicGraph(real_t lattice_a) -> NeighborGraph {
  const size_t node_count = 8;
  NeighborGraph graph(node_count);

  const std::vector<Vector3R> coords = {
      {static_cast<real_t>(0.0), static_cast<real_t>(0.0), static_cast<real_t>(0.0)},
      {lattice_a, static_cast<real_t>(0.0), static_cast<real_t>(0.0)},
      {static_cast<real_t>(0.0), lattice_a, static_cast<real_t>(0.0)},
      {static_cast<real_t>(0.0), static_cast<real_t>(0.0), lattice_a},
      {lattice_a, lattice_a, static_cast<real_t>(0.0)},
      {lattice_a, static_cast<real_t>(0.0), lattice_a},
      {static_cast<real_t>(0.0), lattice_a, lattice_a},
      {lattice_a, lattice_a, lattice_a},
  };

  const auto add_symmetric_edge = [&](size_t i_idx, size_t j_idx) {
    Vector3R const r_ij = coords[j_idx] - coords[i_idx];
    Vector3R const r_ji = coords[i_idx] - coords[j_idx];
    graph.addDirectedEdge(static_cast<correlation::core::AtomID>(i_idx),
                          static_cast<correlation::core::AtomID>(j_idx), lattice_a, r_ij);
    graph.addDirectedEdge(static_cast<correlation::core::AtomID>(j_idx),
                          static_cast<correlation::core::AtomID>(i_idx), lattice_a, r_ji);
  };

  // Edges along X
  add_symmetric_edge(0, 1);
  add_symmetric_edge(2, 4);
  add_symmetric_edge(3, 5);
  add_symmetric_edge(6, 7);

  // Edges along Y
  add_symmetric_edge(0, 2);
  add_symmetric_edge(1, 4);
  add_symmetric_edge(3, 6);
  add_symmetric_edge(5, 7);

  // Edges along Z
  add_symmetric_edge(0, 3);
  add_symmetric_edge(1, 5);
  add_symmetric_edge(2, 6);
  add_symmetric_edge(4, 7);

  return graph;
}

void verifySimpleCubicConnectivity(const NeighborGraph &graph) {
  EXPECT_EQ(graph.nodeCount(), 8);
  EXPECT_TRUE(graph.areConnected(AtomIndex{0}, AtomIndex{1}));
  EXPECT_TRUE(graph.areConnected(AtomIndex{0}, AtomIndex{2}));
  EXPECT_TRUE(graph.areConnected(AtomIndex{0}, AtomIndex{3}));
  EXPECT_FALSE(graph.areConnected(AtomIndex{0}, AtomIndex{7}));
}

void verifySimpleCubicNeighbors(const NeighborGraph &graph, real_t lattice_a) {
  const auto &neighbors_0 = graph.getNeighbors(0);
  ASSERT_EQ(neighbors_0.size(), 3);

  for (const auto &neigh : neighbors_0) {
    EXPECT_NEAR(neigh.distance, lattice_a, static_cast<real_t>(1e-6));
    EXPECT_NEAR(norm(neigh.r_ij), lattice_a, static_cast<real_t>(1e-6));
  }
}

void verifySimpleCubicAdjacency(const NeighborGraph &graph) {
  const auto adj_matrix = graph.getDenseAdjacencyMatrix();
  ASSERT_EQ(adj_matrix.size(), 64);

  const std::vector<bool> expected_row0 = {
      false, true, true, true, false, false, false, false,
  };
  for (size_t col = 0; col < 8; ++col) {
    EXPECT_EQ(adj_matrix[col], expected_row0[col]);
  }
}

TEST_F(NeighborGraphFunctionalTests, VerifySimpleCubicTopology) {
  const auto lattice_a = static_cast<real_t>(3.0);
  const NeighborGraph graph = buildSimpleCubicGraph(lattice_a);

  verifySimpleCubicConnectivity(graph);
  verifySimpleCubicNeighbors(graph, lattice_a);
  verifySimpleCubicAdjacency(graph);
}

TEST_F(NeighborGraphFunctionalTests, VerifySelfInteractionsAndParallelEdges) {
  NeighborGraph graph(3);
  const Vector3<real_t> zero{static_cast<real_t>(0.0), static_cast<real_t>(0.0),
                             static_cast<real_t>(0.0)};
  const Vector3<real_t> vec{static_cast<real_t>(1.5), static_cast<real_t>(0.0),
                            static_cast<real_t>(0.0)};

  // Add self-loop to atom 0
  graph.addDirectedEdge(0, 0, static_cast<real_t>(0.0), zero);

  // Add multiple parallel directed bonds between atom 1 and 2 (e.g. multi-path PBC image bonds)
  graph.addDirectedEdge(1, 2, static_cast<real_t>(1.5), vec);
  graph.addDirectedEdge(1, 2, static_cast<real_t>(4.5), vec * static_cast<real_t>(3.0));

  EXPECT_TRUE(graph.areConnected(AtomIndex{0}, AtomIndex{0}));
  EXPECT_TRUE(graph.areConnected(AtomIndex{1}, AtomIndex{2}));

  const auto &neighbors_0 = graph.getNeighbors(0);
  ASSERT_EQ(neighbors_0.size(), 1);
  EXPECT_EQ(neighbors_0[0].index, 0);
  EXPECT_NEAR(neighbors_0[0].distance, static_cast<real_t>(0.0), static_cast<real_t>(1e-6));

  const auto &neighbors_1 = graph.getNeighbors(1);
  ASSERT_EQ(neighbors_1.size(), 2);
  EXPECT_EQ(neighbors_1[0].index, 2);
  EXPECT_NEAR(neighbors_1[0].distance, static_cast<real_t>(1.5), static_cast<real_t>(1e-6));
  EXPECT_EQ(neighbors_1[1].index, 2);
  EXPECT_NEAR(neighbors_1[1].distance, static_cast<real_t>(4.5), static_cast<real_t>(1e-6));
}

} // namespace
} // namespace correlation::testing
