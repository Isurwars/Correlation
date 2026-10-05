// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/spatial/MotifFinder.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include <gtest/gtest.h>

namespace correlation::calculators {
namespace {
class MotifFinderTests : public ::testing::Test {
public:
  correlation::core::Cell cell;
  correlation::core::NeighborGraph graph;

  void SetUp() override {
    cell = correlation::core::Cell({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0});
  }
};
} // namespace

TEST_F(MotifFinderTests, DetectsSingleTriangle) {
  graph = correlation::core::NeighborGraph(3);

  // 0-1, 1-2, 2-0
  graph.addDirectedEdge(0, 1, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(1, 0, 1.0, {-1.0, 0.0, 0.0});

  graph.addDirectedEdge(1, 2, 1.0, {0.0, 1.0, 0.0});
  graph.addDirectedEdge(2, 1, 1.0, {0.0, -1.0, 0.0});

  graph.addDirectedEdge(2, 0, 1.0, {-1.0, -1.0, 0.0});
  graph.addDirectedEdge(0, 2, 1.0, {1.0, 1.0, 0.0});

  auto rings = MotifFinder::findRings(graph, 6);

  EXPECT_EQ(rings[3], 1); // Exact 1 triangle
  EXPECT_EQ(rings.count(4), 0);
  EXPECT_EQ(rings.count(5), 0);
  EXPECT_EQ(rings.count(6), 0);
}

TEST_F(MotifFinderTests, DetectsSingleSquare) {
  graph = correlation::core::NeighborGraph(4);

  // 0-1, 1-2, 2-3, 3-0
  graph.addDirectedEdge(0, 1, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(1, 0, 1.0, {-1.0, 0.0, 0.0});

  graph.addDirectedEdge(1, 2, 1.0, {0.0, 1.0, 0.0});
  graph.addDirectedEdge(2, 1, 1.0, {0.0, -1.0, 0.0});

  graph.addDirectedEdge(2, 3, 1.0, {-1.0, 0.0, 0.0});
  graph.addDirectedEdge(3, 2, 1.0, {1.0, 0.0, 0.0});

  graph.addDirectedEdge(3, 0, 1.0, {0.0, -1.0, 0.0});
  graph.addDirectedEdge(0, 3, 1.0, {0.0, 1.0, 0.0});

  auto rings = MotifFinder::findRings(graph, 6);

  EXPECT_EQ(rings.count(3), 0);
  EXPECT_EQ(rings[4], 1); // Exact 1 square
  EXPECT_EQ(rings.count(5), 0);
  EXPECT_EQ(rings.count(6), 0);
}

// --- Extreme / Edge-Case Tests ---

TEST_F(MotifFinderTests, EmptyGraphReturnsNoRings) {
  // A graph with no edges at all
  graph = correlation::core::NeighborGraph(5);

  auto rings = MotifFinder::findRings(graph, 6);

  // No edges means no rings of any size
  for (int size = 3; size <= 6; ++size) {
    EXPECT_EQ(rings.count(size), 0) << "Expected no rings of size " << size;
  }
}

TEST_F(MotifFinderTests, IsolatedNodesReturnsNoRings) {
  // Graph with some edges but no closed loops
  graph = correlation::core::NeighborGraph(4);

  // Linear chain: 0-1-2-3 (no cycle)
  graph.addDirectedEdge(0, 1, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(1, 0, 1.0, {-1.0, 0.0, 0.0});
  graph.addDirectedEdge(1, 2, 1.0, {0.0, 1.0, 0.0});
  graph.addDirectedEdge(2, 1, 1.0, {0.0, -1.0, 0.0});
  graph.addDirectedEdge(2, 3, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(3, 2, 1.0, {-1.0, 0.0, 0.0});

  auto rings = MotifFinder::findRings(graph, 6);

  for (int size = 3; size <= 6; ++size) {
    EXPECT_EQ(rings.count(size), 0) << "Expected no rings of size " << size;
  }
}

TEST_F(MotifFinderTests, MaxRingSizeExcludesLargerRings) {
  // Create a square (ring of size 4) but set max_size = 3
  graph = correlation::core::NeighborGraph(4);

  graph.addDirectedEdge(0, 1, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(1, 0, 1.0, {-1.0, 0.0, 0.0});
  graph.addDirectedEdge(1, 2, 1.0, {0.0, 1.0, 0.0});
  graph.addDirectedEdge(2, 1, 1.0, {0.0, -1.0, 0.0});
  graph.addDirectedEdge(2, 3, 1.0, {-1.0, 0.0, 0.0});
  graph.addDirectedEdge(3, 2, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(3, 0, 1.0, {0.0, -1.0, 0.0});
  graph.addDirectedEdge(0, 3, 1.0, {0.0, 1.0, 0.0});

  // max_size = 3 should NOT find the size-4 ring
  auto rings = MotifFinder::findRings(graph, 3);

  EXPECT_EQ(rings.count(3), 0);
  EXPECT_EQ(rings.count(4), 0); // Should be excluded by max_size
}

TEST_F(MotifFinderTests, ExtractCyclesExtractsExactTopology) {
  // Setup a graph containing one 3-ring (0-1-2) and one 4-ring (3-4-5-6)
  graph = correlation::core::NeighborGraph(7);

  // 3-ring: 0-1-2
  graph.addDirectedEdge(0, 1, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(1, 0, 1.0, {-1.0, 0.0, 0.0});
  graph.addDirectedEdge(1, 2, 1.0, {0.0, 1.0, 0.0});
  graph.addDirectedEdge(2, 1, 1.0, {0.0, -1.0, 0.0});
  graph.addDirectedEdge(2, 0, 1.0, {-1.0, -1.0, 0.0});
  graph.addDirectedEdge(0, 2, 1.0, {1.0, 1.0, 0.0});

  // 4-ring: 3-4-5-6
  graph.addDirectedEdge(3, 4, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(4, 3, 1.0, {-1.0, 0.0, 0.0});
  graph.addDirectedEdge(4, 5, 1.0, {0.0, 1.0, 0.0});
  graph.addDirectedEdge(5, 4, 1.0, {0.0, -1.0, 0.0});
  graph.addDirectedEdge(5, 6, 1.0, {-1.0, 0.0, 0.0});
  graph.addDirectedEdge(6, 5, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(6, 3, 1.0, {0.0, -1.0, 0.0});
  graph.addDirectedEdge(3, 6, 1.0, {0.0, 1.0, 0.0});

  auto cycles_3 = MotifFinder::extractCycles(graph, 3);
  ASSERT_EQ(cycles_3.size(), 1);
  EXPECT_EQ(cycles_3[0].size(), 3);

  auto cycles_4 = MotifFinder::extractCycles(graph, 4);
  ASSERT_EQ(cycles_4.size(), 1);
  EXPECT_EQ(cycles_4[0].size(), 4);

  auto cycles_5 = MotifFinder::extractCycles(graph, 5);
  EXPECT_TRUE(cycles_5.empty());
}

TEST_F(MotifFinderTests, MaxSizeLessThanThreeReturnsEmpty) {
  graph = correlation::core::NeighborGraph(3);
  graph.addDirectedEdge(0, 1, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(1, 0, 1.0, {-1.0, 0.0, 0.0});
  graph.addDirectedEdge(1, 2, 1.0, {0.0, 1.0, 0.0});
  graph.addDirectedEdge(2, 1, 1.0, {0.0, -1.0, 0.0});
  graph.addDirectedEdge(2, 0, 1.0, {-1.0, -1.0, 0.0});
  graph.addDirectedEdge(0, 2, 1.0, {1.0, 1.0, 0.0});

  auto rings = MotifFinder::findRings(graph, 2);
  EXPECT_TRUE(rings.empty());

  auto cycles = MotifFinder::extractCycles(graph, 2);
  EXPECT_TRUE(cycles.empty());
}

TEST_F(MotifFinderTests, FranzblauDetectsSimpleHexagon) {
  graph = correlation::core::NeighborGraph(6);
  // Simple hexagon: 0-1-2-3-4-5-0
  for (size_t i = 0; i < 6; ++i) {
    size_t next = (i + 1) % 6;
    graph.addDirectedEdge(i, next, 1.0, {1.0, 0.0, 0.0});
    graph.addDirectedEdge(next, i, 1.0, {-1.0, 0.0, 0.0});
  }

  auto king_rings = MotifFinder::findRings(graph, 6, RingType::King);
  EXPECT_EQ(king_rings[6], 1U);

  auto franzblau_rings = MotifFinder::findFranzblauRings(graph, 6);
  EXPECT_EQ(franzblau_rings[6], 1U);

  auto franzblau_cycles = MotifFinder::extractFranzblauCycles(graph, 6);
  ASSERT_EQ(franzblau_cycles.size(), 1U);
  EXPECT_EQ(franzblau_cycles[0].size(), 6U);
}

TEST_F(MotifFinderTests, FranzblauRejectsCompositeRingWithAlternateGeodesic) {
  // Graph: 6-cycle 0-1-2-3-4-5-0 plus node 6 connected to 0 and 2.
  // 0-1-2 is distance 2. 0-6-2 is also distance 2.
  // Under King, no shortcut shorter than 2 exists, so King accepts the 6-cycle.
  // Under Franzblau, the geodesic between 0 and 2 on the 6-cycle is NOT unique (two geodesics of
  // length 2 exist), so Franzblau strictly rejects the composite 6-cycle.
  graph = correlation::core::NeighborGraph(7);
  for (size_t i = 0; i < 6; ++i) {
    size_t next = (i + 1) % 6;
    graph.addDirectedEdge(i, next, 1.0, {1.0, 0.0, 0.0});
    graph.addDirectedEdge(next, i, 1.0, {-1.0, 0.0, 0.0});
  }
  // Alternate path between 0 and 2 via 6:
  graph.addDirectedEdge(0, 6, 1.0, {0.0, 1.0, 0.0});
  graph.addDirectedEdge(6, 0, 1.0, {0.0, -1.0, 0.0});
  graph.addDirectedEdge(2, 6, 1.0, {0.0, 1.0, 0.0});
  graph.addDirectedEdge(6, 2, 1.0, {0.0, -1.0, 0.0});

  // King rings include both 6-cycles (0-1-2-3-4-5-0 and 6-2-3-4-5-0-6):
  auto king_rings = MotifFinder::findRings(graph, 6, RingType::King);
  EXPECT_EQ(king_rings[6], 2U);
  EXPECT_EQ(king_rings[4], 1U); // 0-1-2-6 square

  // Franzblau primitive rings reject the 6-cycle:
  auto franzblau_rings = MotifFinder::findFranzblauRings(graph, 6);
  EXPECT_EQ(franzblau_rings.count(6), 0U);
  EXPECT_EQ(franzblau_rings[4], 1U); // Only the primitive 4-cycle remains
}

TEST_F(MotifFinderTests, BuildBridgedGraphContractsOxygenBridges) {
  // Create 6-atom cell: Si0 - O1 - Si2 - O3 - Si4 - O5 - Si0
  correlation::core::Cell silica_ring(std::array<real_t, 6>{20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  silica_ring.addAtom("Si", {0.0, 0.0, 0.0}); // 0
  silica_ring.addAtom("O", {1.0, 0.0, 0.0});  // 1
  silica_ring.addAtom("Si", {2.0, 0.0, 0.0}); // 2
  silica_ring.addAtom("O", {2.0, 1.0, 0.0});  // 3
  silica_ring.addAtom("Si", {1.0, 2.0, 0.0}); // 4
  silica_ring.addAtom("O", {0.0, 1.0, 0.0});  // 5

  graph = correlation::core::NeighborGraph(6);
  // Add bonds along the 6-cycle
  for (size_t i = 0; i < 6; ++i) {
    size_t next = (i + 1) % 6;
    graph.addDirectedEdge(i, next, 1.0, {1.0, 0.0, 0.0});
    graph.addDirectedEdge(next, i, 1.0, {-1.0, 0.0, 0.0});
  }

  auto bridged = MotifFinder::buildBridgedGraph(graph, silica_ring, "Si", "O");
  EXPECT_EQ(bridged.nodeCount(), 6U);

  // In the bridged graph, Si0 (node 0) should be connected to Si2 (node 2) and Si4 (node 4)
  EXPECT_EQ(bridged.getNeighbors(0).size(), 2U);
  EXPECT_EQ(bridged.getNeighbors(2).size(), 2U);
  EXPECT_EQ(bridged.getNeighbors(4).size(), 2U);
  // Oxygen atoms (nodes 1, 3, 5) should have 0 neighbors in the former-projected graph
  EXPECT_EQ(bridged.getNeighbors(1).size(), 0U);
  EXPECT_EQ(bridged.getNeighbors(3).size(), 0U);
  EXPECT_EQ(bridged.getNeighbors(5).size(), 0U);

  // Cycle search on bridged graph finds exactly 1 ring of size 3 (3 Si polyhedra)
  auto rings = MotifFinder::findRings(bridged, 6);
  EXPECT_EQ(rings[3], 1U);
}

TEST_F(MotifFinderTests, FilterAlternatingCyclesValidatesSequence) {
  correlation::core::Cell test_cell(std::array<real_t, 6>{20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  test_cell.addAtom("Si", {0.0, 0.0, 0.0}); // 0
  test_cell.addAtom("O", {1.0, 0.0, 0.0});  // 1
  test_cell.addAtom("Si", {2.0, 0.0, 0.0}); // 2
  test_cell.addAtom("O", {3.0, 0.0, 0.0});  // 3
  test_cell.addAtom("Si", {4.0, 0.0, 0.0}); // 4
  test_cell.addAtom("O", {5.0, 0.0, 0.0});  // 5

  std::vector<std::vector<correlation::core::AtomID>> cycles = {
      {0, 1, 2, 3, 4, 5}, // strictly alternating: Si-O-Si-O-Si-O
      {0, 2, 1, 3, 4, 5}, // non-alternating: Si-Si-O-O-Si-O
      {0, 1, 2}           // odd length (3)
  };

  auto filtered = MotifFinder::filterAlternatingCycles(cycles, test_cell, "Si", "O");
  ASSERT_EQ(filtered.size(), 1U);
  EXPECT_EQ(filtered[0], (std::vector<correlation::core::AtomID>{0, 1, 2, 3, 4, 5}));
}

} // namespace correlation::calculators
