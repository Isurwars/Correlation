// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/MotifFinder.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include <gtest/gtest.h>

namespace correlation::calculators {
namespace {
class MotifFinderTests : public ::testing::Test {
public:
  correlation::core::Cell cell;
  correlation::core::NeighborGraph graph;

  void SetUp() override { cell = correlation::core::Cell({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0}); }
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

} // namespace correlation::calculators
