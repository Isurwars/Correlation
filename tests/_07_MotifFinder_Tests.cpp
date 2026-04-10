// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/MotifFinder.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include <gtest/gtest.h>

namespace correlation::calculators {

class _07_MotifFinder_Tests : public ::testing::Test {
protected:
  correlation::core::Cell cell;
  correlation::core::NeighborGraph graph;

  void SetUp() override {
    cell = correlation::core::Cell({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0},
                                   {0.0, 0.0, 20.0});
  }
};

TEST_F(_07_MotifFinder_Tests, DetectsSingleTriangle) {
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

TEST_F(_07_MotifFinder_Tests, DetectsSingleSquare) {
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

} // namespace correlation::calculators
