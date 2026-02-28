// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>

#include "NeighborGraph.hpp"
#include "calculators/MDCalculator.hpp"

class MDCalculatorTest : public ::testing::Test {
protected:
  NeighborGraph graph;

  void SetUp() override {
    // Generate a simple triangle (3-ring)
    graph = NeighborGraph(3);
    graph.addDirectedEdge(0, 1, 1.0, {1, 0, 0});
    graph.addDirectedEdge(1, 0, 1.0, {-1, 0, 0});

    graph.addDirectedEdge(1, 2, 1.0, {0, 1, 0});
    graph.addDirectedEdge(2, 1, 1.0, {0, -1, 0});

    graph.addDirectedEdge(2, 0, 1.414, {-1, -1, 0});
    graph.addDirectedEdge(0, 2, 1.414, {1, 1, 0});
  }
};

TEST_F(MDCalculatorTest, ComputeMotif) {
  size_t max_ring_size = 5;
  Histogram f_motif = MDCalculator::calculate(graph, max_ring_size);

  ASSERT_FALSE(f_motif.partials.empty());
  ASSERT_EQ(f_motif.bins.size(), max_ring_size - 2);

  // Bins should be 3, 4, 5
  EXPECT_EQ(f_motif.bins[0], 3.0);
  EXPECT_EQ(f_motif.bins[1], 4.0);
  EXPECT_EQ(f_motif.bins[2], 5.0);

  auto &partial = f_motif.partials["Rings"];
  EXPECT_EQ(partial.size(), max_ring_size - 2);

  // We expect 1 triangle
  EXPECT_EQ(partial[0], 1.0);
  // No 4-rings or 5-rings
  EXPECT_EQ(partial[1], 0.0);
  EXPECT_EQ(partial[2], 0.0);
}

TEST_F(MDCalculatorTest, InvalidMaxRingSize) {
  EXPECT_THROW(MDCalculator::calculate(graph, 2), std::invalid_argument);
}
