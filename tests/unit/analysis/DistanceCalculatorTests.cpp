// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/DistanceCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {

using namespace correlation::calculators;
using namespace correlation::core;

TEST(DistanceCalculatorTests, ComputesPairwiseDistancesAndNeighborGraph) {
  // Construct a cubic cell of size 10.0
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});

  // Add 2 atoms: Si (0,0,0) and O (1.5, 0.0, 0.0)
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("O", {1.5, 0.0, 0.0});

  // Setup inputs
  double cutoff_sq = 4.0; // Cutoff distance = 2.0
  std::vector<std::vector<double>> bond_cutoffs_sq = {
      {4.0, 4.0}, // Si-Si, Si-O
      {4.0, 4.0}  // O-Si, O-O
  };

  size_t num_elements = cell.elements().size();
  DistanceTensor out_distances(num_elements, std::vector<std::vector<double>>(num_elements));
  NeighborGraph out_graph(2);

  // Act
  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs_sq, true, out_distances, out_graph);

  // Assert: Check distances
  // Si (index 0, type 0), O (index 1, type 1)
  ASSERT_EQ(out_distances[0][1].size(), 1);
  EXPECT_DOUBLE_EQ(out_distances[0][1][0], 1.5);

  // Assert: Check neighbor graph connections
  EXPECT_TRUE(out_graph.areConnected(0, 1));
  EXPECT_TRUE(out_graph.areConnected(1, 0));

  const auto& neighbors = out_graph.getNeighbors(0);
  ASSERT_EQ(neighbors.size(), 1);
  EXPECT_EQ(neighbors[0].index, 1);
  EXPECT_DOUBLE_EQ(neighbors[0].distance, 1.5);
}

} // namespace correlation::testing
