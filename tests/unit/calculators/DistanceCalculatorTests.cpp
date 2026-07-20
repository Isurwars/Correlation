// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
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
  real_t const cutoff_sq = 4.0; // Cutoff distance = 2.0
  std::vector<std::vector<real_t>> const bond_cutoffs_sq = {
      {4.0, 4.0}, // Si-Si, Si-O
      {4.0, 4.0}  // O-Si, O-O
  };

  size_t const num_elements = cell.elements().size();
  DistanceTensor out_distances(num_elements, std::vector<std::vector<real_t>>(num_elements));
  NeighborGraph out_graph(2);

  // Act
  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs_sq, true, out_distances, out_graph);

  // Assert: Check distances
  // Si (index 0, type 0), O (index 1, type 1)
  ASSERT_EQ(out_distances[0][1].size(), 1);
  EXPECT_NEAR(out_distances[0][1][0], 1.5, correlation::is_single_precision ? 1e-5 : 1e-9);

  // Assert: Check neighbor graph connections
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{1}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{1}, AtomIndex{0}));

  const auto &neighbors = out_graph.getNeighbors(0);
  ASSERT_EQ(neighbors.size(), 1);
  EXPECT_EQ(neighbors[0].index, 1);
  EXPECT_NEAR(neighbors[0].distance, 1.5, correlation::is_single_precision ? 1e-5 : 1e-9);
}

// --- Extreme / Edge-Case Tests ---

TEST(DistanceCalculatorTests, DistanceAcrossPeriodicBoundary) {
  // Atom near box edge: distance should be computed across PBC
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.5, 5.0, 5.0}); // Near left edge
  cell.addAtom("Si", {9.5, 5.0, 5.0}); // Near right edge
  // PBC distance = 1.0 (not 9.0)

  real_t const cutoff_sq = 4.0; // cutoff = 2.0
  std::vector<std::vector<real_t>> const bond_cutoffs_sq = {{4.0}};

  size_t const num_elements = cell.elements().size();
  DistanceTensor out_distances(num_elements, std::vector<std::vector<real_t>>(num_elements));
  NeighborGraph out_graph(2);

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs_sq, true, out_distances, out_graph);

  // Should find one pair at distance 1.0
  ASSERT_GE(out_distances[0][0].size(), 1);
  EXPECT_NEAR(out_distances[0][0][0], 1.0, correlation::is_single_precision ? 1e-6 : 1e-9);

  // Neighbor graph should reflect the bond
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{1}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{1}, AtomIndex{0}));
}

TEST(DistanceCalculatorTests, SingleAtomProducesNoDistances) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});

  real_t const cutoff_sq = 25.0;
  std::vector<std::vector<real_t>> const bond_cutoffs_sq = {{25.0}};

  size_t const num_elements = cell.elements().size();
  DistanceTensor out_distances(num_elements, std::vector<std::vector<real_t>>(num_elements));
  NeighborGraph out_graph(1);

  // With ignore_periodic_self_interactions = true, a single atom has no pairs
  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs_sq, true, out_distances, out_graph);

  EXPECT_TRUE(out_distances[0][0].empty());
  EXPECT_TRUE(out_graph.getNeighbors(0).empty());
}

TEST(DistanceCalculatorTests, NonOrthogonalCell) {
  // Triclinic cell
  Cell cell({5.0, 5.0, 5.0, 60.0, 60.0, 60.0});
  cell.addAtom("Ar", {0.0, 0.0, 0.0});
  cell.addAtom("Ar", {2.5, 0.0, 0.0});

  real_t const cutoff_sq = 9.0; // cutoff = 3.0
  std::vector<std::vector<real_t>> const bond_cutoffs_sq = {{9.0}};

  size_t const num_elements = cell.elements().size();
  DistanceTensor out_distances(num_elements, std::vector<std::vector<real_t>>(num_elements));
  NeighborGraph out_graph(2);

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs_sq, true, out_distances, out_graph);

  // Should find at least one pair
  EXPECT_GE(out_distances[0][0].size(), 1);
  // The computed distance should be the actual shortest distance under PBC
  for (real_t const distance : out_distances[0][0]) {
    EXPECT_GT(distance, 0.0);
    EXPECT_LE(distance, 3.0); // Within cutoff
  }
}

TEST(DistanceCalculatorTests, AtomsOutsideCutoff) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("O", {4.0, 0.0, 0.0}); // Distance = 4.0

  // Cutoff = 2.0, so this pair should NOT be found
  real_t const cutoff_sq = 4.0;
  std::vector<std::vector<real_t>> const bond_cutoffs_sq = {{4.0, 4.0}, {4.0, 4.0}};

  size_t const num_elements = cell.elements().size();
  DistanceTensor out_distances(num_elements, std::vector<std::vector<real_t>>(num_elements));
  NeighborGraph out_graph(2);

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs_sq, true, out_distances, out_graph);

  // Distance = 4.0 >= cutoff_sq = 4.0, so no pair found
  EXPECT_TRUE(out_distances[0][1].empty());
  EXPECT_TRUE(out_distances[1][0].empty());
  EXPECT_FALSE(out_graph.areConnected(AtomIndex{0}, AtomIndex{1}));
}

TEST(DistanceCalculatorTests, ThrowsOnInvalidCutoff) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  std::vector<std::vector<real_t>> const bond_cutoffs_sq = {{4.0}};
  size_t const num_elements = cell.elements().size();
  DistanceTensor out_distances(num_elements, std::vector<std::vector<real_t>>(num_elements));
  NeighborGraph out_graph(1);

  EXPECT_THROW(DistanceCalculator::compute(cell, -1.0, bond_cutoffs_sq, true, out_distances, out_graph),
               std::invalid_argument);
  EXPECT_THROW(DistanceCalculator::compute(cell, 0.0, bond_cutoffs_sq, true, out_distances, out_graph),
               std::invalid_argument);
}

} // namespace correlation::testing
