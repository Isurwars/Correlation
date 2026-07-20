// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/AngleCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"
#include "math/Constants.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {

using namespace correlation::calculators;
using namespace correlation::core;
using namespace correlation::math;

TEST(AngleCalculatorTests, ComputesCorrect90DegreeAngle) {
  // Construct a cubic cell with lattice vectors of length 10.0
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});

  // Add 3 atoms: O (0,0,0), H1 (1,0,0), H2 (0,1,0)
  cell.addAtom("O", {0.0, 0.0, 0.0});
  cell.addAtom("H", {1.0, 0.0, 0.0});
  cell.addAtom("H", {0.0, 1.0, 0.0});

  // O is index 0 (type 0)
  // H1 is index 1 (type 1)
  // H2 is index 2 (type 1)

  // Construct neighbor graph manually
  NeighborGraph graph(3);
  // Add connections: O -> H1 and O -> H2
  graph.addDirectedEdge(0, 1, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(0, 2, 1.0, {0.0, 1.0, 0.0});

  // Initialize AngleTensor
  size_t num_elements = cell.elements().size();
  AngleTensor out_angles(num_elements, std::vector<std::vector<std::vector<real_t>>>(
                                           num_elements, std::vector<std::vector<real_t>>(num_elements)));

  // Act
  AngleCalculator::compute(cell, graph, out_angles);

  // Assert
  // Angle formed by H1(type 1) - O(type 0) - H2(type 1) should be pi / 2 (90 degrees)
  // Indices: [type1][type_central][type2] -> [1][0][1]
  ASSERT_EQ(out_angles[1][0][1].size(), 1);
  EXPECT_NEAR(out_angles[1][0][1][0], correlation::math::pi / 2.0, 1e-6);
}

// --- Extreme / Edge-Case Tests ---

TEST(AngleCalculatorTests, ComputesCorrect180DegreeAngle) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  // Linear arrangement: H1 - O - H2 along x-axis (180 degrees)
  cell.addAtom("O", {5.0, 5.0, 5.0});
  cell.addAtom("H", {6.0, 5.0, 5.0});
  cell.addAtom("H", {4.0, 5.0, 5.0});

  NeighborGraph graph(3);
  graph.addDirectedEdge(0, 1, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(0, 2, 1.0, {-1.0, 0.0, 0.0});

  size_t num_elements = cell.elements().size();
  AngleTensor out_angles(num_elements, std::vector<std::vector<std::vector<real_t>>>(
                                           num_elements, std::vector<std::vector<real_t>>(num_elements)));

  AngleCalculator::compute(cell, graph, out_angles);

  ASSERT_EQ(out_angles[1][0][1].size(), 1);
  EXPECT_NEAR(out_angles[1][0][1][0], correlation::math::pi, 1e-6);
}

TEST(AngleCalculatorTests, ComputesCorrect60DegreeAngle) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  // Equilateral triangle arrangement: 60 degrees
  cell.addAtom("C", {5.0, 5.0, 5.0});                        // Center
  cell.addAtom("C", {6.0, 5.0, 5.0});                        // +x
  cell.addAtom("C", {5.5, 5.0 + std::sqrt(3.0) / 2.0, 5.0}); // 60° offset

  NeighborGraph graph(3);
  graph.addDirectedEdge(0, 1, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(0, 2, 1.0, {0.5, std::sqrt(3.0) / 2.0, 0.0});

  size_t num_elements = cell.elements().size();
  AngleTensor out_angles(num_elements, std::vector<std::vector<std::vector<real_t>>>(
                                           num_elements, std::vector<std::vector<real_t>>(num_elements)));

  AngleCalculator::compute(cell, graph, out_angles);

  // C-C-C angle should be 60 degrees = pi/3
  ASSERT_EQ(out_angles[0][0][0].size(), 1);
  EXPECT_NEAR(out_angles[0][0][0][0], correlation::math::pi / 3.0, 1e-6);
}

TEST(AngleCalculatorTests, SingleNeighborProducesNoAngles) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("O", {5.0, 5.0, 5.0});
  cell.addAtom("H", {6.0, 5.0, 5.0});

  // O has only 1 neighbor — no angle can be formed
  NeighborGraph graph(2);
  graph.addDirectedEdge(0, 1, 1.0, {1.0, 0.0, 0.0});

  size_t num_elements = cell.elements().size();
  AngleTensor out_angles(num_elements, std::vector<std::vector<std::vector<real_t>>>(
                                           num_elements, std::vector<std::vector<real_t>>(num_elements)));

  AngleCalculator::compute(cell, graph, out_angles);

  // No angles should be produced (need at least 2 neighbors on one atom)
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
      for (size_t k = 0; k < num_elements; ++k) {
        EXPECT_TRUE(out_angles[i][j][k].empty());
      }
    }
  }
}

TEST(AngleCalculatorTests, OverlappingAtomsProduceNoAngles) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("O", {0.0, 0.0, 0.0});
  cell.addAtom("H", {0.0, 0.0, 0.0});
  cell.addAtom("H", {0.0, 0.0, 0.0});

  NeighborGraph graph(3);
  graph.addDirectedEdge(0, 1, 0.0, {0.0, 0.0, 0.0});
  graph.addDirectedEdge(0, 2, 0.0, {0.0, 0.0, 0.0});

  size_t num_elements = cell.elements().size();
  AngleTensor out_angles(num_elements, std::vector<std::vector<std::vector<real_t>>>(
                                           num_elements, std::vector<std::vector<real_t>>(num_elements)));

  AngleCalculator::compute(cell, graph, out_angles);

  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
      for (size_t k = 0; k < num_elements; ++k) {
        EXPECT_TRUE(out_angles[i][j][k].empty());
      }
    }
  }
}

} // namespace correlation::testing
