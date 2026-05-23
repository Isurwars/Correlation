// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
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
  AngleTensor out_angles(
      num_elements,
      std::vector<std::vector<std::vector<double>>>(
          num_elements, std::vector<std::vector<double>>(num_elements)));

  // Act
  AngleCalculator::compute(cell, graph, out_angles);

  // Assert
  // Angle formed by H1(type 1) - O(type 0) - H2(type 1) should be pi / 2 (90 degrees)
  // Indices: [type1][type_central][type2] -> [1][0][1]
  ASSERT_EQ(out_angles[1][0][1].size(), 1);
  EXPECT_NEAR(out_angles[1][0][1][0], correlation::math::pi / 2.0, 1e-6);
}

} // namespace correlation::testing
