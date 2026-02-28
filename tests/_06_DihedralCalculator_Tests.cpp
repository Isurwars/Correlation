// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>
#include <vector>

#include "../include/Cell.hpp"
#include "../include/NeighborGraph.hpp"
#include "../include/calculators/DihedralCalculator.hpp"

class _06_DihedralCalculator_Tests : public ::testing::Test {
protected:
  Cell cell;
  NeighborGraph graph;

  void SetUp() override {
    // Setup a non-periodic cell large enough
    cell = Cell({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0});
  }
};

TEST_F(_06_DihedralCalculator_Tests, ComputesCorrect90DegreeDihedral) {
  // 4 atoms forming a 90-degree dihedral:
  // A: (1, 0, 0)
  // B: (0, 0, 0)
  // C: (0, 1, 0)
  // D: (0, 1, 1)

  cell.addAtom("C", {1.0, 0.0, 0.0}); // A
  cell.addAtom("C", {0.0, 0.0, 0.0}); // B
  cell.addAtom("C", {0.0, 1.0, 0.0}); // C
  cell.addAtom("C", {0.0, 1.0, 1.0}); // D

  graph = NeighborGraph(4);
  // Bind A-B (dist 1.0)
  graph.addDirectedEdge(
      0, 1, 1.0, cell.atoms()[1].position() - cell.atoms()[0].position());
  graph.addDirectedEdge(
      1, 0, 1.0, cell.atoms()[0].position() - cell.atoms()[1].position());

  // Bind B-C (dist 1.0)
  graph.addDirectedEdge(
      1, 2, 1.0, cell.atoms()[2].position() - cell.atoms()[1].position());
  graph.addDirectedEdge(
      2, 1, 1.0, cell.atoms()[1].position() - cell.atoms()[2].position());

  // Bind C-D (dist 1.0)
  graph.addDirectedEdge(
      2, 3, 1.0, cell.atoms()[3].position() - cell.atoms()[2].position());
  graph.addDirectedEdge(
      3, 2, 1.0, cell.atoms()[2].position() - cell.atoms()[3].position());

  calculators::DihedralTensor dict;
  dict.resize(1, std::vector<std::vector<std::vector<std::vector<double>>>>(
                     1, std::vector<std::vector<std::vector<double>>>(
                            1, std::vector<std::vector<double>>(
                                   1, std::vector<double>()))));

  calculators::DihedralCalculator::compute(cell, graph, dict);

  const auto &angles = dict[0][0][0][0];
  // Expected: either 1 angle if we only do A-B-C-D, or 2 if we do A-B-C-D and
  // D-C-B-A. The implementation only stores A-B-C-D once for the same element
  // types if we didn't add the reverse condition. Wait, the implementation
  // handles reverse based on element types, but for same elements, A-B-C-D is
  // pushed, and D-C-B-A is pushed if types differ. Since types are same, it
  // pushes once per sequence. However, the loop will traverse B=1, C=2 (which
  // gives A=0, D=3 -> A-B-C-D). And it will NOT traverse B=2, C=1 because B < C
  // restriction. So there should be exactly 1 angle.

  ASSERT_EQ(angles.size(), 1);

  // Test the angle: expected pi/2.
  EXPECT_NEAR(angles[0], M_PI / 2.0, 1e-5);
}

TEST_F(_06_DihedralCalculator_Tests, ComputesCorrect0DegreeDihedral) {
  // A: (1, 1, 0)
  // B: (0, 1, 0)
  // C: (0, -1, 0)
  // D: (1, -1, 0)

  cell.addAtom("C", {1.0, 1.0, 0.0});  // A
  cell.addAtom("C", {0.0, 1.0, 0.0});  // B
  cell.addAtom("C", {0.0, -1.0, 0.0}); // C
  cell.addAtom("C", {1.0, -1.0, 0.0}); // D

  graph = NeighborGraph(4);
  graph.addDirectedEdge(
      0, 1, 1.0, cell.atoms()[1].position() - cell.atoms()[0].position());
  graph.addDirectedEdge(
      1, 0, 1.0, cell.atoms()[0].position() - cell.atoms()[1].position());
  graph.addDirectedEdge(
      1, 2, 2.0, cell.atoms()[2].position() - cell.atoms()[1].position());
  graph.addDirectedEdge(
      2, 1, 2.0, cell.atoms()[1].position() - cell.atoms()[2].position());
  graph.addDirectedEdge(
      2, 3, 1.0, cell.atoms()[3].position() - cell.atoms()[2].position());
  graph.addDirectedEdge(
      3, 2, 1.0, cell.atoms()[2].position() - cell.atoms()[3].position());

  calculators::DihedralTensor dict;
  dict.resize(1, std::vector<std::vector<std::vector<std::vector<double>>>>(
                     1, std::vector<std::vector<std::vector<double>>>(
                            1, std::vector<std::vector<double>>(
                                   1, std::vector<double>()))));

  calculators::DihedralCalculator::compute(cell, graph, dict);

  const auto &angles = dict[0][0][0][0];
  ASSERT_EQ(angles.size(), 1);
  EXPECT_NEAR(angles[0], 0.0, 1e-5);
}

TEST_F(_06_DihedralCalculator_Tests, ComputesCorrect180DegreeDihedral) {
  // A: (1, 1, 0)
  // B: (0, 1, 0)
  // C: (0, -1, 0)
  // D: (-1, -1, 0)

  cell.addAtom("C", {1.0, 1.0, 0.0});   // A
  cell.addAtom("C", {0.0, 1.0, 0.0});   // B
  cell.addAtom("C", {0.0, -1.0, 0.0});  // C
  cell.addAtom("C", {-1.0, -1.0, 0.0}); // D

  graph = NeighborGraph(4);
  graph.addDirectedEdge(
      0, 1, 1.0, cell.atoms()[1].position() - cell.atoms()[0].position());
  graph.addDirectedEdge(
      1, 0, 1.0, cell.atoms()[0].position() - cell.atoms()[1].position());
  graph.addDirectedEdge(
      1, 2, 2.0, cell.atoms()[2].position() - cell.atoms()[1].position());
  graph.addDirectedEdge(
      2, 1, 2.0, cell.atoms()[1].position() - cell.atoms()[2].position());
  graph.addDirectedEdge(
      2, 3, 1.0, cell.atoms()[3].position() - cell.atoms()[2].position());
  graph.addDirectedEdge(
      3, 2, 1.0, cell.atoms()[2].position() - cell.atoms()[3].position());

  calculators::DihedralTensor dict;
  dict.resize(1, std::vector<std::vector<std::vector<std::vector<double>>>>(
                     1, std::vector<std::vector<std::vector<double>>>(
                            1, std::vector<std::vector<double>>(
                                   1, std::vector<double>()))));

  calculators::DihedralCalculator::compute(cell, graph, dict);

  const auto &angles = dict[0][0][0][0];
  ASSERT_EQ(angles.size(), 1);
  EXPECT_NEAR(std::abs(angles[0]), M_PI, 1e-5);
}
