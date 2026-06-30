// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "calculators/DihedralCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"
#include "math/Constants.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace {
class DihedralCalculatorTests : public ::testing::Test {
public:
  correlation::core::Cell cell;
  correlation::core::NeighborGraph graph;

protected:
  void SetUp() override {
    // Setup a non-periodic cell large enough
    cell = correlation::core::Cell({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0});
  }
};
} // namespace

TEST_F(DihedralCalculatorTests, ComputesCorrect90DegreeDihedral) {
  // 4 atoms forming a 90-degree dihedral:
  cell.addAtom("C", {1.0, 0.0, 0.0}); // A
  cell.addAtom("C", {0.0, 0.0, 0.0}); // B
  cell.addAtom("C", {0.0, 1.0, 0.0}); // C
  cell.addAtom("C", {0.0, 1.0, 1.0}); // D

  graph = correlation::core::NeighborGraph(4);
  // Bind A-B (dist 1.0)
  graph.addDirectedEdge(0, 1, 1.0, cell.atoms()[1].position() - cell.atoms()[0].position());
  graph.addDirectedEdge(1, 0, 1.0, cell.atoms()[0].position() - cell.atoms()[1].position());

  // Bind B-C (dist 1.0)
  graph.addDirectedEdge(1, 2, 1.0, cell.atoms()[2].position() - cell.atoms()[1].position());
  graph.addDirectedEdge(2, 1, 1.0, cell.atoms()[1].position() - cell.atoms()[2].position());

  // Bind C-D (dist 1.0)
  graph.addDirectedEdge(2, 3, 1.0, cell.atoms()[3].position() - cell.atoms()[2].position());
  graph.addDirectedEdge(3, 2, 1.0, cell.atoms()[2].position() - cell.atoms()[3].position());

  correlation::analysis::StructureAnalyzer::DihedralTensor dict;
  dict.resize(1, std::vector<std::vector<std::vector<std::vector<double>>>>(
                     1, std::vector<std::vector<std::vector<double>>>(
                            1, std::vector<std::vector<double>>(1, std::vector<double>()))));

  correlation::calculators::DihedralCalculator::compute(cell, graph, dict);

  const auto &angles = dict[0][0][0][0];
  ASSERT_EQ(angles.size(), 1);

  // Test the angle: expected pi/2.
  EXPECT_NEAR(angles[0], correlation::math::pi / 2.0, 1e-12);
}

TEST_F(DihedralCalculatorTests, ComputesCorrect0DegreeDihedral) {
  // A: (1, 1, 0)
  // B: (0, 1, 0)
  // C: (0, -1, 0)
  // D: (1, -1, 0)

  cell.addAtom("C", {1.0, 1.0, 0.0});  // A
  cell.addAtom("C", {0.0, 1.0, 0.0});  // B
  cell.addAtom("C", {0.0, -1.0, 0.0}); // C
  cell.addAtom("C", {1.0, -1.0, 0.0}); // D

  graph = correlation::core::NeighborGraph(4);
  graph.addDirectedEdge(0, 1, 1.0, cell.atoms()[1].position() - cell.atoms()[0].position());
  graph.addDirectedEdge(1, 0, 1.0, cell.atoms()[0].position() - cell.atoms()[1].position());
  graph.addDirectedEdge(1, 2, 2.0, cell.atoms()[2].position() - cell.atoms()[1].position());
  graph.addDirectedEdge(2, 1, 2.0, cell.atoms()[1].position() - cell.atoms()[2].position());
  graph.addDirectedEdge(2, 3, 1.0, cell.atoms()[3].position() - cell.atoms()[2].position());
  graph.addDirectedEdge(3, 2, 1.0, cell.atoms()[2].position() - cell.atoms()[3].position());

  correlation::analysis::StructureAnalyzer::DihedralTensor dict;
  dict.resize(1, std::vector<std::vector<std::vector<std::vector<double>>>>(
                     1, std::vector<std::vector<std::vector<double>>>(
                            1, std::vector<std::vector<double>>(1, std::vector<double>()))));

  correlation::calculators::DihedralCalculator::compute(cell, graph, dict);

  const auto &angles = dict[0][0][0][0];
  ASSERT_EQ(angles.size(), 1);
  EXPECT_NEAR(angles[0], 0.0, 1e-12);
}

TEST_F(DihedralCalculatorTests, ComputesCorrect180DegreeDihedral) {
  // A: (1, 1, 0)
  // B: (0, 1, 0)
  // C: (0, -1, 0)
  // D: (-1, -1, 0)

  cell.addAtom("C", {1.0, 1.0, 0.0});   // A
  cell.addAtom("C", {0.0, 1.0, 0.0});   // B
  cell.addAtom("C", {0.0, -1.0, 0.0});  // C
  cell.addAtom("C", {-1.0, -1.0, 0.0}); // D

  graph = correlation::core::NeighborGraph(4);
  graph.addDirectedEdge(0, 1, 1.0, cell.atoms()[1].position() - cell.atoms()[0].position());
  graph.addDirectedEdge(1, 0, 1.0, cell.atoms()[0].position() - cell.atoms()[1].position());
  graph.addDirectedEdge(1, 2, 2.0, cell.atoms()[2].position() - cell.atoms()[1].position());
  graph.addDirectedEdge(2, 1, 2.0, cell.atoms()[1].position() - cell.atoms()[2].position());
  graph.addDirectedEdge(2, 3, 1.0, cell.atoms()[3].position() - cell.atoms()[2].position());
  graph.addDirectedEdge(3, 2, 1.0, cell.atoms()[2].position() - cell.atoms()[3].position());

  correlation::analysis::StructureAnalyzer::DihedralTensor dict;
  dict.resize(1, std::vector<std::vector<std::vector<std::vector<double>>>>(
                     1, std::vector<std::vector<std::vector<double>>>(
                            1, std::vector<std::vector<double>>(1, std::vector<double>()))));

  correlation::calculators::DihedralCalculator::compute(cell, graph, dict);

  const auto &angles = dict[0][0][0][0];
  ASSERT_EQ(angles.size(), 1);
  EXPECT_NEAR(std::abs(angles[0]), correlation::math::pi, 1e-12);
}

TEST_F(DihedralCalculatorTests, HandlesCoincidentCentralBondSafely) {
  cell.addAtom("C", {1.0, 0.0, 0.0}); // A
  cell.addAtom("C", {0.0, 0.0, 0.0}); // B
  cell.addAtom("C", {0.0, 0.0, 0.0}); // C
  cell.addAtom("C", {0.0, 0.0, 1.0}); // D

  graph = correlation::core::NeighborGraph(4);
  graph.addDirectedEdge(0, 1, 1.0, {-1.0, 0.0, 0.0});
  graph.addDirectedEdge(1, 0, 1.0, {1.0, 0.0, 0.0});
  graph.addDirectedEdge(1, 2, 0.0, {0.0, 0.0, 0.0});
  graph.addDirectedEdge(2, 1, 0.0, {0.0, 0.0, 0.0});
  graph.addDirectedEdge(2, 3, 1.0, {0.0, 0.0, 1.0});
  graph.addDirectedEdge(3, 2, 1.0, {0.0, 0.0, -1.0});

  correlation::analysis::StructureAnalyzer::DihedralTensor dict;
  dict.resize(1, std::vector<std::vector<std::vector<std::vector<double>>>>(
                     1, std::vector<std::vector<std::vector<double>>>(
                            1, std::vector<std::vector<double>>(1, std::vector<double>()))));

  correlation::calculators::DihedralCalculator::compute(cell, graph, dict);

  const auto &angles = dict[0][0][0][0];
  EXPECT_TRUE(angles.empty());
}
