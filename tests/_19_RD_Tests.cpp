// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/StructureAnalyzer.hpp"
#include "calculators/RDCalculator.hpp"
#include "core/NeighborGraph.hpp"
#include "core/Trajectory.hpp"
#include "readers/FileReader.hpp"

#include <fstream>
#include <gtest/gtest.h>

namespace correlation::analysis {

class _19_RD_Tests : public ::testing::Test {
protected:
  correlation::core::NeighborGraph graph;

  void SetUp() override {
    // Generate a simple triangle (3-ring)
    graph = correlation::core::NeighborGraph(3);
    graph.addDirectedEdge(0, 1, 1.0, {1, 0, 0});
    graph.addDirectedEdge(1, 0, 1.0, {-1, 0, 0});

    graph.addDirectedEdge(1, 2, 1.0, {0, 1, 0});
    graph.addDirectedEdge(2, 1, 1.0, {0, -1, 0});

    graph.addDirectedEdge(2, 0, 1.414, {-1, -1, 0});
    graph.addDirectedEdge(0, 2, 1.414, {1, 1, 0});
  }
};

TEST_F(_19_RD_Tests, ComputeMotif) {
  size_t max_ring_size = 5;
  Histogram f_motif =
      correlation::calculators::RDCalculator::calculate(graph, max_ring_size);

  EXPECT_EQ(f_motif.x_label, "Ring Size");
  EXPECT_EQ(f_motif.bins.size(), 3);

  // Check bins are up to 5
  EXPECT_EQ(f_motif.bins[0], 3.0);
  EXPECT_EQ(f_motif.bins[1], 4.0);
  EXPECT_EQ(f_motif.bins[2], 5.0);

  auto &partial = f_motif.partials.at("Rings");
  EXPECT_EQ(partial.size(), 3);
  EXPECT_EQ(partial[0], 1.0);
  EXPECT_EQ(partial[1], 0.0);
  EXPECT_EQ(partial[2], 0.0);
}

TEST_F(_19_RD_Tests, InvalidMaxRingSize) {
  EXPECT_THROW(correlation::calculators::RDCalculator::calculate(graph, 2),
               std::invalid_argument);
}

TEST_F(_19_RD_Tests, CelluloseRingDistribution) {
  std::string cellulose_path = "../../examples/Cellulose/Cellulose.md";
  std::ifstream f(cellulose_path);
  if (!f.good()) {
    cellulose_path = "../examples/Cellulose/Cellulose.md";
    std::ifstream f2(cellulose_path);
    if (!f2.good()) {
      cellulose_path = "examples/Cellulose/Cellulose.md";
      std::ifstream f3(cellulose_path);
      if (!f3.good()) {
        GTEST_SKIP() << "Cellulose.md example file not found. Skipping test.";
        return;
      }
    }
  }

  correlation::readers::FileType type =
      correlation::readers::determineFileType(cellulose_path);
  correlation::core::Trajectory traj =
      correlation::readers::readTrajectory(cellulose_path, type);
  if (traj.getFrames().empty()) {
    GTEST_SKIP() << "No frames found in Cellulose.md";
    return;
  }

  traj.precomputeBondCutoffs();
  auto cutoffs = traj.getBondCutoffsSQ();

  // Custom bond cutoff: The covalent radius of O is 0.73, sum is 1.46.
  // The default bond factor is 1.2, making the cutoff 1.752.
  // We need the O-O cutoff to be smaller than the H-bond distance (1.478)
  // Let's set the O-O cutoff squared to 1.45^2 = 2.1025.

  // Find the index for Oxygen
  int o_idx = -1;
  const auto &elements = traj.getFrames()[0].elements();
  for (size_t i = 0; i < elements.size(); ++i) {
    if (elements[i].symbol == "O") {
      o_idx = i;
      break;
    }
  }

  if (o_idx != -1) {
    cutoffs[o_idx][o_idx] = 1.45 * 1.45;
  }
  traj.setBondCutoffsSQ(cutoffs);

  correlation::core::Cell frame = traj.getFrames()[0];
  StructureAnalyzer analyzer(frame, 3.0, traj.getBondCutoffsSQ(), true);
  const correlation::core::NeighborGraph &graph_cellulose =
      analyzer.neighborGraph();

  size_t max_ring_size = 10;

  Histogram f_motif = correlation::calculators::RDCalculator::calculate(
      graph_cellulose, max_ring_size);

  EXPECT_EQ(f_motif.x_label, "Ring Size");
  ASSERT_EQ(f_motif.bins.size(), max_ring_size - 2);

  auto &partial = f_motif.partials.at("Rings");
  EXPECT_EQ(partial.size(), max_ring_size - 2);

  double sum_rings = 0.0;
  for (size_t i = 0; i < partial.size(); ++i) {
    if (i == 3) { // Bin 3 corresponds to ring size 6 (since bins start at 3)
      EXPECT_NEAR(partial[i], 1.0, 1e-6)
          << "Ring size 6 should have all the counts.";
    } else {
      EXPECT_NEAR(partial[i], 0.0, 1e-6)
          << "Ring size " << (i + 3) << " should be empty.";
    }
    sum_rings += partial[i];
  }
  EXPECT_NEAR(sum_rings, 1.0, 1e-6);
}

} // namespace correlation::analysis

// -------------------------------------------------------------------------- //
// ----------------------------- Main function ------------------------------ //
