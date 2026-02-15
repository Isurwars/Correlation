// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/Cell.hpp"
#include "../include/Trajectory.hpp"
#include "../include/TrajectoryAnalyzer.hpp"
#include <gtest/gtest.h>

TEST(Test05_TrajectoryAnalyzer, BasicUsage) {
  // Create a dummy cell
  Cell cell;
  cell.setLatticeParameters({10, 10, 10, 90, 90, 90});
  cell.addAtom("H", {0.5, 0.5, 0.5});
  cell.addAtom("H", {1.5, 0.5, 0.5}); // Distance 1.0

  // Create a trajectory with multiple frames
  Trajectory trajectory;
  trajectory.addFrame(cell);
  trajectory.addFrame(cell);

  double neighbor_cutoff = 2.0;
  std::vector<std::vector<double>> bond_cutoffs = {
      {1.1}}; // Assuming H-H is index 0-0, simplified

  TrajectoryAnalyzer analyzer(trajectory, neighbor_cutoff, bond_cutoffs);

  EXPECT_EQ(analyzer.getAnalyzers().size(), 2);
  EXPECT_EQ(analyzer.getNeighborCutoff(), neighbor_cutoff);
  EXPECT_EQ(analyzer.getBondCutoffsSQ()[0][0], 1.1);

  // Verify that StructureAnalyzers were created and ran
  const auto &frame_analyzers = analyzer.getAnalyzers();
  for (const auto &sa : frame_analyzers) {
    // Just checking if we can access them and they are not null
    EXPECT_TRUE(sa != nullptr);
    // In a real test we might check neighbor counts if we knew indices
  }
}
