/**
 * @file BaselineRegressionTests.cpp
 * @brief Golden baseline regression test suite for distance, graph and ring calculations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "../../CrystalTestHelper.hpp"
#include "calculators/spatial/DistanceCalculator.hpp"
#include "calculators/spatial/MotifFinder.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"
#include "math/Precision.hpp"
#include "mlip/GraphDescriptors.hpp"
#include "mlip/PeriodicGraphBuilder.hpp"

#include <gtest/gtest.h>

namespace correlation::calculators {
namespace {

correlation::analysis::BondCutoffMatrix makeCutoffMatrix(real_t min_sq, real_t max_sq) {
  correlation::analysis::BondCutoffMatrix matrix(1);
  matrix[0][0] = correlation::analysis::BondCutoffRange{
      .min_sq = min_sq,
      .max_sq = max_sq,
  };
  return matrix;
}

core::NeighborGraph buildDiamondTestGraph(int n_cells) {
  const auto cell = testing::crystals::createDiamondCell(5.43, "Si", n_cells, n_cells, n_cells);
  const real_t cutoff = 3.0;
  const real_t cutoff_sq = cutoff * cutoff;

  auto bond_cutoffs = makeCutoffMatrix(0.0, 2.8 * 2.8);
  core::NeighborGraph graph(cell.atomCount());
  calculators::DistanceCalculationConfig config{
      .r_max = cutoff,
      .r_bin_width = 0.0,
      .num_bins = 0,
  };

  calculators::DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, graph, nullptr,
                                           config);
  return graph;
}

TEST(BaselineRegressionTests, SmallCellDistanceBaseline) {
  const auto cell = testing::crystals::createSimpleCubicCell(4.0, "Ar", 3, 3, 3);
  const real_t cutoff = 6.0;
  const real_t cutoff_sq = cutoff * cutoff;
  const auto bond_cutoffs = makeCutoffMatrix(0.0, 4.5 * 4.5);
  const size_t num_bins = 60;
  const calculators::DistanceCalculationConfig config{
      .r_max = cutoff,
      .r_bin_width = static_cast<real_t>(0.1),
      .num_bins = num_bins,
  };

  core::NeighborGraph graph(cell.atomCount());
  calculators::RawHistogramTensor histograms;
  calculators::DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, graph, &histograms,
                                           config);

  size_t total_edges = 0;
  for (size_t i = 0; i < graph.nodeCount(); ++i) {
    total_edges += graph.getNeighbors(i).size();
  }

  real_t total_hist_sum = 0.0;
  for (const auto &row : histograms) {
    for (const auto &col : row) {
      for (real_t count : col) {
        total_hist_sum += count;
      }
    }
  }

  // Golden baseline assertions
  EXPECT_EQ(graph.nodeCount(), 27);
  EXPECT_EQ(total_edges, 162); // 6 neighbors per atom in SC lattice (27 * 6 = 162)
  EXPECT_DOUBLE_EQ(total_hist_sum, 243.0);
}

TEST(BaselineRegressionTests, MidCellDistanceBaseline) {
  const auto cell = testing::crystals::createFCCCell(3.61, "Cu", 4, 4, 4);
  const real_t cutoff = 5.0;
  const real_t cutoff_sq = cutoff * cutoff;
  const auto bond_cutoffs = makeCutoffMatrix(0.0, 2.9 * 2.9);
  const size_t num_bins = 100;
  const calculators::DistanceCalculationConfig config{
      .r_max = cutoff,
      .r_bin_width = static_cast<real_t>(0.05),
      .num_bins = num_bins,
  };

  core::NeighborGraph graph(cell.atomCount());
  calculators::RawHistogramTensor histograms;
  calculators::DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, graph, &histograms,
                                           config);

  size_t total_edges = 0;
  for (size_t i = 0; i < graph.nodeCount(); ++i) {
    total_edges += graph.getNeighbors(i).size();
  }

  real_t total_hist_sum = 0.0;
  for (const auto &row : histograms) {
    for (const auto &col : row) {
      for (real_t count : col) {
        total_hist_sum += count;
      }
    }
  }

  EXPECT_EQ(graph.nodeCount(), 256);
  EXPECT_EQ(total_edges, 3072); // 12 nn per atom in FCC lattice (256 * 12 = 3072)
  EXPECT_DOUBLE_EQ(total_hist_sum, 5376.0);
}

TEST(BaselineRegressionTests, DiamondMotifRingsBaseline) {
  const auto graph = buildDiamondTestGraph(2); // 8 * 2^3 = 64 atoms
  EXPECT_EQ(graph.nodeCount(), 64);

  // Each atom in diamond Si has coordination 4
  for (size_t i = 0; i < graph.nodeCount(); ++i) {
    EXPECT_EQ(graph.getNeighbors(i).size(), 4);
  }

  const auto king_rings = calculators::MotifFinder::findRings(graph, 6, analysis::RingType::King);
  const auto franz_rings =
      calculators::MotifFinder::findRings(graph, 6, analysis::RingType::Franzblau);

  // No 3, 4, or 5-membered rings in diamond
  EXPECT_EQ(king_rings.count(3), 0);
  EXPECT_EQ(king_rings.count(4), 0);
  EXPECT_EQ(king_rings.count(5), 0);
  EXPECT_EQ(king_rings.at(6), 128);

  EXPECT_EQ(franz_rings.count(3), 0);
  EXPECT_EQ(franz_rings.count(4), 0);
  EXPECT_EQ(franz_rings.count(5), 0);
  EXPECT_EQ(franz_rings.at(6), 128);
}

TEST(BaselineRegressionTests, GraphDescriptorsRingBaseline) {
  const auto cell = testing::crystals::createDiamondCell(5.43, "Si", 2, 2, 2);
  const auto graph = mlip::PeriodicGraphBuilder::buildGraph(cell, static_cast<real_t>(2.8));
  const size_t max_size = 6;

  const auto desc = mlip::GraphDescriptors::computeRingStatisticsDescriptor(graph, max_size);
  ASSERT_EQ(desc.size(), graph.atom_count * max_size);

  real_t total_ring_desc = 0.0;
  for (real_t value : desc) {
    total_ring_desc += value;
  }
  EXPECT_DOUBLE_EQ(total_ring_desc, 768.0);
}

} // namespace
} // namespace correlation::calculators
