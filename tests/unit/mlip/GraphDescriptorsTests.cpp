/**
 * @file GraphDescriptorsTests.cpp
 * @brief Unit tests for topological, structural, and spectral graph descriptors.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/spatial/MotifFinder.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"
#include "mlip/GraphDescriptors.hpp"
#include "mlip/PeriodicGraphBuilder.hpp"

#include <gtest/gtest.h>

#include "../../CrystalTestHelper.hpp"

namespace correlation::mlip {

namespace {

correlation::core::Cell makeFccCell(real_t lattice_a = static_cast<real_t>(3.615)) {
  return correlation::testing::crystals::createFCCCell(lattice_a, "Cu", 2, 2, 2);
}

correlation::core::Cell makeBccCell(real_t lattice_a = static_cast<real_t>(2.866)) {
  return correlation::testing::crystals::createBCCCell(lattice_a, "Fe", 3, 3, 3);
}

/// Legacy reference: one King ring search per ring size (pre-optimization behavior).
std::vector<real_t> referenceRingDescriptor(const PeriodicGraphData &graph, size_t max_size) {
  std::vector<real_t> desc(graph.atom_count * max_size, static_cast<real_t>(0.0));
  correlation::core::NeighborGraph neighbor_graph(graph.atom_count);
  for (size_t edge = 0; edge < graph.edge_count; ++edge) {
    neighbor_graph.addDirectedEdge(
        static_cast<correlation::core::AtomID>(graph.edge_index_flat[edge]),
        static_cast<correlation::core::AtomID>(graph.edge_index_flat[graph.edge_count + edge]),
        graph.edge_distances.size() == graph.edge_count ? graph.edge_distances[edge]
                                                        : static_cast<real_t>(0.0),
        correlation::math::Vector3<real_t>{0.0, 0.0, 0.0});
  }
  neighbor_graph.finalize();
  for (size_t size = 3; size <= max_size; ++size) {
    for (const auto &cycle :
         correlation::calculators::MotifFinder::extractCycles(neighbor_graph, size)) {
      for (const auto atom_id : cycle) {
        desc[static_cast<size_t>(atom_id) * max_size + (size - 1)] += static_cast<real_t>(1.0);
      }
    }
  }
  return desc;
}

void expectRingDescriptorMatchesReference(const correlation::core::Cell &cell, real_t cutoff,
                                          size_t max_size) {
  const auto graph = PeriodicGraphBuilder::buildGraph(cell, cutoff);
  const auto actual = GraphDescriptors::computeRingStatisticsDescriptor(graph, max_size);
  const auto expected = referenceRingDescriptor(graph, max_size);
  ASSERT_EQ(actual.size(), expected.size());
  for (size_t i = 0; i < actual.size(); ++i) {
    EXPECT_DOUBLE_EQ(actual[i], expected[i]) << "index " << i;
  }
}

} // anonymous namespace

TEST(GraphDescriptorsTests, EmptyGraphReturnsEmptyVectors) {
  const PeriodicGraphData empty_graph;
  EXPECT_TRUE(GraphDescriptors::computeRingStatisticsDescriptor(empty_graph).empty());
  EXPECT_TRUE(GraphDescriptors::computeCNADescriptor(empty_graph).empty());
  EXPECT_TRUE(GraphDescriptors::computeCoordinationEmbedding(empty_graph).empty());
  EXPECT_TRUE(GraphDescriptors::computeGraphSpectrum(empty_graph, 5).empty());
}

TEST(GraphDescriptorsTests, CoordinationEmbeddingMatchesGraphDegrees) {
  const auto cell = makeFccCell();
  // Cutoff covering first shell (a / sqrt(2) ~ 2.556 A)
  const auto graph = PeriodicGraphBuilder::buildGraph(cell, static_cast<real_t>(2.8));

  ASSERT_EQ(graph.atom_count, 32);
  const auto coord = GraphDescriptors::computeCoordinationEmbedding(graph);
  ASSERT_EQ(coord.size(), 32);

  // In FCC unit cell under PBC, each Cu atom has 12 nearest neighbors
  for (size_t i = 0; i < 32; ++i) {
    EXPECT_EQ(coord[i], static_cast<real_t>(12.0));
  }
}

TEST(GraphDescriptorsTests, CNAClassifiesFccMotifs) {
  const auto cell = makeFccCell();
  const auto graph = PeriodicGraphBuilder::buildGraph(cell, static_cast<real_t>(2.8));

  ASSERT_EQ(graph.atom_count, 32);
  const auto cna_labels = GraphDescriptors::computeCNADescriptor(graph);
  ASSERT_EQ(cna_labels.size(), 32);

  // All 32 atoms in ideal FCC lattice must be classified as FCC (1)
  for (size_t i = 0; i < 32; ++i) {
    EXPECT_EQ(cna_labels[i], static_cast<int>(CNALabel::FCC));
  }
}

TEST(GraphDescriptorsTests, CNAClassifiesBccMotifs) {
  const auto cell = makeBccCell();
  const auto graph = PeriodicGraphBuilder::buildGraph(cell, static_cast<real_t>(3.2));

  ASSERT_EQ(graph.atom_count, 54);
  const auto cna_labels = GraphDescriptors::computeCNADescriptor(graph);
  ASSERT_EQ(cna_labels.size(), 54);

  for (size_t i = 0; i < 54; ++i) {
    EXPECT_EQ(cna_labels[i], static_cast<int>(CNALabel::BCC));
  }
}

TEST(GraphDescriptorsTests, RingStatisticsDescriptorExtraction) {
  const auto cell = makeFccCell();
  const auto graph = PeriodicGraphBuilder::buildGraph(cell, static_cast<real_t>(2.8));

  const size_t max_ring_size = 4;
  const auto rings = GraphDescriptors::computeRingStatisticsDescriptor(graph, max_ring_size);

  ASSERT_EQ(rings.size(), graph.atom_count * max_ring_size);
  // Triangle rings (size 3, index 2) exist abundantly in FCC {111} close-packed planes
  for (size_t i = 0; i < graph.atom_count; ++i) {
    EXPECT_GT(rings[i * max_ring_size + 2], 0.0);
  }
}

TEST(GraphDescriptorsTests, RingDescriptorMatchesPerSizeReferenceFcc) {
  expectRingDescriptorMatchesReference(makeFccCell(), static_cast<real_t>(2.8), 6);
}

TEST(GraphDescriptorsTests, RingDescriptorMatchesPerSizeReferenceDiamond) {
  const auto cell = correlation::testing::crystals::createDiamondCell(5.43, "Si", 2, 2, 2);
  expectRingDescriptorMatchesReference(cell, static_cast<real_t>(2.8), 8);
}

TEST(GraphDescriptorsTests, RingDescriptorMatchesPerSizeReferenceBcc) {
  expectRingDescriptorMatchesReference(makeBccCell(), static_cast<real_t>(3.2), 5);
}

TEST(GraphDescriptorsTests, GraphSpectrumReturnsSortedEigenvalues) {
  const auto cell = makeFccCell();
  const auto graph = PeriodicGraphBuilder::buildGraph(cell, static_cast<real_t>(2.8));

  const size_t num_k = 3;
  const auto spectrum = GraphDescriptors::computeGraphSpectrum(graph, num_k);

  ASSERT_EQ(spectrum.size(), num_k);
  for (size_t i = 1; i < spectrum.size(); ++i) {
    EXPECT_GE(spectrum[i - 1], spectrum[i]);
  }
}

TEST(GraphDescriptorsTests, PopulateDescriptorsEnrichesGraphBuffers) {
  const auto cell = makeFccCell();
  auto graph = PeriodicGraphBuilder::buildGraph(cell, static_cast<real_t>(2.8));

  EXPECT_TRUE(graph.cna_labels.empty());
  EXPECT_TRUE(graph.coordination_desc.empty());
  EXPECT_TRUE(graph.ring_desc.empty());

  GraphDescriptors::populateDescriptors(graph, 4);

  EXPECT_EQ(graph.cna_labels.size(), graph.atom_count);
  EXPECT_EQ(graph.coordination_desc.size(), graph.atom_count);
  EXPECT_EQ(graph.ring_desc.size(), graph.atom_count * 4);
}

} // namespace correlation::mlip
