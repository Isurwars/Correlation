/**
 * @file GraphDescriptorsTests.cpp
 * @brief Unit tests for topological, structural, and spectral graph descriptors.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/Cell.hpp"
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
