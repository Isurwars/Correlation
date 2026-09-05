/**
 * @file PeriodicGraphBuilderTests.cpp
 * @brief Unit tests for PeriodicGraphBuilder and TorchGNNModel.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/Cell.hpp"
#include "mlip/PeriodicGraphBuilder.hpp"
#include "mlip/TorchGNNModel.hpp"

#include <gtest/gtest.h>

namespace correlation::mlip {

TEST(PeriodicGraphBuilderTests, AtomicNumberResolution) {
  EXPECT_EQ(PeriodicGraphBuilder::getAtomicNumber("H"), 1);
  EXPECT_EQ(PeriodicGraphBuilder::getAtomicNumber("He"), 2);
  EXPECT_EQ(PeriodicGraphBuilder::getAtomicNumber("C"), 6);
  EXPECT_EQ(PeriodicGraphBuilder::getAtomicNumber("N"), 7);
  EXPECT_EQ(PeriodicGraphBuilder::getAtomicNumber("O"), 8);
  EXPECT_EQ(PeriodicGraphBuilder::getAtomicNumber("Si"), 14);
  EXPECT_EQ(PeriodicGraphBuilder::getAtomicNumber("Fe"), 26);
  EXPECT_EQ(PeriodicGraphBuilder::getAtomicNumber("Au"), 79);
  EXPECT_EQ(PeriodicGraphBuilder::getAtomicNumber("Og"), 118);
  EXPECT_EQ(PeriodicGraphBuilder::getAtomicNumber("NonExistent"), 0);
}

TEST(PeriodicGraphBuilderTests, EmptyCell) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  auto graph = PeriodicGraphBuilder::buildGraph(cell, 5.0, false);

  EXPECT_EQ(graph.atom_count, 0U);
  EXPECT_EQ(graph.edge_count, 0U);
  EXPECT_TRUE(graph.positions_flat.empty());
  EXPECT_TRUE(graph.atomic_numbers.empty());
  EXPECT_TRUE(graph.edge_index_flat.empty());
}

TEST(PeriodicGraphBuilderTests, SimpleCubicFirstNeighborShell) {
  // Simple cubic unit cell: 1 atom at origin, a = 4.0 Å
  // First neighbor shell distance = 4.0 Å, coordination number = 6
  correlation::core::Cell cell({4.0, 4.0, 4.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  // Cutoff 4.1 Å encompasses the 6 periodic nearest neighbors
  auto graph = PeriodicGraphBuilder::buildGraph(cell, 4.1, false);

  EXPECT_EQ(graph.atom_count, 1U);
  EXPECT_EQ(graph.edge_count, 6U);
  ASSERT_EQ(graph.atomic_numbers.size(), 1U);
  EXPECT_EQ(graph.atomic_numbers[0], 14);

  ASSERT_EQ(graph.edge_index_flat.size(), 12U);  // 2 * 6
  ASSERT_EQ(graph.edge_shifts_flat.size(), 18U); // 6 * 3

  // All edges should originate from node 0 and point to node 0 across periodic boundaries
  for (size_t edge = 0; edge < 6; ++edge) {
    EXPECT_EQ(graph.edge_index_flat[static_cast<size_t>(0 * 6) + edge], 0);
    EXPECT_EQ(graph.edge_index_flat[static_cast<size_t>(1 * 6) + edge], 0);

    const real_t s_x = graph.edge_shifts_flat[edge * 3 + 0];
    const real_t s_y = graph.edge_shifts_flat[edge * 3 + 1];
    const real_t s_z = graph.edge_shifts_flat[edge * 3 + 2];
    const real_t shift_norm_sq = s_x * s_x + s_y * s_y + s_z * s_z;
    EXPECT_NEAR(shift_norm_sq, static_cast<real_t>(1.0), 1e-5);
  }
}

TEST(PeriodicGraphBuilderTests, DimerWithinCutoff) {
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  cell.addAtom("C", {0.0, 0.0, 0.0});
  cell.addAtom("O", {1.2, 0.0, 0.0});

  auto graph = PeriodicGraphBuilder::buildGraph(cell, 2.0, false);

  EXPECT_EQ(graph.atom_count, 2U);
  EXPECT_EQ(graph.edge_count, 2U); // 0 -> 1 and 1 -> 0

  ASSERT_EQ(graph.atomic_numbers.size(), 2U);
  EXPECT_EQ(graph.atomic_numbers[0], 6);
  EXPECT_EQ(graph.atomic_numbers[1], 8);

  ASSERT_EQ(graph.edge_index_flat.size(), 4U);
  ASSERT_EQ(graph.edge_vectors_flat.size(), 6U);
  ASSERT_EQ(graph.edge_distances.size(), 2U);

  // Check distances
  EXPECT_NEAR(graph.edge_distances[0], static_cast<real_t>(1.2), 1e-5);
  EXPECT_NEAR(graph.edge_distances[1], static_cast<real_t>(1.2), 1e-5);

  // Edge 0: 0 -> 1 (vector +1.2, 0, 0), Edge 1: 1 -> 0 (vector -1.2, 0, 0)
  const int64_t src0 = graph.edge_index_flat[0 * 2 + 0];
  const int64_t dst0 = graph.edge_index_flat[1 * 2 + 0];
  const int64_t src1 = graph.edge_index_flat[0 * 2 + 1];
  const int64_t dst1 = graph.edge_index_flat[1 * 2 + 1];

  EXPECT_TRUE((src0 == 0 && dst0 == 1 && src1 == 1 && dst1 == 0) || (src0 == 1 && dst0 == 0 && src1 == 0 && dst1 == 1));

  if (src0 == 0 && dst0 == 1) {
    EXPECT_NEAR(graph.edge_vectors_flat[0 * 3 + 0], static_cast<real_t>(1.2), 1e-5);
    EXPECT_NEAR(graph.edge_vectors_flat[1 * 3 + 0], static_cast<real_t>(-1.2), 1e-5);
  } else {
    EXPECT_NEAR(graph.edge_vectors_flat[0 * 3 + 0], static_cast<real_t>(-1.2), 1e-5);
    EXPECT_NEAR(graph.edge_vectors_flat[1 * 3 + 0], static_cast<real_t>(1.2), 1e-5);
  }
}

TEST(PeriodicGraphBuilderTests, SelfLoopsToggle) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  // Small cutoff that has no periodic images within 2.0 Å
  auto graph_no_loops = PeriodicGraphBuilder::buildGraph(cell, 2.0, false);
  EXPECT_EQ(graph_no_loops.edge_count, 0U);

  auto graph_with_loops = PeriodicGraphBuilder::buildGraph(cell, 2.0, true);
  EXPECT_EQ(graph_with_loops.edge_count, 1U);
  EXPECT_EQ(graph_with_loops.edge_index_flat[0], 0);
  EXPECT_EQ(graph_with_loops.edge_index_flat[1], 0);
  EXPECT_NEAR(graph_with_loops.edge_distances[0], static_cast<real_t>(0.0), 1e-6);
}

TEST(PeriodicGraphBuilderTests, CutoffEnvelopeMath) {
  const real_t cutoff_radius = 5.0;
  EXPECT_NEAR(PeriodicGraphBuilder::computeCutoffEnvelope(0.0, cutoff_radius), static_cast<real_t>(1.0), 1e-6);
  EXPECT_NEAR(PeriodicGraphBuilder::computeCutoffEnvelope(5.0, cutoff_radius), static_cast<real_t>(0.0), 1e-6);
  EXPECT_NEAR(PeriodicGraphBuilder::computeCutoffEnvelope(6.0, cutoff_radius), static_cast<real_t>(0.0), 1e-6);

  const real_t mid = PeriodicGraphBuilder::computeCutoffEnvelope(2.5, cutoff_radius);
  EXPECT_GT(mid, static_cast<real_t>(0.0));
  EXPECT_LT(mid, static_cast<real_t>(1.0));
}

TEST(PeriodicGraphBuilderTests, BesselBasisFeatures) {
  const real_t cutoff_radius = 5.0;
  const size_t num_basis = 8;

  auto basis_inside = PeriodicGraphBuilder::computeBesselBasis(1.5, cutoff_radius, num_basis);
  EXPECT_EQ(basis_inside.size(), num_basis);
  for (size_t i = 0; i < num_basis; ++i) {
    EXPECT_TRUE(std::isfinite(basis_inside[i]));
  }

  auto basis_outside = PeriodicGraphBuilder::computeBesselBasis(5.5, cutoff_radius, num_basis);
  EXPECT_EQ(basis_outside.size(), num_basis);
  for (size_t i = 0; i < num_basis; ++i) {
    EXPECT_NEAR(basis_outside[i], static_cast<real_t>(0.0), 1e-6);
  }
}

TEST(PeriodicGraphBuilderTests, GaussianRBFExpansion) {
  const size_t num_basis = 5;
  auto rbf = PeriodicGraphBuilder::computeGaussianRBF(1.0, {.start = 0.0, .stop = 4.0, .num_basis = num_basis});
  EXPECT_EQ(rbf.size(), num_basis);
  // Center 1 is at 1.0, so rbf[1] should be exactly 1.0 (diff = 0)
  EXPECT_NEAR(rbf[1], static_cast<real_t>(1.0), 1e-5);
  EXPECT_LT(rbf[0], static_cast<real_t>(1.0));
  EXPECT_LT(rbf[2], static_cast<real_t>(1.0));
}

} // namespace correlation::mlip
