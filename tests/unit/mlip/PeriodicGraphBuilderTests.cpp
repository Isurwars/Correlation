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

    const real_t s_x = graph.edge_shifts_flat[edge * static_cast<size_t>(3 + 0)];
    const real_t s_y = graph.edge_shifts_flat[edge * static_cast<size_t>(3 + 1)];
    const real_t s_z = graph.edge_shifts_flat[edge * static_cast<size_t>(3 + 2)];
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
  // Edge 0: 0 -> 1, Edge 1: 1 -> 0 (or vice versa)
  const int64_t src0 = graph.edge_index_flat[0 * 2 + 0];
  const int64_t dst0 = graph.edge_index_flat[1 * 2 + 0];
  const int64_t src1 = graph.edge_index_flat[0 * 2 + 1];
  const int64_t dst1 = graph.edge_index_flat[1 * 2 + 1];

  EXPECT_TRUE((src0 == 0 && dst0 == 1 && src1 == 1 && dst1 == 0) || (src0 == 1 && dst0 == 0 && src1 == 0 && dst1 == 1));
}

TEST(TorchGNNModelTests, InterfaceAndAccessors) {
  TorchGNNModel model("non_existent_model.pt", "cpu", 6.0);

  EXPECT_EQ(model.cutoff(), static_cast<real_t>(6.0));
  model.setCutoff(4.5);
  EXPECT_EQ(model.cutoff(), static_cast<real_t>(4.5));
  EXPECT_TRUE(model.getModelName().contains("non_existent_model.pt"));

#ifndef CORRELATION_HAS_LIBTORCH
  EXPECT_FALSE(model.isLoaded());
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  auto output = model.evaluate(cell);
  EXPECT_EQ(output.forces.size(), 1U);
#endif
}

} // namespace correlation::mlip
