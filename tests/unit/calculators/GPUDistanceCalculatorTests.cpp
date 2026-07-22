/**
 * @file GPUDistanceCalculatorTests.cpp
 * @brief Unit tests for GPU distance calculation in float and double precision.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/GPUDistanceCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::calculators::gpu {

TEST(GPUDistanceCalculatorTests, HasGPUDeviceCheck) { EXPECT_NO_THROW(has_gpu_device()); }

TEST(GPUDistanceCalculatorTests, FloatPrecisionDistanceComputation) {
  correlation::core::Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("O", {1.5, 0.0, 0.0});

  float const cutoff_sq = 4.0F;
  std::vector<std::vector<float>> const bond_cutoffs_sq = {{{4.0F, 4.0F}, {4.0F, 4.0F}}};

  size_t const num_elements = cell.elements().size();
  DistanceTensor out_distances(num_elements, std::vector<std::vector<real_t>>(num_elements));
  correlation::core::NeighborGraph out_graph(2);

  if (has_gpu_device()) {
    EXPECT_NO_THROW(compute_distances_gpu<float>(cell, cutoff_sq, bond_cutoffs_sq, true, out_distances, out_graph));

    ASSERT_EQ(out_distances[0][1].size(), 1);
    EXPECT_NEAR(out_distances[0][1][0], 1.5, 1e-4);

    EXPECT_TRUE(out_graph.areConnected(correlation::core::AtomIndex{0}, correlation::core::AtomIndex{1}));
    const auto &neighbors = out_graph.getNeighbors(0);
    ASSERT_EQ(neighbors.size(), 1);
    EXPECT_EQ(neighbors[0].index, 1);
    EXPECT_NEAR(neighbors[0].distance, 1.5, 1e-4);
  }
}

TEST(GPUDistanceCalculatorTests, DoublePrecisionDistanceComputation) {
  correlation::core::Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.5, 5.0, 5.0});
  cell.addAtom("Si", {9.5, 5.0, 5.0});

  double const cutoff_sq = 4.0;
  std::vector<std::vector<double>> const bond_cutoffs_sq = {{4.0}};

  size_t const num_elements = cell.elements().size();
  DistanceTensor out_distances(num_elements, std::vector<std::vector<real_t>>(num_elements));
  correlation::core::NeighborGraph out_graph(2);

  if (has_gpu_device()) {
    EXPECT_NO_THROW(compute_distances_gpu<double>(cell, cutoff_sq, bond_cutoffs_sq, true, out_distances, out_graph));

    ASSERT_GE(out_distances[0][0].size(), 1);
    EXPECT_NEAR(out_distances[0][0][0], 1.0, 1e-7);

    EXPECT_TRUE(out_graph.areConnected(correlation::core::AtomIndex{0}, correlation::core::AtomIndex{1}));
  }
}

} // namespace correlation::calculators::gpu
