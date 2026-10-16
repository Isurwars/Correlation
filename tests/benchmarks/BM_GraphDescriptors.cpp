/**
 * @file BM_GraphDescriptors.cpp
 * @brief Microbenchmarks for MLIP GraphDescriptors (ring statistics descriptor).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "CrystalTestHelper.hpp"
#include "math/Precision.hpp"
#include "mlip/GraphDescriptors.hpp"
#include "mlip/PeriodicGraphBuilder.hpp"

#include <benchmark/benchmark.h>

namespace correlation::benchmarks {

namespace {

void BM_GraphDescriptors_RingStatistics(benchmark::State &state) {
  // Diamond supercell with 8 * 3^3 = 216 atoms
  const auto cell = testing::crystals::createDiamondCell(5.43, "Si", 3, 3, 3);
  const auto graph = mlip::PeriodicGraphBuilder::buildGraph(cell, static_cast<real_t>(2.8));
  const size_t max_ring_size = 6;

  for (auto _ : state) {
    auto desc = mlip::GraphDescriptors::computeRingStatisticsDescriptor(graph, max_ring_size);
    benchmark::DoNotOptimize(desc.size());
    benchmark::ClobberMemory();
  }
  state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(graph.atom_count));
}

} // namespace

BENCHMARK(BM_GraphDescriptors_RingStatistics);

} // namespace correlation::benchmarks
