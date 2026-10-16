/**
 * @file BM_MotifFinder.cpp
 * @brief Microbenchmarks for MotifFinder ring searching algorithms (King & Franzblau).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "CrystalTestHelper.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "calculators/spatial/DistanceCalculator.hpp"
#include "calculators/spatial/MotifFinder.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"
#include "math/Precision.hpp"

#include <benchmark/benchmark.h>

namespace correlation::benchmarks {

namespace {

core::NeighborGraph buildDiamondGraph(int n_cells) {
  // lat_a = 5.43 A, nn distance = 2.35 A
  const auto cell = testing::crystals::createDiamondCell(5.43, "Si", n_cells, n_cells, n_cells);
  const real_t cutoff = 3.0;
  const real_t cutoff_sq = cutoff * cutoff;

  correlation::analysis::BondCutoffMatrix bond_cutoffs(1);
  bond_cutoffs[0][0] = correlation::analysis::BondCutoffRange{
      .min_sq = 0.0,
      .max_sq = 2.8 * 2.8,
  };

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

void BM_MotifFinder_KingRings_Diamond(benchmark::State &state) {
  // Diamond supercell with 8 * 4^3 = 512 atoms
  const auto graph = buildDiamondGraph(4);
  const size_t max_ring_size = 6;

  for (auto _ : state) {
    auto rings =
        calculators::MotifFinder::findRings(graph, max_ring_size, analysis::RingType::King);
    benchmark::DoNotOptimize(rings.size());
    benchmark::ClobberMemory();
  }
  state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(graph.nodeCount()));
}

void BM_MotifFinder_FranzblauRings_Diamond(benchmark::State &state) {
  // Diamond supercell with 8 * 4^3 = 512 atoms
  const auto graph = buildDiamondGraph(4);
  const size_t max_ring_size = 6;

  for (auto _ : state) {
    auto rings =
        calculators::MotifFinder::findRings(graph, max_ring_size, analysis::RingType::Franzblau);
    benchmark::DoNotOptimize(rings.size());
    benchmark::ClobberMemory();
  }
  state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(graph.nodeCount()));
}

void BM_MotifFinder_ExtractCycles_Diamond(benchmark::State &state) {
  // Extracts exact cycle paths
  const auto graph = buildDiamondGraph(4);
  const size_t target_size = 6;

  for (auto _ : state) {
    auto cycles =
        calculators::MotifFinder::extractCycles(graph, target_size, analysis::RingType::King);
    benchmark::DoNotOptimize(cycles.size());
    benchmark::ClobberMemory();
  }
  state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(graph.nodeCount()));
}

} // namespace

BENCHMARK(BM_MotifFinder_KingRings_Diamond);
BENCHMARK(BM_MotifFinder_FranzblauRings_Diamond);
BENCHMARK(BM_MotifFinder_ExtractCycles_Diamond);

} // namespace correlation::benchmarks
