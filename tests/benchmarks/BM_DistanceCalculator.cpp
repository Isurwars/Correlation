/**
 * @file BM_DistanceCalculator.cpp
 * @brief Microbenchmarks for pairwise DistanceCalculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "CrystalTestHelper.hpp"
#include "calculators/spatial/DistanceCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"
#include "math/Precision.hpp"

#include <benchmark/benchmark.h>

namespace correlation::benchmarks {

namespace {

correlation::analysis::BondCutoffMatrix makeCutoffMatrix(real_t min_sq, real_t max_sq) {
  correlation::analysis::BondCutoffMatrix matrix(1);
  matrix[0][0] = correlation::analysis::BondCutoffRange{
      .min_sq = min_sq,
      .max_sq = max_sq,
  };
  return matrix;
}

void benchDistanceSmallCell(benchmark::State &state) {
  // 3x3x3 SC cell (27 atoms) with large cutoff triggering small-cell multi-image branch (k=1)
  const auto cell = testing::crystals::createSimpleCubicCell(4.0, "Ar", 3, 3, 3);
  const real_t cutoff = 6.0;
  const real_t cutoff_sq = cutoff * cutoff;
  const auto bond_cutoffs = makeCutoffMatrix(0.0, 4.5 * 4.5);
  const calculators::DistanceCalculationConfig config{
      .r_max = cutoff,
      .r_bin_width = static_cast<real_t>(0.05),
      .num_bins = static_cast<size_t>(cutoff / 0.05),
  };

  for ([[maybe_unused]] auto iter : state) {
    core::NeighborGraph graph(cell.atomCount());
    calculators::RawHistogramTensor histograms;
    calculators::DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, graph,
                                             &histograms, config);
    benchmark::DoNotOptimize(graph.nodeCount());
    benchmark::DoNotOptimize(histograms.size());
    benchmark::ClobberMemory();
  }
  state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(cell.atomCount()));
}

void benchDistanceMidCell(benchmark::State &state) {
  // 6x6x6 FCC supercell (864 atoms)
  const auto cell = testing::crystals::createFCCCell(3.61, "Cu", 6, 6, 6);
  const real_t cutoff = 5.0;
  const real_t cutoff_sq = cutoff * cutoff;
  const auto bond_cutoffs = makeCutoffMatrix(0.0, 3.0 * 3.0);
  const calculators::DistanceCalculationConfig config{
      .r_max = cutoff,
      .r_bin_width = static_cast<real_t>(0.05),
      .num_bins = static_cast<size_t>(cutoff / 0.05),
  };

  for ([[maybe_unused]] auto iter : state) {
    core::NeighborGraph graph(cell.atomCount());
    calculators::RawHistogramTensor histograms;
    calculators::DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, graph,
                                             &histograms, config);
    benchmark::DoNotOptimize(graph.nodeCount());
    benchmark::DoNotOptimize(histograms.size());
    benchmark::ClobberMemory();
  }
  state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(cell.atomCount()));
}

void benchDistanceLargeCell(benchmark::State &state) {
  // 14x14x14 FCC supercell (10,976 atoms)
  const auto cell = testing::crystals::createFCCCell(3.61, "Cu", 14, 14, 14);
  const real_t cutoff = 5.0;
  const real_t cutoff_sq = cutoff * cutoff;
  const auto bond_cutoffs = makeCutoffMatrix(0.0, 3.0 * 3.0);
  const calculators::DistanceCalculationConfig config{
      .r_max = cutoff,
      .r_bin_width = static_cast<real_t>(0.05),
      .num_bins = static_cast<size_t>(cutoff / 0.05),
  };

  for ([[maybe_unused]] auto iter : state) {
    core::NeighborGraph graph(cell.atomCount());
    calculators::RawHistogramTensor histograms;
    calculators::DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, graph,
                                             &histograms, config);
    benchmark::DoNotOptimize(graph.nodeCount());
    benchmark::DoNotOptimize(histograms.size());
    benchmark::ClobberMemory();
  }
  state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(cell.atomCount()));
}

} // namespace

BENCHMARK(benchDistanceSmallCell);
BENCHMARK(benchDistanceMidCell);
BENCHMARK(benchDistanceLargeCell);

} // namespace correlation::benchmarks
