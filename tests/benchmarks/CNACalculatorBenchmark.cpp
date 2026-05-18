// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <benchmark/benchmark.h>
#include "calculators/CNACalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "core/Cell.hpp"

#include <random>

using namespace correlation::calculators;
using namespace correlation::analysis;

/// @brief Benchmarks Common Neighbor Analysis.
///
/// CNA is computationally demanding because it must iterate over all neighbor
/// pairs of each atom and classify their shared-neighbor topology. This
/// benchmark creates a dense random cell and measures the CNA classification
/// throughput, isolating it from the neighbor-search cost.
static void BM_CNACalculator(benchmark::State& state) {
    correlation::core::Cell cell;
    double L = 30.0;
    cell.setLatticeParameters({L, L, L, 90.0, 90.0, 90.0});

    int n_atoms = state.range(0);
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(0.0, L);

    for (int i = 0; i < n_atoms; ++i) {
        cell.addAtom("Ar", {dist(gen), dist(gen), dist(gen)});
    }

    double cutoff = 3.5;
    std::vector<std::vector<double>> bond_cutoffs_sq = {{cutoff * cutoff}};

    // Pre-compute structure analysis (neighbor lists)
    StructureAnalyzer analyzer(cell, cutoff, bond_cutoffs_sq, false);
    DistributionFunctions df(cell);
    df.setStructureAnalyzer(&analyzer);

    AnalysisSettings settings;

    CNACalculator calc;

    for (auto _ : state) {
        calc.calculateFrame(df, settings);
        benchmark::DoNotOptimize(df);
    }

    state.SetComplexityN(state.range(0));
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            n_atoms);
}

BENCHMARK(BM_CNACalculator)
    ->RangeMultiplier(2)
    ->Range(512, 4096)
    ->Unit(benchmark::kMillisecond)
    ->Complexity(benchmark::oN);
