// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <benchmark/benchmark.h>
#include "calculators/ClusterCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "core/Cell.hpp"

#include <random>

using namespace correlation::calculators;
using namespace correlation::analysis;

static void BM_ClusterCalculator(benchmark::State& state) {
    correlation::core::Cell cell;
    double L = 50.0;
    cell.setLatticeParameters({L, L, L, 90.0, 90.0, 90.0});
    
    // Add atoms
    int n_atoms = state.range(0);
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(0.0, L);
    
    for (int i = 0; i < n_atoms; ++i) {
        cell.addAtom("Ar", {dist(gen), dist(gen), dist(gen)});
    }
    
    DistributionFunctions df(cell);
    AnalysisSettings settings;
    
    double cutoff = 2.5; // Enough to create some clusters
    StructureAnalyzer analyzer(cell, cutoff, {{cutoff * cutoff}}, false);
    df.setStructureAnalyzer(&analyzer);
    
    ClusterCalculator calc;
    
    for (auto _ : state) {
        calc.calculateFrame(df, settings);
        benchmark::DoNotOptimize(df);
    }
    
    state.SetComplexityN(state.range(0));
}

// Register the benchmark
BENCHMARK(BM_ClusterCalculator)
    ->RangeMultiplier(2)
    ->Range(1024, 8192)
    ->Complexity(benchmark::oN);
