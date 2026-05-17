// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <benchmark/benchmark.h>
#include "analysis/StructureAnalyzer.hpp"
#include "core/Cell.hpp"

#include <random>

using namespace correlation::analysis;

static void BM_StructureAnalyzer(benchmark::State& state) {
    correlation::core::Cell cell;
    double L = 30.0;
    cell.setLatticeParameters({L, L, L, 90.0, 90.0, 90.0});
    
    // Add atoms
    int n_atoms = state.range(0);
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(0.0, L);
    
    for (int i = 0; i < n_atoms; ++i) {
        cell.addAtom("Ar", {dist(gen), dist(gen), dist(gen)});
    }
    
    double cutoff = 3.5; 
    std::vector<std::vector<double>> bond_cutoffs = {{cutoff * cutoff}};
    
    for (auto _ : state) {
        // This measures the time it takes to instantiate StructureAnalyzer
        // which now heavily relies on tbb::task_group for parallelism
        StructureAnalyzer analyzer(cell, cutoff, bond_cutoffs, false);
        benchmark::DoNotOptimize(analyzer);
    }
    
    state.SetComplexityN(state.range(0));
}

// Register the benchmark
BENCHMARK(BM_StructureAnalyzer)
    ->RangeMultiplier(2)
    ->Range(512, 4096)
    ->Unit(benchmark::kMillisecond)
    ->Complexity(benchmark::oN);
