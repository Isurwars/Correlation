// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <benchmark/benchmark.h>
#include "calculators/SteinhardtCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "core/Cell.hpp"

#include <random>

using namespace correlation::calculators;
using namespace correlation::analysis;

/// @brief Benchmarks Steinhardt Bond-Orientational Order Parameters.
///
/// Steinhardt calculations are extremely math-intensive, involving spherical
/// harmonics (Y_l^m) and Wigner 3j symbols for every neighbor pair. This
/// benchmark isolates the Steinhardt computation from the neighbor-search.
static void BM_SteinhardtCalculator(benchmark::State& state) {
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

    SteinhardtCalculator calc;

    for (auto _ : state) {
        calc.calculateFrame(df, settings);
        benchmark::DoNotOptimize(df);
    }

    state.SetComplexityN(state.range(0));
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            n_atoms);
}

BENCHMARK(BM_SteinhardtCalculator)
    ->RangeMultiplier(2)
    ->Range(256, 2048)
    ->Unit(benchmark::kMillisecond)
    ->Complexity(benchmark::oN);
