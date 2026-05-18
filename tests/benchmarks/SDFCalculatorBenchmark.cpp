// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <benchmark/benchmark.h>
#include "calculators/SDFCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "core/Cell.hpp"

#include <random>

using namespace correlation::calculators;
using namespace correlation::analysis;

/// @brief Benchmarks the Spatial Distribution Function (3D grid binning).
///
/// The SDF calculator allocates a 3D density grid and bins every atom into it.
/// This benchmark exercises the grid allocation, PBC wrapping, and density
/// accumulation for increasing atom counts to measure memory throughput and
/// binning efficiency.
static void BM_SDFCalculator(benchmark::State& state) {
    correlation::core::Cell cell;
    double L = 30.0;
    cell.setLatticeParameters({L, L, L, 90.0, 90.0, 90.0});

    int n_atoms = state.range(0);
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(0.0, L);

    for (int i = 0; i < n_atoms; ++i) {
        cell.addAtom("Ar", {dist(gen), dist(gen), dist(gen)});
    }

    DistributionFunctions df(cell);
    AnalysisSettings settings;
    settings.r_bin_width = 0.5; // Grid resolution for SDF

    SDFCalculator calc;

    for (auto _ : state) {
        calc.calculateFrame(df, settings);
        benchmark::DoNotOptimize(df);
    }

    state.SetComplexityN(state.range(0));
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            n_atoms);
}

BENCHMARK(BM_SDFCalculator)
    ->RangeMultiplier(2)
    ->Range(1024, 8192)
    ->Unit(benchmark::kMillisecond)
    ->Complexity(benchmark::oN);
