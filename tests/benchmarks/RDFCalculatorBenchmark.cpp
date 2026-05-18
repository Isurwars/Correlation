// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <benchmark/benchmark.h>
#include "calculators/RDFCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "core/Cell.hpp"

#include <random>

using namespace correlation::calculators;
using namespace correlation::analysis;

/// @brief Benchmarks the Radial Distribution Function computation.
///
/// Sets up a cell with random Argon atoms, pre-computes the StructureAnalyzer
/// (neighbor lists), and then times only the RDF histogram assembly via
/// RDFCalculator::calculateFrame. This isolates the binning and normalisation
/// cost from the neighbor-search overhead.
static void BM_RDFCalculator(benchmark::State& state) {
    correlation::core::Cell cell;
    double L = 40.0;
    cell.setLatticeParameters({L, L, L, 90.0, 90.0, 90.0});

    int n_atoms = state.range(0);
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(0.0, L);

    for (int i = 0; i < n_atoms; ++i) {
        cell.addAtom("Ar", {dist(gen), dist(gen), dist(gen)});
    }

    double cutoff = 10.0;
    std::vector<std::vector<double>> bond_cutoffs_sq = {{cutoff * cutoff}};

    // Pre-compute structure analysis (neighbor lists)
    StructureAnalyzer analyzer(cell, cutoff, bond_cutoffs_sq, false);
    DistributionFunctions df(cell);
    df.setStructureAnalyzer(&analyzer);

    AnalysisSettings settings;
    settings.r_max = cutoff;
    settings.r_bin_width = 0.02;

    RDFCalculator calc;

    for (auto _ : state) {
        calc.calculateFrame(df, settings);
        benchmark::DoNotOptimize(df);
    }

    state.SetComplexityN(state.range(0));
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            n_atoms);
}

BENCHMARK(BM_RDFCalculator)
    ->RangeMultiplier(2)
    ->Range(512, 4096)
    ->Unit(benchmark::kMillisecond)
    ->Complexity(benchmark::oN);
