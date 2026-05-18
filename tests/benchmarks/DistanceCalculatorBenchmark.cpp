// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <benchmark/benchmark.h>
#include "calculators/DistanceCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include <random>

using namespace correlation::calculators;

/// @brief Benchmarks the core pairwise distance computation with PBC.
///
/// This measures DistanceCalculator::compute, the lowest-level hot-path that
/// drives almost every structural analysis calculator. The benchmark sweeps
/// atom counts from 512 to 8192 to characterise the cell-list performance
/// under increasing density.
static void BM_DistanceCalculator(benchmark::State& state) {
    correlation::core::Cell cell;
    double L = 40.0;
    cell.setLatticeParameters({L, L, L, 90.0, 90.0, 90.0});

    int n_atoms = state.range(0);
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(0.0, L);

    for (int i = 0; i < n_atoms; ++i) {
        cell.addAtom("Ar", {dist(gen), dist(gen), dist(gen)});
    }

    double cutoff = 6.0;
    double cutoff_sq = cutoff * cutoff;
    std::vector<std::vector<double>> bond_cutoffs_sq = {{cutoff_sq}};
    size_t num_elements = cell.elements().size();

    for (auto _ : state) {
        DistanceTensor out_distances(
            num_elements, std::vector<std::vector<double>>(num_elements));
        correlation::core::NeighborGraph out_graph(cell.atoms().size());

        DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs_sq, false,
                                    out_distances, out_graph);

        benchmark::DoNotOptimize(out_distances);
        benchmark::DoNotOptimize(out_graph);
    }

    state.SetComplexityN(state.range(0));
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            n_atoms);
}

BENCHMARK(BM_DistanceCalculator)
    ->RangeMultiplier(2)
    ->Range(512, 8192)
    ->Unit(benchmark::kMillisecond)
    ->Complexity(benchmark::oNSquared);
