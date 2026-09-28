// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

/**
 * @file ThreadDeterminismTests.cpp
 * @brief Multi-threaded determinism tests verifying oneTBB parallel consistency across thread
 * counts.
 */

#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/order/SteinhardtCalculator.hpp"
#include "calculators/scattering/StructureFactorCalculator.hpp"
#include "calculators/spatial/CNCalculator.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/Precision.hpp"

#include "../../CrystalTestHelper.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <map>
#include <string>
#include <tbb/task_arena.h>
#include <utility>
#include <vector>

namespace correlation::testing {

namespace {

using correlation::analysis::DistributionFunctions;
using correlation::analysis::StructureAnalyzer;
using correlation::calculators::CNCalculator;
using correlation::calculators::SteinhardtCalculator;
using correlation::calculators::StructureFactorCalculator;
using correlation::core::Cell;
using correlation::core::Trajectory;

/// Executes a callable inside an isolated tbb::task_arena with a specified thread count.
template <typename Func> decltype(auto) runInArena(int num_threads, Func &&func) {
  tbb::task_arena arena(num_threads);
  return arena.execute(std::forward<Func>(func));
}

/// Asserts that two histogram partials are numerically identical within allowed relative tolerance.
void assertHistogramsMatch(const std::vector<real_t> &baseline, const std::vector<real_t> &target,
                           double max_rel_err) {
  ASSERT_EQ(baseline.size(), target.size());
  for (std::size_t i = 0; i < baseline.size(); ++i) {
    auto const baseline_val = static_cast<double>(baseline[i]);
    auto const target_val = static_cast<double>(target[i]);
    double const abs_err = std::abs(baseline_val - target_val);
    double const denom = std::max(std::abs(baseline_val), 1e-12);
    double const rel_err = abs_err / denom;
    EXPECT_LE(rel_err, max_rel_err) << "Bin " << i << ": Baseline=" << baseline_val
                                    << ", Target=" << target_val << ", RelError=" << rel_err;
  }
}

/// Asserts that two discrete histogram partials are bitwise identical across threads.
void assertDiscreteHistogramsEqual(const std::vector<real_t> &baseline,
                                   const std::vector<real_t> &target, int thread_count) {
  ASSERT_EQ(baseline.size(), target.size());
  for (std::size_t i = 0; i < baseline.size(); ++i) {
    EXPECT_DOUBLE_EQ(static_cast<double>(target[i]), static_cast<double>(baseline[i]))
        << "Mismatch at bin " << i << " with " << thread_count << " threads";
  }
}

class ThreadDeterminismTests : public ::testing::Test {
protected:
  static constexpr double REDUCTION_REL_TOL = 1e-6;
  std::vector<int> thread_counts_{1, 2, 4};
};

} // namespace

TEST_F(ThreadDeterminismTests, CoordinationNumberDeterminismAcrossThreads) {
  // Dense BCC Iron supercell (3x3x3 unit cells = 54 atoms)
  constexpr auto LAT_A = static_cast<real_t>(2.8665);
  Cell cell = crystals::createBCCCell(LAT_A, "Fe", 3, 3, 3);

  constexpr auto CUTOFF = static_cast<real_t>(3.1);

  std::vector<real_t> baseline_cn;

  for (const int threads : thread_counts_) {
    Trajectory traj;
    traj.addFrame(cell);
    traj.precomputeBondCutoffs();

    auto hist = runInArena(threads, [&] {
      StructureAnalyzer const analyzer(cell, CUTOFF, traj.getBondCutoffsSQ());
      return CNCalculator::calculate(cell, &analyzer);
    });

    ASSERT_TRUE(hist.partials.contains("Fe-Fe"));
    const auto &partial = hist.partials.at("Fe-Fe");

    if (baseline_cn.empty()) {
      baseline_cn = partial;
    } else {
      assertDiscreteHistogramsEqual(baseline_cn, partial, threads);
    }
  }
}

TEST_F(ThreadDeterminismTests, RDFDeterminismAcrossThreads) {
  // FCC Copper supercell (3x3x3 unit cells = 108 atoms)
  constexpr auto LAT_A = static_cast<real_t>(3.615);
  Cell cell = crystals::createFCCCell(LAT_A, "Cu", 3, 3, 3);

  constexpr auto MAX_R = static_cast<real_t>(5.0);

  std::vector<real_t> baseline_gr;

  for (const int threads : thread_counts_) {
    Trajectory traj;
    traj.addFrame(cell);
    traj.precomputeBondCutoffs();

    DistributionFunctions dists(cell, MAX_R, traj.getBondCutoffsSQ());

    runInArena(threads, [&] {
      dists.calculateRDF({.r_max = MAX_R, .r_bin_width = static_cast<real_t>(0.02)});
    });

    const auto &hist = dists.getHistogram("g_r");
    ASSERT_TRUE(hist.partials.contains("Cu-Cu"));
    const auto &partial = hist.partials.at("Cu-Cu");

    if (baseline_gr.empty()) {
      baseline_gr = partial;
    } else {
      assertHistogramsMatch(baseline_gr, partial, REDUCTION_REL_TOL);
    }
  }
}

TEST_F(ThreadDeterminismTests, PADDeterminismAcrossThreads) {
  // Diamond Silicon supercell (2x2x2 = 64 atoms)
  constexpr auto LAT_A = static_cast<real_t>(5.43);
  Cell cell = crystals::createDiamondCell(LAT_A, "Si", 2, 2, 2);

  constexpr auto CUTOFF = static_cast<real_t>(2.60);

  std::vector<real_t> baseline_pad;

  for (const int threads : thread_counts_) {
    Trajectory traj;
    traj.addFrame(cell);
    traj.precomputeBondCutoffs();

    DistributionFunctions dists(cell, CUTOFF, traj.getBondCutoffsSQ());

    runInArena(threads, [&] { dists.calculatePAD(0.2); });

    const auto &hist = dists.getHistogram("PAD");
    ASSERT_TRUE(hist.partials.contains("Si-Si-Si"));
    const auto &partial = hist.partials.at("Si-Si-Si");

    if (baseline_pad.empty()) {
      baseline_pad = partial;
    } else {
      assertHistogramsMatch(baseline_pad, partial, REDUCTION_REL_TOL);
    }
  }
}

TEST_F(ThreadDeterminismTests, StructureFactorDeterminismAcrossThreads) {
  // BCC Iron supercell (3x3x3 unit cells = 54 atoms)
  constexpr auto LAT_A = static_cast<real_t>(2.8665);
  Cell cell = crystals::createBCCCell(LAT_A, "Fe", 3, 3, 3);

  std::vector<real_t> baseline_sq;

  for (const int threads : thread_counts_) {
    DistributionFunctions dists(cell);
    StructureFactorCalculator calc;
    analysis::AnalysisSettings settings;
    settings.q_max = static_cast<real_t>(6.0);
    settings.q_bin_width = static_cast<real_t>(0.05);

    runInArena(threads, [&] { calc.calculateFrame(dists, settings); });

    const auto &hist = dists.getHistogram("S_q");
    ASSERT_TRUE(hist.partials.contains("Total"));
    const auto &partial = hist.partials.at("Total");

    if (baseline_sq.empty()) {
      baseline_sq = partial;
    } else {
      assertHistogramsMatch(baseline_sq, partial, REDUCTION_REL_TOL);
    }
  }
}

TEST_F(ThreadDeterminismTests, SteinhardtOrderParametersDeterminismAcrossThreads) {
  // FCC Copper supercell (2x2x2 = 32 atoms)
  constexpr auto LAT_A = static_cast<real_t>(3.615);
  Cell cell = crystals::createFCCCell(LAT_A, "Cu", 2, 2, 2);

  constexpr auto CUTOFF = static_cast<real_t>(2.8);

  std::vector<real_t> baseline_q6;

  for (const int threads : thread_counts_) {
    Trajectory traj;
    traj.addFrame(cell);
    traj.precomputeBondCutoffs();

    auto hists = runInArena(threads, [&] {
      StructureAnalyzer const analyzer(cell, CUTOFF, traj.getBondCutoffsSQ());
      return SteinhardtCalculator::calculate(cell, &analyzer);
    });

    ASSERT_TRUE(hists.contains("Q6"));
    const auto &q6_data = hists.at("Q6").partials.at("Total");

    if (baseline_q6.empty()) {
      baseline_q6 = q6_data;
    } else {
      assertHistogramsMatch(baseline_q6, q6_data, REDUCTION_REL_TOL);
    }
  }
}

} // namespace correlation::testing
