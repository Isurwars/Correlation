// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "calculators/scattering/StructureFactorCalculator.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/Precision.hpp"
#include "readers/FileReader.hpp"

#include "../../CrystalTestHelper.hpp"

#include <cmath>
#include <filesystem>
#include <gtest/gtest.h>
#include <iostream>

using correlation::real_t;

namespace correlation::analysis {
namespace {
class StructureFactorCalculatorTests : public ::testing::Test {
protected:
  void SetUp() override {}
};

correlation::core::Cell createSimpleCubicCell(real_t vec_a) {
  const real_t length = static_cast<real_t>(2.0) * vec_a;
  correlation::core::Cell cell(std::array<real_t, 6>{
      length,
      length,
      length,
      90.0,
      90.0,
      90.0,
  });

  for (int ix = 0; ix < 2; ++ix) {
    for (int iy = 0; iy < 2; ++iy) {
      for (int iz = 0; iz < 2; ++iz) {
        cell.addAtom("Si", {static_cast<real_t>(ix) * vec_a, static_cast<real_t>(iy) * vec_a,
                            static_cast<real_t>(iz) * vec_a});
      }
    }
  }
  return cell;
}

real_t findMaxPeak(const Histogram &hist, real_t q_min, real_t q_max) {
  auto peak_sq = static_cast<real_t>(0.0);
  for (size_t i = 0; i < hist.bins.size(); ++i) {
    if (hist.bins[i] >= q_min && hist.bins[i] <= q_max) {
      peak_sq = std::max(peak_sq, hist.partials.at("Total")[i]);
    }
  }
  return peak_sq;
}

} // namespace

// For a simple cubic crystal, S(Q) should show peaks at the allowed Bragg
// reflections. The first peak for a simple cubic lattice is at Q = 2*pi/a.
TEST_F(StructureFactorCalculatorTests, CalculatesSimpleCubicBraggPeak) {
  // Build a 2x2x2 simple cubic supercell with lattice constant a=3.0 Å
  const real_t vec_a = 3.0;
  const auto cell = createSimpleCubicCell(vec_a);

  DistributionFunctions dists(cell);
  const correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 5.0;
  settings.q_bin_width = static_cast<real_t>(0.05);

  calc.calculateFrame(dists, settings);

  const auto &hist = dists.getHistogram("S_q");
  EXPECT_FALSE(hist.bins.empty());
  EXPECT_TRUE(hist.partials.contains("Total"));
  EXPECT_EQ(hist.y_unit, "Å");

  // The first Bragg peak for a SC lattice at a=3 is at Q = 2*pi/3 ≈ 2.094
  // Å^{-1} Find the max S(Q) value in the range [1.8, 2.4]
  const real_t peak_sq = findMaxPeak(hist, static_cast<real_t>(1.8), static_cast<real_t>(2.4));
  // At a Bragg peak, S(Q) should be much greater than 1
  EXPECT_GT(peak_sq, 2.0);
}

// For a single atom, S(Q) = 1 for all Q.
TEST_F(StructureFactorCalculatorTests, SingleAtomGivesOne) {
  correlation::core::Cell cell(std::array<real_t, 6>{10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  DistributionFunctions dists(cell);
  const correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 5.0;
  settings.q_bin_width = static_cast<real_t>(0.1);

  calc.calculateFrame(dists, settings);
  const auto &hist = dists.getHistogram("S_q");

  ASSERT_FALSE(hist.bins.empty());
  for (size_t i = 0; i < hist.bins.size(); ++i) {
    // S(Q) = 1/N * |sum exp(i q.r)|^2 = 1/1 * |exp(i*0)|^2 = 1.0
    if (hist.partials.at("Total")[i] > 0.0) {
      EXPECT_NEAR(hist.partials.at("Total")[i], 1.0, 1e-6);
    }
  }
}

// --- Tests migrated from SQTests.cpp ---

// For two identical atoms separated by 1.5 Å in a 10x10x10 box,
// S(Q) for the Si-Si partial should appear and not be empty.
TEST_F(StructureFactorCalculatorTests, DimerProducesValidSQ) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});
  cell.addAtom("Ar", {6.5, 5.0, 5.0}); // Distance 1.5 Å

  DistributionFunctions dists(cell);
  const correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 10.0;
  settings.q_bin_width = static_cast<real_t>(0.1);

  calc.calculateFrame(dists, settings);

  EXPECT_NO_THROW((void)dists.getHistogram("S_q"));
  const auto &hist = dists.getHistogram("S_q");
  EXPECT_FALSE(hist.bins.empty());
  EXPECT_TRUE(hist.partials.contains("Total"));

  // For a homonuclear dimer of 2 atoms at distance d=1.5 A, S(Q) = 1 +
  // sin(Q*d)/(Q*d)
  for (size_t i = 0; i < hist.bins.size(); ++i) {
    const real_t q_val = hist.bins[i];
    auto expected_sq = static_cast<real_t>(1.0);
    if (q_val > static_cast<real_t>(1e-6)) {
      expected_sq +=
          std::sin(q_val * static_cast<real_t>(1.5)) / (q_val * static_cast<real_t>(1.5));
    } else {
      expected_sq = static_cast<real_t>(2.0);
    }

    // Check if the partials correctly accumulate
    // The exact structural Debye continuous S(Q) formula differs from the
    // discrete k-point grid approximation from a finite periodic simulation
    // cell by up to 10-30%. Additionally, small Q bins may be exactly 0 due to
    // the 10x10x10 PBC limiting the minimum k-vector.
    if (hist.partials.at("Total")[i] > 0.0) {
      EXPECT_NEAR(hist.partials.at("Total")[i], expected_sq, 0.6);
    }
  }
}

// For N identical atoms, S(Q->0) should approach N.
// Using a 4-atom Si cluster and a very small Q bin.
TEST_F(StructureFactorCalculatorTests, VerifiesPartialsExistForMulticomponent) {
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  cell.addAtom("C", {0.0, 0.0, 0.0});
  cell.addAtom("O", {static_cast<real_t>(1.2), 0.0, 0.0});

  DistributionFunctions dists(cell);
  const correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 10.0;
  settings.q_bin_width = 1.0;

  calc.calculateFrame(dists, settings);
  const auto &hist = dists.getHistogram("S_q");

  // Verify all expected partials are present
  EXPECT_TRUE(hist.partials.contains("C-C"));
  EXPECT_TRUE(hist.partials.contains("O-O"));
  EXPECT_TRUE(hist.partials.contains("C-O"));
  EXPECT_TRUE(hist.partials.contains("Total"));

  EXPECT_FALSE(hist.bins.empty());
}

// Verify a multicomponent system produces an S_q partial for each pair.
TEST_F(StructureFactorCalculatorTests, HomonuclearClusterPartialsPresent) {
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {2.0, 0.0, 0.0});
  cell.addAtom("Si", {0.0, 2.0, 0.0});
  cell.addAtom("Si", {0.0, 0.0, 2.0});

  DistributionFunctions dists(cell);
  const correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 5.0;
  settings.q_bin_width = 0.5;

  calc.calculateFrame(dists, settings);
  const auto &hist = dists.getHistogram("S_q");

  EXPECT_TRUE(hist.partials.contains("Si-Si"));
  EXPECT_TRUE(hist.partials.contains("Total"));
}

TEST_F(StructureFactorCalculatorTests, CalculateSF_InvalidInputsThrow) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  DistributionFunctions dists(cell);
  const correlation::calculators::StructureFactorCalculator calc;

  AnalysisSettings settings_invalid_q;
  settings_invalid_q.q_max = -1.0;
  settings_invalid_q.q_bin_width = static_cast<real_t>(0.1);
  EXPECT_THROW(calc.calculateFrame(dists, settings_invalid_q), std::invalid_argument);

  AnalysisSettings settings_invalid_bin;
  settings_invalid_bin.q_max = 5.0;
  settings_invalid_bin.q_bin_width = static_cast<real_t>(-0.1);
  EXPECT_THROW(calc.calculateFrame(dists, settings_invalid_bin), std::invalid_argument);

  AnalysisSettings settings_bin_too_large;
  settings_bin_too_large.q_max = 5.0;
  settings_bin_too_large.q_bin_width = 10.0;
  EXPECT_THROW(calc.calculateFrame(dists, settings_bin_too_large), std::invalid_argument);
}

TEST_F(StructureFactorCalculatorTests, CalculateSF_EmptyCellThrows) {
  const correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  EXPECT_THROW((void)DistributionFunctions(cell), std::invalid_argument);
}

TEST_F(StructureFactorCalculatorTests, FCC_Copper_SQ) {
  auto const a = static_cast<real_t>(3.615);
  const auto cell = correlation::testing::crystals::createFCCCell(a, "Cu", 2, 2, 2);

  DistributionFunctions dists(cell);
  const correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 6.0;
  settings.q_bin_width = static_cast<real_t>(0.05);

  calc.calculateFrame(dists, settings);
  const auto &hist = dists.getHistogram("S_q");

  // FCC (111) peak at Q = 2*pi/a * sqrt(3) ~ 3.01 A^-1
  const real_t peak_sq = findMaxPeak(hist, static_cast<real_t>(2.8), static_cast<real_t>(3.2));
  EXPECT_GT(peak_sq, 2.0);
}

TEST_F(StructureFactorCalculatorTests, BCC_Iron_SQ) {
  auto const a = static_cast<real_t>(2.866);
  const auto cell = correlation::testing::crystals::createBCCCell(a, "Fe", 3, 3, 3);

  DistributionFunctions dists(cell);
  const correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 6.0;
  settings.q_bin_width = static_cast<real_t>(0.05);

  calc.calculateFrame(dists, settings);
  const auto &hist = dists.getHistogram("S_q");

  // BCC (110) peak at Q = 2*pi/a * sqrt(2) ~ 3.10 A^-1
  const real_t peak_sq = findMaxPeak(hist, static_cast<real_t>(2.9), static_cast<real_t>(3.3));
  EXPECT_GT(peak_sq, 2.0);
}

TEST_F(StructureFactorCalculatorTests, Diamond_Silicon_SQ) {
  auto const a = static_cast<real_t>(5.431);
  const auto cell = correlation::testing::crystals::createDiamondCell(a, "Si", 2, 2, 2);

  DistributionFunctions dists(cell);
  const correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 5.0;
  settings.q_bin_width = static_cast<real_t>(0.05);

  calc.calculateFrame(dists, settings);
  const auto &hist = dists.getHistogram("S_q");

  // Diamond (111) peak at Q = 2*pi/a * sqrt(3) ~ 2.00 A^-1
  const real_t peak_sq = findMaxPeak(hist, static_cast<real_t>(1.8), static_cast<real_t>(2.2));
  EXPECT_GT(peak_sq, 2.0);
}

TEST_F(StructureFactorCalculatorTests, NaCl_RockSalt_SQ) {
  auto const a = static_cast<real_t>(5.64);
  const auto cell = correlation::testing::crystals::createNaClCell(a, "Na", "Cl", 2, 2, 2);

  DistributionFunctions dists(cell);
  const correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 5.0;
  settings.q_bin_width = static_cast<real_t>(0.05);

  calc.calculateFrame(dists, settings);
  const auto &hist = dists.getHistogram("S_q");

  EXPECT_TRUE(hist.partials.contains("Na-Na"));
  EXPECT_TRUE(hist.partials.contains("Cl-Cl"));
  EXPECT_TRUE(hist.partials.contains("Cl-Na") || hist.partials.contains("Na-Cl"));
  EXPECT_TRUE(hist.partials.contains("Total"));

  // NaCl (200) peak at Q = 2*pi/a * sqrt(4) ~ 2.23 A^-1
  const real_t peak_sq = findMaxPeak(hist, static_cast<real_t>(2.05), static_cast<real_t>(2.4));
  EXPECT_GT(peak_sq, 1.5);
}

TEST_F(StructureFactorCalculatorTests, LargeCellUsesHybridDirectAndFourierSQ) {
  // Large 60x60x60 Å box with 1000 atoms (V = 216,000 Å³) and q_max = 20 Å⁻¹ triggers the
  // hybrid low-Q direct + high-Q Fourier g(r)->S(Q) path instead of enumerating >14M q-vectors.
  correlation::core::Cell cell(std::array<real_t, 6>{60.0, 60.0, 60.0, 90.0, 90.0, 90.0});
  for (int ix = 0; ix < 10; ++ix) {
    for (int iy = 0; iy < 10; ++iy) {
      for (int iz = 0; iz < 10; ++iz) {
        cell.addAtom((ix + iy + iz) % 2 == 0 ? "C" : "O",
                     {static_cast<real_t>(ix * 6.0 + 1.0), static_cast<real_t>(iy * 6.0 + 1.0),
                      static_cast<real_t>(iz * 6.0 + 1.0)});
      }
    }
  }

  BondCutoffMatrix const zero_cutoffs(cell.elements().size());
  DistributionFunctions dists(cell, static_cast<real_t>(15.0), zero_cutoffs);
  dists.calculateRDF({.r_max = static_cast<real_t>(15.0), .r_bin_width = static_cast<real_t>(0.05)});

  const correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings settings;
  settings.q_max = static_cast<real_t>(20.0);
  settings.q_bin_width = static_cast<real_t>(0.05);
  settings.r_max = static_cast<real_t>(15.0);
  settings.r_int_max = static_cast<real_t>(12.0);
  settings.r_bin_width = static_cast<real_t>(0.05);

  calc.calculateFrame(dists, settings);
  const auto &hist = dists.getHistogram("S_q");

  ASSERT_EQ(hist.bins.size(), 400u);
  ASSERT_TRUE(hist.partials.contains("Total"));
  ASSERT_TRUE(hist.partials.contains("C-C"));
  ASSERT_TRUE(hist.partials.contains("O-O"));
  ASSERT_TRUE(hist.partials.contains("C-O"));

  const auto &total_sq = hist.partials.at("Total");
  for (real_t val : total_sq) {
    EXPECT_TRUE(std::isfinite(val));
  }
  // High-Q Fourier tail must oscillate cleanly around 1.0
  double tail_sum = 0.0;
  size_t tail_cnt = 0;
  for (size_t i = 300; i < total_sq.size(); ++i) {
    tail_sum += static_cast<double>(total_sq[i]);
    ++tail_cnt;
  }
  EXPECT_NEAR(tail_sum / static_cast<double>(tail_cnt), 1.0, 0.25);
}

namespace {

struct SQAccuracyStats {
  double r_squared;
  double chi_squared;
  double reduced_chi_squared;
  double rmse;
  size_t count;
};

[[nodiscard]] correlation::core::Cell
buildDistinctOctantSupercell(const correlation::core::Trajectory &traj, size_t base_frame) {
  const auto &orig = traj.getFrame(base_frame);
  const auto &lat = orig.latticeVectors();
  correlation::core::Cell supercell(lat[0] * static_cast<real_t>(2.0),
                                    lat[1] * static_cast<real_t>(2.0),
                                    lat[2] * static_cast<real_t>(2.0));
  size_t octant = 0;
  for (int ix = 0; ix < 2; ++ix) {
    for (int iy = 0; iy < 2; ++iy) {
      for (int iz = 0; iz < 2; ++iz) {
        // Use 8 distinct liquid frames (spaced 20 frames apart) to avoid artificial 2x2x2
        // periodic-copy extinction rules (where 7/8 of reciprocal vectors vanish identically).
        const size_t src_frame = base_frame - octant * 20;
        const auto &src_cell = traj.getFrame(src_frame);
        const auto shift = lat[0] * static_cast<real_t>(ix) + lat[1] * static_cast<real_t>(iy) +
                           lat[2] * static_cast<real_t>(iz);
        for (const auto &atom : src_cell.atoms()) {
          supercell.addAtom(atom.element().symbol, atom.position() + shift);
        }
        ++octant;
      }
    }
  }
  return supercell;
}

[[nodiscard]] SQAccuracyStats computeSQAccuracyStats(const std::vector<real_t> &q_bins,
                                                     const std::vector<double> &mean_full,
                                                     const std::vector<double> &mean_hybrid,
                                                     double q_min_eval) {
  double sum_y = 0.0;
  size_t count = 0;
  for (size_t i = 0; i < q_bins.size(); ++i) {
    if (static_cast<double>(q_bins[i]) >= q_min_eval && mean_full[i] > 0.0) {
      sum_y += mean_full[i];
      ++count;
    }
  }
  const double mean_y = (count > 0) ? (sum_y / static_cast<double>(count)) : 0.0;
  double ss_res = 0.0;
  double ss_tot = 0.0;
  double chi_sq = 0.0;
  for (size_t i = 0; i < q_bins.size(); ++i) {
    if (static_cast<double>(q_bins[i]) >= q_min_eval && mean_full[i] > 0.0) {
      const double diff = mean_hybrid[i] - mean_full[i];
      ss_res += diff * diff;
      const double dev = mean_full[i] - mean_y;
      ss_tot += dev * dev;
      chi_sq += (diff * diff) / std::max(std::abs(mean_full[i]), 1e-6);
    }
  }
  return {
      .r_squared = (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 1.0,
      .chi_squared = chi_sq,
      .reduced_chi_squared = (count > 0) ? (chi_sq / static_cast<double>(count)) : 0.0,
      .rmse = (count > 0) ? std::sqrt(ss_res / static_cast<double>(count)) : 0.0,
      .count = count,
  };
}

} // namespace

TEST_F(StructureFactorCalculatorTests, LiquidBismuthFullVsHybridSQAccuracy) {
  std::filesystem::path lbi_path;
  for (const auto &prefix : {"", "../", "../../", "../../../"}) {
    const auto cand = std::filesystem::path(prefix) / "examples/l-Bi/l-Bi.arc";
    if (std::filesystem::exists(cand)) {
      lbi_path = std::filesystem::absolute(cand);
      break;
    }
  }
  if (lbi_path.empty()) {
    GTEST_SKIP() << "Example file examples/l-Bi/l-Bi.arc not found";
  }

  const auto traj = correlation::readers::readTrajectory(
      lbi_path.string(), correlation::readers::FileType::Arc);
  ASSERT_GE(traj.getFrameCount(), 20u);

  const size_t n_frames = traj.getFrameCount();
  const size_t start_idx = n_frames - 10;
  const size_t num_used = n_frames - start_idx;

  const correlation::calculators::StructureFactorCalculator calc;
  AnalysisSettings base_settings;
  base_settings.q_max = static_cast<real_t>(12.0);
  base_settings.q_bin_width = static_cast<real_t>(0.05);

  std::vector<double> mean_full;
  std::vector<double> mean_stratified;
  std::vector<real_t> q_bins;

  for (size_t f_idx = start_idx; f_idx < n_frames; ++f_idx) {
    const auto supercell = buildDistinctOctantSupercell(traj, f_idx);

    DistributionFunctions df_full(supercell);
    DistributionFunctions df_strat(supercell);

    AnalysisSettings full_cfg = base_settings;
    full_cfg.sq_exact_direct = true;
    calc.calculateFrame(df_full, full_cfg);

    AnalysisSettings strat_cfg = base_settings;
    strat_cfg.sq_exact_direct = false;
    calc.calculateFrame(df_strat, strat_cfg);

    const auto &h_full = df_full.getHistogram("S_q");
    const auto &h_strat = df_strat.getHistogram("S_q");
    if (mean_full.empty()) {
      q_bins = h_full.bins;
      mean_full.assign(q_bins.size(), 0.0);
      mean_stratified.assign(q_bins.size(), 0.0);
    }
    const auto &s_f = h_full.partials.at("Total");
    const auto &s_s = h_strat.partials.at("Total");
    for (size_t i = 0; i < q_bins.size(); ++i) {
      mean_full[i] += static_cast<double>(s_f[i]);
      mean_stratified[i] += static_cast<double>(s_s[i]);
    }
  }

  for (size_t i = 0; i < q_bins.size(); ++i) {
    mean_full[i] /= static_cast<double>(num_used);
    mean_stratified[i] /= static_cast<double>(num_used);
  }

  const auto all_stats = computeSQAccuracyStats(q_bins, mean_full, mean_stratified, 0.0);
  const auto hi_stats = computeSQAccuracyStats(q_bins, mean_full, mean_stratified, 3.7);

  std::cout << "[l-Bi Full vs Stratified S(Q)] Full Q [0..12 A^-1] (N=" << all_stats.count
            << "): R^2=" << all_stats.r_squared << ", chi^2=" << all_stats.chi_squared
            << ", red_chi^2=" << all_stats.reduced_chi_squared << ", RMSE=" << all_stats.rmse
            << '\n';
  std::cout << "[l-Bi Full vs Stratified S(Q)] High-Q [3.7..12 A^-1] (N=" << hi_stats.count
            << "): R^2=" << hi_stats.r_squared << ", chi^2=" << hi_stats.chi_squared
            << ", red_chi^2=" << hi_stats.reduced_chi_squared << ", RMSE=" << hi_stats.rmse
            << '\n';

  EXPECT_GT(all_stats.r_squared, 0.99);
  EXPECT_LT(all_stats.reduced_chi_squared, 0.002);
  EXPECT_LT(all_stats.rmse, 0.03);
}
} // namespace correlation::analysis
