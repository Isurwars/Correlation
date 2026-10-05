// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "analysis/TrajectoryAnalyzer.hpp"
#include "calculators/spatial/RDFCalculator.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include "../../CrystalTestHelper.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <iterator>
#include <numbers>
#include <vector>

namespace correlation::analysis {

namespace {
// Test fixture for DistributionFunctions tests.
class RDFCalculatorTests : public ::testing::Test {
protected:
  void SetUp() override {
    // A simple cubic cell containing two atoms
    cell_ = correlation::core::Cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
    cell_.addAtom("Ar", {5.0, 5.0, 5.0});
    cell_.addAtom("Ar", {6.5, 5.0, 5.0}); // Distance 1.5
  }

  void updateTrajectory() {
    trajectory_ = correlation::core::Trajectory();
    trajectory_.addFrame(cell_);
    trajectory_.precomputeBondCutoffs();
  }

  void updateTrajectory(const correlation::core::Cell &cell) {
    trajectory_ = correlation::core::Trajectory();
    trajectory_.addFrame(cell);
    trajectory_.precomputeBondCutoffs();
  }

public:
  correlation::core::Cell cell_;
  correlation::core::Trajectory trajectory_;
};

[[nodiscard]] std::pair<real_t, real_t> findPeakInRange(const std::vector<real_t> &bins,
                                                        const std::vector<real_t> &values,
                                                        real_t r_min, real_t r_max) {
  real_t max_val = 0;
  real_t peak_r = 0;
  for (size_t i = 0; i < bins.size() && i < values.size(); ++i) {
    if (bins[i] >= r_min && bins[i] <= r_max && values[i] > max_val) {
      max_val = values[i];
      peak_r = bins[i];
    }
  }
  return {peak_r, max_val};
}

void verifyAshcroftSums(const Histogram &g_r, const Histogram &g_r_total) {
  const auto &g_ar_ar = g_r.partials.at("Ar-Ar");
  const auto &g_xe_xe = g_r.partials.at("Xe-Xe");
  const auto &g_ar_xe = g_r.partials.at("Ar-Xe");
  const auto &g_total = g_r.partials.at("Total");

  const auto &g_tot_ar_ar = g_r_total.partials.at("Ar-Ar");
  const auto &g_tot_xe_xe = g_r_total.partials.at("Xe-Xe");
  const auto &g_tot_ar_xe = g_r_total.partials.at("Ar-Xe");
  const auto &g_tot_total = g_r_total.partials.at("Total");

  ASSERT_EQ(g_total.size(), g_ar_ar.size());
  for (size_t i = 0; i < g_total.size(); ++i) {
    double const sum_g_partials = g_ar_ar[i] + g_xe_xe[i] + g_ar_xe[i];
    EXPECT_NEAR(g_total[i], sum_g_partials, 1e-6);

    double const sum_g_tot_partials = g_tot_ar_ar[i] + g_tot_xe_xe[i] + g_tot_ar_xe[i];
    EXPECT_NEAR(g_tot_total[i], sum_g_tot_partials, 1e-6);
  }
}

void verifyHrHistogram(const DistributionFunctions &dists) {
  const auto &h_r = dists.getHistogram("H_r");
  EXPECT_EQ(h_r.y_unit, "counts");
  EXPECT_EQ(h_r.title, "H(r) — Distance Histogram");
  ASSERT_TRUE(h_r.partials.contains("Ar-Ar"));
}

void verifyGrHistograms(const DistributionFunctions &dists) {
  const auto &g_unw = dists.getHistogram("g_r_unweighted");
  EXPECT_EQ(g_unw.title, "g(r) — Unweighted Radial Distribution Function");
  EXPECT_EQ(g_unw.y_unit, "");
  ASSERT_TRUE(g_unw.partials.contains("Ar-Ar"));

  const auto &g_r = dists.getHistogram("g_r");
  EXPECT_EQ(g_r.y_unit, "");

  const auto &g_r_tot = dists.getHistogram("G_r");
  EXPECT_EQ(g_r_tot.y_unit, "Å⁻²");
}
} // namespace

TEST_F(RDFCalculatorTests, DefaultConstructorWorks) {
  updateTrajectory();
  ASSERT_NO_THROW(DistributionFunctions const dists(cell_, 5.0, trajectory_.getBondCutoffsSQ()));
}

TEST_F(RDFCalculatorTests, MoveConstructorWorks) {
  updateTrajectory();
  DistributionFunctions df_source(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  df_source.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });

  DistributionFunctions const df_dest(std::move(df_source));

  EXPECT_NO_THROW(df_dest.getHistogram("g_r"));
  EXPECT_EQ(df_dest.cell().atomCount(), 2);
}

TEST_F(RDFCalculatorTests, MoveAssignmentWorks) {
  updateTrajectory();
  DistributionFunctions df_source(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  df_source.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });

  DistributionFunctions df_dest(cell_);
  df_dest = std::move(df_source);

  EXPECT_NO_THROW(df_dest.getHistogram("g_r"));
  EXPECT_EQ(df_dest.cell().atomCount(), 2);
}

TEST_F(RDFCalculatorTests, AccessorsWork) {
  updateTrajectory();
  DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // cell()
  EXPECT_EQ(dists.cell().atomCount(), 2);

  // getAvailableHistograms() - initially empty or minimal
  auto hist_names = dists.getAvailableHistograms();
  EXPECT_TRUE(hist_names.empty());

  // Calculate something
  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });
  hist_names = dists.getAvailableHistograms();
  EXPECT_FALSE(hist_names.empty());
  EXPECT_NE(std::ranges::find(hist_names, "g_r"), hist_names.end());

  // getHistogram()
  EXPECT_NO_THROW(dists.getHistogram("g_r"));
  EXPECT_THROW(dists.getHistogram("NonExistent"), std::out_of_range);

  // getAllHistograms()
  const auto &all_hists = dists.getAllHistograms();
  EXPECT_EQ(all_hists.size(), 5);
  EXPECT_TRUE(all_hists.count("g_r"));
  EXPECT_TRUE(all_hists.count("g_r_unweighted"));
  EXPECT_TRUE(all_hists.count("H_r"));
  EXPECT_TRUE(all_hists.count("J_r"));
  EXPECT_TRUE(all_hists.count("G_r"));
}

TEST_F(RDFCalculatorTests, CalculateRDF) {
  updateTrajectory();
  DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // Invalid parameters
  EXPECT_THROW(dists.calculateRDF({
                   .r_max = 5.0,
                   .r_bin_width = 0.0,
               }),
               std::invalid_argument); // Zero bin width
  EXPECT_THROW(dists.calculateRDF({
                   .r_max = 0.0,
                   .r_bin_width = 0.1,
               }),
               std::invalid_argument); // Zero r_max
  EXPECT_THROW(dists.calculateRDF({
                   .r_max = 51.0,
                   .r_bin_width = 0.1,
               }),
               std::invalid_argument); // r_max > max_cutoff_radius (50.0)
  EXPECT_THROW(dists.calculateRDF({
                   .r_max = 16.0,
                   .r_bin_width = 0.1,
                   .max_radius = 15.0,
               }),
               std::invalid_argument); // Custom max_radius exceeded
  EXPECT_THROW(static_cast<void>(correlation::calculators::RDFCalculator::calculate(cell_, nullptr,
                                                                                    {}, 51.0, 0.1)),
               std::invalid_argument); // Direct calculate exceeding max_cutoff_radius

  // Valid calculation with tight bins for numerical accuracy
  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.001,
  });
  const auto &hist = dists.getHistogram("g_r");
  const auto &total = hist.partials.at("Ar-Ar");

  // High precision peak location
  auto max_it = std::ranges::max_element(total);
  size_t const peak_idx = std::distance(total.begin(), max_it);
  real_t const peak_r = hist.bins[peak_idx];

  // Bin size is 0.001. Bin containing 1.50 is index 1500 (center 1.5005) or
  // 1499 (1.4995)
  EXPECT_NEAR(peak_r, 1.5, 0.001);
}

TEST_F(RDFCalculatorTests, CalculateCoordinationNumber) {
  // Use a setup where we know neighbors exactly
  correlation::core::Cell cn_cell({10, 10, 10, 90, 90, 90});
  cn_cell.addAtom("Si", {5.0, 5.0, 5.0});
  cn_cell.addAtom("O", {6.6, 5.0, 5.0}); // 1.6 dist
  cn_cell.addAtom("O", {3.4, 5.0, 5.0}); // 1.6 dist
  // Si has 2 O neighbors at 1.6.

  updateTrajectory(cn_cell);
  DistributionFunctions dists(cn_cell, 2.0, trajectory_.getBondCutoffsSQ());

  dists.calculateCoordinationNumber();

  const auto &hist = dists.getHistogram("CN");
  const auto &sio_cn = hist.partials.at("Si-O");

  // Si has 2 O neighbors. So bin 2 should be 1.
  // Ensure sufficient size
  ASSERT_GT(sio_cn.size(), 2);
  EXPECT_EQ(sio_cn[2], 1);

  const auto &osi_cn = hist.partials.at("O-Si");
  // Each O has 1 Si neighbor. There are 2 Os. So bin 1 should be 2.
  ASSERT_GT(osi_cn.size(), 1);
  EXPECT_EQ(osi_cn[1], 2);
}

TEST_F(RDFCalculatorTests, Smoothing) {
  updateTrajectory();
  DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });

  // Checks "smooth" single
  ASSERT_NO_THROW(dists.smooth("g_r", 0.2));
  const auto &hist = dists.getHistogram("g_r");
  EXPECT_FALSE(hist.smoothed_partials.empty());

  // Check "smoothAll"
  // Add another histogram
  dists.calculateCoordinationNumber();
  dists.smoothAll(0.2);
  const auto &cn_hist = dists.getHistogram("CN");
  EXPECT_FALSE(cn_hist.smoothed_partials.empty());
}

TEST_F(RDFCalculatorTests, SetStructureAnalyzer) {
  updateTrajectory();
  // Create an external analyzer
  StructureAnalyzer const analyzer(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  DistributionFunctions dists(cell_, 0.0, {});
  // Should depend on analyzer for RDF
  dists.setStructureAnalyzer(&analyzer);

  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });
  EXPECT_NO_THROW(dists.getHistogram("g_r"));
  EXPECT_FALSE(dists.getHistogram("g_r").partials.empty());
}

TEST_F(RDFCalculatorTests, AddAndScale) {
  updateTrajectory();
  DistributionFunctions df1(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  df1.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });

  DistributionFunctions df2(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  df2.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });

  // Add
  df1.add(df2);
  const auto &h_1 = df1.getHistogram("g_r").partials.at("Ar-Ar");

  // Peak should be doubled roughly (since they are identical)
  // Actually add() sums the bins.
  // If both have 1 count at peak, sum is 2.

  // Scale
  df1.scale(0.5);
  const auto &h1_scaled = df1.getHistogram("g_r").partials.at("Ar-Ar");

  // Should be back to original magnitude
  // We check peak value
  auto max_it = std::ranges::max_element(h1_scaled);
  real_t const peak = *max_it;

  // Single frame RDF peak value depends on volume and density, but it's
  // consistent. Let's compare with a fresh one
  DistributionFunctions df_ref(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  df_ref.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });
  real_t const ref_peak =
      *std::ranges::max_element(df_ref.getHistogram("g_r").partials.at("Ar-Ar"));

  EXPECT_NEAR(peak, ref_peak, 1e-4);
}

TEST_F(RDFCalculatorTests, ComputeMean) {
  updateTrajectory();
  TrajectoryAnalyzer const analyzer(trajectory_, 5.0, trajectory_.getBondCutoffsSQ());

  AnalysisSettings settings;
  settings.r_max = 5.0;
  settings.r_bin_width = 0.1;
  settings.smoothing = false;
  settings.active_calculators["RDF"] = true;

  auto df_mean = DistributionFunctions::computeMean(trajectory_, analyzer, 0, settings);
  ASSERT_TRUE(df_mean != nullptr);
  EXPECT_NO_THROW(df_mean->getHistogram("g_r"));
}

TEST_F(RDFCalculatorTests, HandlesMissingPartialInAdd) {
  // Build cell 1: pure Ar
  correlation::core::Cell c_1({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  c_1.addAtom("Ar", {0.0, 0.0, 0.0});
  c_1.addAtom("Ar", {2.0, 0.0, 0.0});

  // Build cell 2: Ar and Xe
  correlation::core::Cell c_2({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  c_2.addAtom("Ar", {0.0, 0.0, 0.0});
  c_2.addAtom("Xe", {2.5, 0.0, 0.0});

  correlation::core::Trajectory t_1;
  t_1.addFrame(c_1);
  t_1.precomputeBondCutoffs();

  correlation::core::Trajectory t_2;
  t_2.addFrame(c_2);
  t_2.precomputeBondCutoffs();

  DistributionFunctions df1(c_1, 5.0, t_1.getBondCutoffsSQ());
  DistributionFunctions df2(c_2, 5.0, t_2.getBondCutoffsSQ());

  df1.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });
  df2.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });

  // df1 initially only has Ar-Ar (plus Total)
  const auto &g_r_df1_before = df1.getHistogram("g_r");
  EXPECT_TRUE(g_r_df1_before.partials.count("Ar-Ar"));
  EXPECT_FALSE(g_r_df1_before.partials.count("Ar-Xe"));

  // Act
  df1.add(df2);

  // df1 should now have acquired Ar-Xe and Xe-Xe partials
  const auto &g_r_df1_after = df1.getHistogram("g_r");
  EXPECT_TRUE(g_r_df1_after.partials.count("Ar-Ar"));
  EXPECT_TRUE(g_r_df1_after.partials.count("Ar-Xe"));
  EXPECT_TRUE(g_r_df1_after.partials.count("Xe-Xe"));
}

TEST_F(RDFCalculatorTests, VerifyAshcroftWeightsAreCorrect) {
  // Build cell with 3 Ar atoms and 1 Xe atom (total 4 atoms)
  // Concentration: x_Ar = 0.75, x_Xe = 0.25
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {1.0, 1.0, 1.0});
  cell.addAtom("Ar", {2.0, 2.0, 2.0});
  cell.addAtom("Ar", {3.0, 3.0, 3.0});
  cell.addAtom("Xe", {4.0, 4.0, 4.0});

  correlation::core::Trajectory traj;
  traj.addFrame(cell);
  traj.precomputeBondCutoffs();

  DistributionFunctions dists(cell, 5.0, traj.getBondCutoffsSQ());

  // 1. Verify calculated Ashcroft weights
  const auto &weights = dists.getAshcroftWeights();
  double const expected_w_ar_ar = 0.75 * 0.75;       // 0.5625
  double const expected_w_xe_xe = 0.25 * 0.25;       // 0.0625
  double const expected_w_ar_xe = 2.0 * 0.75 * 0.25; // 0.375 (doubled!)

  EXPECT_NEAR(weights.at("Ar-Ar"), expected_w_ar_ar, 1e-6);
  EXPECT_NEAR(weights.at("Xe-Xe"), expected_w_xe_xe, 1e-6);
  EXPECT_NEAR(weights.at("Ar-Xe"), expected_w_ar_xe, 1e-6);

  // 2. Verify RDF total is the sum of weighted partials
  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });
  const auto &g_r = dists.getHistogram("g_r");
  const auto &g_r_total = dists.getHistogram("G_r");
  verifyAshcroftSums(g_r, g_r_total);
}

TEST_F(RDFCalculatorTests, FCC_Copper_RDF) {
  real_t const lat_a = 3.615;
  auto cell = correlation::testing::crystals::createFCCCell(lat_a, "Cu", 2, 2, 2);
  updateTrajectory(cell);

  DistributionFunctions dists(cell, 5.0, trajectory_.getBondCutoffsSQ());
  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.01,
  });

  const auto &hist = dists.getHistogram("g_r");
  const auto &cu_cu = hist.partials.at("Cu-Cu");

  const auto [peak_r1, max_val1] = findPeakInRange(hist.bins, cu_cu, 2.0, 3.0);
  EXPECT_NEAR(peak_r1, lat_a / std::numbers::sqrt2, 0.02);
  EXPECT_GT(max_val1, 1.0);
}

TEST_F(RDFCalculatorTests, BCC_Iron_RDF) {
  real_t const lat_a = 2.866;
  auto cell = correlation::testing::crystals::createBCCCell(lat_a, "Fe", 3, 3, 3);
  updateTrajectory(cell);

  DistributionFunctions dists(cell, 4.5, trajectory_.getBondCutoffsSQ());
  dists.calculateRDF({
      .r_max = 4.5,
      .r_bin_width = 0.01,
  });

  const auto &hist = dists.getHistogram("g_r");
  const auto &fe_fe = hist.partials.at("Fe-Fe");

  const auto [peak_r1, max_val1] = findPeakInRange(hist.bins, fe_fe, 2.0, 2.7);
  EXPECT_NEAR(peak_r1, lat_a * std::sqrt(3.0) / 2.0, 0.02);
  EXPECT_GT(max_val1, 1.0);
}

TEST_F(RDFCalculatorTests, Diamond_Silicon_RDF) {
  real_t const lat_a = 5.431;
  auto cell = correlation::testing::crystals::createDiamondCell(lat_a, "Si", 2, 2, 2);
  updateTrajectory(cell);

  DistributionFunctions dists(cell, 5.0, trajectory_.getBondCutoffsSQ());
  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.01,
  });

  const auto &hist = dists.getHistogram("g_r");
  const auto &si_si = hist.partials.at("Si-Si");

  const auto [peak_r1, max_val1] = findPeakInRange(hist.bins, si_si, 2.0, 2.6);
  EXPECT_NEAR(peak_r1, lat_a * std::sqrt(3.0) / 4.0, 0.02);
  EXPECT_GT(max_val1, 1.0);
}

TEST_F(RDFCalculatorTests, NaCl_RockSalt_RDF) {
  real_t const lat_a = 5.64;
  auto cell = correlation::testing::crystals::createNaClCell(lat_a, "Na", "Cl", 2, 2, 2);
  updateTrajectory(cell);

  DistributionFunctions dists(cell, 5.0, trajectory_.getBondCutoffsSQ());
  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.01,
  });

  const auto &hist = dists.getHistogram("g_r");
  const auto &nacl = hist.partials.at("Na-Cl");
  const auto &nana = hist.partials.at("Na-Na");

  const auto [peak_nacl, max_nacl] = findPeakInRange(hist.bins, nacl, 2.5, 3.2);
  EXPECT_NEAR(peak_nacl, lat_a / 2.0, 0.02);
  EXPECT_GT(max_nacl, 1.0);

  const auto [peak_nana, max_nana] = findPeakInRange(hist.bins, nana, 3.6, 4.4);
  EXPECT_NEAR(peak_nana, lat_a / std::numbers::sqrt2, 0.02);
  EXPECT_GT(max_nana, 1.0);
}

TEST_F(RDFCalculatorTests, VerifyRawAndUnweightedHistograms) {
  updateTrajectory();
  DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.1,
  });

  ASSERT_TRUE(dists.getAllHistograms().contains("H_r"));
  ASSERT_TRUE(dists.getAllHistograms().contains("g_r_unweighted"));

  const auto &h_r = dists.getHistogram("H_r");
  real_t raw_counts_sum = 0;
  for (const real_t val : h_r.partials.at("Ar-Ar")) {
    raw_counts_sum += val;
  }
  EXPECT_GT(raw_counts_sum, 0.0);

  verifyHrHistogram(dists);
  verifyGrHistograms(dists);
}

TEST_F(RDFCalculatorTests, AddAccumulatesWithMismatchedPartialSizes) {
  updateTrajectory();
  DistributionFunctions df1(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  DistributionFunctions df2(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  Histogram h_1;
  h_1.bins = {0.1, 0.2, 0.3};
  h_1.partials["Ar-Ar"] = {1.0, 2.0, 3.0};
  h_1.compute_count = 1;

  Histogram h_2;
  h_2.bins = {0.1, 0.2, 0.3};
  h_2.partials["Ar-Ar"] = {10.0, 20.0}; // Shorter partial size
  h_2.compute_count = 1;

  df1.addHistogram("test_hist", std::move(h_1));
  df2.addHistogram("test_hist", std::move(h_2));

  // Should safely accumulate min(3, 2) = 2 elements without throwing or reading out-of-bounds
  EXPECT_NO_THROW(df1.add(df2));

  const auto &res = df1.getHistogram("test_hist");
  EXPECT_EQ(res.compute_count, 2);
  ASSERT_EQ(res.partials.at("Ar-Ar").size(), 3);
  EXPECT_DOUBLE_EQ(res.partials.at("Ar-Ar")[0], 11.0);
  EXPECT_DOUBLE_EQ(res.partials.at("Ar-Ar")[1], 22.0);
  EXPECT_DOUBLE_EQ(res.partials.at("Ar-Ar")[2], 3.0);
}

TEST_F(RDFCalculatorTests, PrimitiveCellRDFComputesSelfImages) {
  // A primitive FCC cell with 1 atom (e.g., Cu a = 3.615)
  // Distance to nearest neighbor is a / sqrt(2) ~ 2.556 A
  const auto lat_a = static_cast<real_t>(3.615);
  const real_t half_a = lat_a / static_cast<real_t>(2.0);
  const correlation::math::Vector3<real_t> v_a{static_cast<real_t>(0.0), half_a, half_a};
  const correlation::math::Vector3<real_t> v_b{half_a, static_cast<real_t>(0.0), half_a};
  const correlation::math::Vector3<real_t> v_c{half_a, half_a, static_cast<real_t>(0.0)};
  correlation::core::Cell prim_cell(v_a, v_b, v_c);
  prim_cell.addAtom("Cu", {0.0, 0.0, 0.0});

  updateTrajectory(prim_cell);
  DistributionFunctions dists(prim_cell, 5.0, trajectory_.getBondCutoffsSQ());
  dists.calculateRDF({
      .r_max = 5.0,
      .r_bin_width = 0.02,
  });

  const auto &hist = dists.getHistogram("g_r");
  ASSERT_TRUE(hist.partials.contains("Cu-Cu"));
  const auto &cu_cu = hist.partials.at("Cu-Cu");

  const auto [peak_r, max_val] = findPeakInRange(hist.bins, cu_cu, 2.0, 3.0);
  EXPECT_NEAR(peak_r, lat_a / std::numbers::sqrt2, 0.03);
  EXPECT_GT(max_val, 1.0);

  const auto &h_r = dists.getHistogram("H_r");
  ASSERT_TRUE(h_r.partials.contains("Cu-Cu"));
  real_t total_counts = 0.0;
  for (const real_t val : h_r.partials.at("Cu-Cu")) {
    total_counts += val;
  }
  EXPECT_GT(total_counts, 0.0);
}

} // namespace correlation::analysis
