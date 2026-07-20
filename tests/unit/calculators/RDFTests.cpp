// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "analysis/TrajectoryAnalyzer.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <algorithm>
#include <gtest/gtest.h>
#include <iterator>
#include <vector>

namespace correlation::analysis {

namespace {
// Test fixture for DistributionFunctions tests.
class RDFTests : public ::testing::Test {
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
} // namespace

TEST_F(RDFTests, DefaultConstructorWorks) {
  updateTrajectory();
  ASSERT_NO_THROW(DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ()));
}

TEST_F(RDFTests, MoveConstructorWorks) {
  updateTrajectory();
  DistributionFunctions dfSource(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  dfSource.calculateRDF(5.0, 0.1);

  DistributionFunctions const dfDest(std::move(dfSource));

  EXPECT_NO_THROW(dfDest.getHistogram("g_r"));
  EXPECT_EQ(dfDest.cell().atomCount(), 2);
}

TEST_F(RDFTests, MoveAssignmentWorks) {
  updateTrajectory();
  DistributionFunctions dfSource(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  dfSource.calculateRDF(5.0, 0.1);

  DistributionFunctions dfDest(cell_, 0.0, std::vector<std::vector<real_t>>{});
  dfDest = std::move(dfSource);

  EXPECT_NO_THROW(dfDest.getHistogram("g_r"));
  EXPECT_EQ(dfDest.cell().atomCount(), 2);
}


TEST_F(RDFTests, AccessorsWork) {
  updateTrajectory();
  DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // cell()
  EXPECT_EQ(dists.cell().atomCount(), 2);

  // getAvailableHistograms() - initially empty or minimal
  auto histNames = dists.getAvailableHistograms();
  EXPECT_TRUE(histNames.empty());

  // Calculate something
  dists.calculateRDF(5.0, 0.1);
  histNames = dists.getAvailableHistograms();
  EXPECT_FALSE(histNames.empty());
  EXPECT_NE(std::find(histNames.begin(), histNames.end(), "g_r"), histNames.end());

  // getHistogram()
  EXPECT_NO_THROW(dists.getHistogram("g_r"));
  EXPECT_THROW(dists.getHistogram("NonExistent"), std::out_of_range);

  // getAllHistograms()
  const auto &allHists = dists.getAllHistograms();
  EXPECT_EQ(allHists.size(), 3);
  EXPECT_TRUE(allHists.count("g_r"));
  EXPECT_TRUE(allHists.count("J_r"));
  EXPECT_TRUE(allHists.count("G_r"));
}


TEST_F(RDFTests, CalculateRDF) {
  updateTrajectory();
  DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // Invalid parameters
  EXPECT_THROW(dists.calculateRDF(5.0, 0.0),
               std::invalid_argument);                               // Zero bin width
  EXPECT_THROW(dists.calculateRDF(0.0, 0.1), std::invalid_argument); // Zero r_max

  // Valid calculation with tight bins for numerical accuracy
  dists.calculateRDF(5.0, 0.001);
  const auto &hist = dists.getHistogram("g_r");
  const auto &total = hist.partials.at("Ar-Ar");

  // High precision peak location
  auto max_it = std::max_element(total.begin(), total.end());
  size_t const peak_idx = std::distance(total.begin(), max_it);
  double const peak_r = hist.bins[peak_idx];

  // Bin size is 0.001. Bin containing 1.50 is index 1500 (center 1.5005) or
  // 1499 (1.4995)
  EXPECT_NEAR(peak_r, 1.5, 0.001);
}

TEST_F(RDFTests, CalculateCoordinationNumber) {
  // Use a setup where we know neighbors exactly
  correlation::core::Cell cnCall({10, 10, 10, 90, 90, 90});
  cnCall.addAtom("Si", {5, 5, 5});
  cnCall.addAtom("O", {6, 5, 5}); // 1.0 dist
  cnCall.addAtom("O", {4, 5, 5}); // 1.0 dist
  // Si has 2 O neighbors at 1.0.

  updateTrajectory(cnCall);
  DistributionFunctions dists(cnCall, 2.0, trajectory_.getBondCutoffsSQ());

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

TEST_F(RDFTests, Smoothing) {
  updateTrajectory();
  DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  dists.calculateRDF(5.0, 0.1);

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

TEST_F(RDFTests, SetStructureAnalyzer) {
  updateTrajectory();
  // Create an external analyzer
  StructureAnalyzer const analyzer(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  DistributionFunctions dists(cell_, 0.0, {});
  // Should depend on analyzer for RDF
  dists.setStructureAnalyzer(&analyzer);

  dists.calculateRDF(5.0, 0.1);
  EXPECT_NO_THROW(dists.getHistogram("g_r"));
  EXPECT_FALSE(dists.getHistogram("g_r").partials.empty());
}


TEST_F(RDFTests, AddAndScale) {
  updateTrajectory();
  DistributionFunctions df1(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  df1.calculateRDF(5.0, 0.1);

  DistributionFunctions df2(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  df2.calculateRDF(5.0, 0.1);

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
  double const peak = *max_it;

  // Single frame RDF peak value depends on volume and density, but it's
  // consistent. Let's compare with a fresh one
  DistributionFunctions dfRef(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  dfRef.calculateRDF(5.0, 0.1);
  double const refPeak = *std::max_element(dfRef.getHistogram("g_r").partials.at("Ar-Ar").begin(),
                                           dfRef.getHistogram("g_r").partials.at("Ar-Ar").end());

  EXPECT_NEAR(peak, refPeak, 1e-6);
}

TEST_F(RDFTests, ComputeMean) {
  updateTrajectory();
  TrajectoryAnalyzer const analyzer(trajectory_, 5.0, trajectory_.getBondCutoffsSQ());

  AnalysisSettings settings;
  settings.r_max = 5.0;
  settings.r_bin_width = 0.1;
  settings.smoothing = false;

  auto dfMean = DistributionFunctions::computeMean(trajectory_, analyzer, 0, settings);
  ASSERT_TRUE(dfMean != nullptr);
  EXPECT_NO_THROW(dfMean->getHistogram("g_r"));
}

TEST_F(RDFTests, HandlesMissingPartialInAdd) {
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

  df1.calculateRDF(5.0, 0.1);
  df2.calculateRDF(5.0, 0.1);

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

TEST_F(RDFTests, VerifyAshcroftWeightsAreCorrect) {
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
  double const expected_w_ArAr = 0.75 * 0.75;       // 0.5625
  double const expected_w_XeXe = 0.25 * 0.25;       // 0.0625
  double const expected_w_ArXe = 2.0 * 0.75 * 0.25; // 0.375 (doubled!)

  EXPECT_NEAR(weights.at("Ar-Ar"), expected_w_ArAr, 1e-6);
  EXPECT_NEAR(weights.at("Xe-Xe"), expected_w_XeXe, 1e-6);
  EXPECT_NEAR(weights.at("Ar-Xe"), expected_w_ArXe, 1e-6);

  // 2. Verify RDF total is the sum of weighted partials
  dists.calculateRDF(5.0, 0.1);
  const auto &g_r = dists.getHistogram("g_r");
  const auto &G_r = dists.getHistogram("G_r");

  const auto &g_ArAr = g_r.partials.at("Ar-Ar");
  const auto &g_XeXe = g_r.partials.at("Xe-Xe");
  const auto &g_ArXe = g_r.partials.at("Ar-Xe");
  const auto &g_Total = g_r.partials.at("Total");

  const auto &G_ArAr = G_r.partials.at("Ar-Ar");
  const auto &G_XeXe = G_r.partials.at("Xe-Xe");
  const auto &G_ArXe = G_r.partials.at("Ar-Xe");
  const auto &G_Total = G_r.partials.at("Total");

  ASSERT_EQ(g_Total.size(), g_ArAr.size());
  for (size_t i = 0; i < g_Total.size(); ++i) {
    double const sum_g_partials = g_ArAr[i] + g_XeXe[i] + g_ArXe[i];
    EXPECT_NEAR(g_Total[i], sum_g_partials, 1e-6);

    double const sum_G_partials = G_ArAr[i] + G_XeXe[i] + G_ArXe[i];
    EXPECT_NEAR(G_Total[i], sum_G_partials, 1e-6);
  }
}
} // namespace correlation::analysis
