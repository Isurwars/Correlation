// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <algorithm>
#include <gtest/gtest.h>
#include <iterator>
#include <vector>

#include "../include/Cell.hpp"
#include "../include/DistributionFunctions.hpp"
#include "../include/PhysicalData.hpp"
#include "../include/Trajectory.hpp"
#include "../include/TrajectoryAnalyzer.hpp"

// Test fixture for DistributionFunctions tests.
class Test06_DistributionFunctions : public ::testing::Test {
protected:
  void SetUp() override {
    // A simple cubic cell containing two atoms
    cell_ = Cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
    cell_.addAtom("Ar", {5.0, 5.0, 5.0});
    cell_.addAtom("Ar", {6.5, 5.0, 5.0}); // Distance 1.5
  }

  void updateTrajectory() {
    trajectory_ = Trajectory();
    trajectory_.addFrame(cell_);
    trajectory_.precomputeBondCutoffs();
  }

  void updateTrajectory(const Cell &c) {
    trajectory_ = Trajectory();
    trajectory_.addFrame(c);
    trajectory_.precomputeBondCutoffs();
  }

  Cell cell_{};
  Trajectory trajectory_;
};

//---------------------------------------------------------------------------//
//----------------------------- Constructors --------------------------------//
//---------------------------------------------------------------------------//

TEST_F(Test06_DistributionFunctions, DefaultConstructorWorks) {
  updateTrajectory();
  ASSERT_NO_THROW(
      DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffsSQ()));
}

TEST_F(Test06_DistributionFunctions, MoveConstructorWorks) {
  updateTrajectory();
  DistributionFunctions dfSource(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  dfSource.calculateRDF(5.0, 0.1);

  DistributionFunctions dfDest(std::move(dfSource));

  EXPECT_NO_THROW(dfDest.getHistogram("g(r)"));
  EXPECT_EQ(dfDest.cell().atomCount(), 2);
}

TEST_F(Test06_DistributionFunctions, MoveAssignmentWorks) {
  updateTrajectory();
  DistributionFunctions dfSource(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  dfSource.calculateRDF(5.0, 0.1);

  DistributionFunctions dfDest(cell_, 0.0, std::vector<std::vector<double>>{});
  dfDest = std::move(dfSource);

  EXPECT_NO_THROW(dfDest.getHistogram("g(r)"));
  EXPECT_EQ(dfDest.cell().atomCount(), 2);
}

//---------------------------------------------------------------------------//
//------------------------------- Accessors ---------------------------------//
//---------------------------------------------------------------------------//

TEST_F(Test06_DistributionFunctions, AccessorsWork) {
  updateTrajectory();
  DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // cell()
  EXPECT_EQ(df.cell().atomCount(), 2);

  // getAvailableHistograms() - initially empty or minimal
  auto histNames = df.getAvailableHistograms();
  EXPECT_TRUE(histNames.empty());

  // Calculate something
  df.calculateRDF(5.0, 0.1);
  histNames = df.getAvailableHistograms();
  EXPECT_FALSE(histNames.empty());
  EXPECT_NE(std::find(histNames.begin(), histNames.end(), "g(r)"),
            histNames.end());

  // getHistogram()
  EXPECT_NO_THROW(df.getHistogram("g(r)"));
  EXPECT_THROW(df.getHistogram("NonExistent"), std::out_of_range);

  // getAllHistograms()
  const auto &allHists = df.getAllHistograms();
  EXPECT_EQ(allHists.size(), 3);
  EXPECT_TRUE(allHists.count("g(r)"));
  EXPECT_TRUE(allHists.count("J(r)"));
  EXPECT_TRUE(allHists.count("G(r)"));
}

//---------------------------------------------------------------------------//
//--------------------------- Calculation Methods ---------------------------//
//---------------------------------------------------------------------------//

TEST_F(Test06_DistributionFunctions, CalculateRDF) {
  updateTrajectory();
  DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // Invalid parameters
  EXPECT_THROW(df.calculateRDF(5.0, 0.0),
               std::invalid_argument); // Zero bin width
  EXPECT_THROW(df.calculateRDF(0.0, 0.1), std::invalid_argument); // Zero r_max

  // Valid calculation
  df.calculateRDF(5.0, 0.1);
  const auto &hist = df.getHistogram("g(r)");
  const auto &total = hist.partials.at("Ar-Ar");

  // Peak at 1.5 (+- bin width)
  auto max_it = std::max_element(total.begin(), total.end());
  size_t peak_idx = std::distance(total.begin(), max_it);
  double peak_r = hist.bins[peak_idx];
  EXPECT_NEAR(peak_r, 1.5, 0.1);
}

TEST_F(Test06_DistributionFunctions, CalculateCoordinationNumber) {
  // Use a setup where we know neighbors exactly
  Cell cnCall({10, 10, 10, 90, 90, 90});
  cnCall.addAtom("Si", {5, 5, 5});
  cnCall.addAtom("O", {6, 5, 5}); // 1.0 dist
  cnCall.addAtom("O", {4, 5, 5}); // 1.0 dist
  // Si has 2 O neighbors at 1.0.

  updateTrajectory(cnCall);
  DistributionFunctions df(cnCall, 2.0, trajectory_.getBondCutoffsSQ());

  df.calculateCoordinationNumber();

  const auto &hist = df.getHistogram("CN");
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

TEST_F(Test06_DistributionFunctions, CalculatePAD) {
  // Water molecule angle 104.5ish
  Cell water({10, 10, 10, 90, 90, 90});
  water.addAtom("O", {5, 5, 5});
  water.addAtom("H", {6, 5, 5});
  double angRad = 104.5 * constants::deg2rad;
  water.addAtom("H", {5 + std::cos(angRad), 5 + std::sin(angRad), 5.0});

  updateTrajectory(water);
  DistributionFunctions df(water, 2.0, trajectory_.getBondCutoffsSQ());

  df.calculatePAD(1.0);
  const auto &hist = df.getHistogram("f(theta)");
  const auto &hoh = hist.partials.at("H-O-H");

  auto max_it = std::max_element(hoh.begin(), hoh.end());
  size_t idx = std::distance(hoh.begin(), max_it);
  double angle = hist.bins[idx];
  EXPECT_NEAR(angle, 104.5, 2.0);
}

TEST_F(Test06_DistributionFunctions, CalculateSQ) {
  updateTrajectory();
  DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // Need RDF first
  df.calculateRDF(5.0, 0.05);
  df.calculateSQ(10.0, 0.1, 5.0);

  EXPECT_NO_THROW(df.getHistogram("S(Q)"));
  const auto &hist = df.getHistogram("S(Q)");
  EXPECT_FALSE(hist.bins.empty());

  // Check peak position for 1.5A distance -> Q = 2pi/r ~ 4.18
  const auto &total_sq = hist.partials.at("Total");

  // For a dimer, S(Q) oscillates around 1. S(Q) = 1 + sin(Qr)/(Qr).
  // At Q = 4.18 (2pi/r), Qr = 2pi, sin(2pi)=0. So S(Q) should be approx 1.
  // The previous test expected a peak (large value), which is wrong for a
  // 2-atom system.

  // Search for value at Q ~ 4.18
  // Find bin closest to 4.18
  double target_Q = 4.18;
  double min_diff = 1000.0;
  double sq_val = 0.0;

  for (size_t i = 0; i < total_sq.size(); ++i) {
    if (std::abs(hist.bins[i] - target_Q) < min_diff) {
      min_diff = std::abs(hist.bins[i] - target_Q);
      sq_val = total_sq[i];
    }
  }

  // S(Q) should be near 1.0
  EXPECT_NEAR(sq_val, 1.0, 0.2);
}

TEST_F(Test06_DistributionFunctions, CalculateXRD) {
  updateTrajectory();
  DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // Needs RDF first
  df.calculateRDF(5.0, 0.05);

  // Calling calculateXRD
  df.calculateXRD(1.5406, 5.0, 90.0, 0.5);

  EXPECT_NO_THROW(df.getHistogram("XRD"));
  const auto &hist = df.getHistogram("XRD");
  EXPECT_FALSE(hist.bins.empty());

  // Check bins range
  EXPECT_GE(hist.bins.front(), 5.0);
  EXPECT_LE(hist.bins.back(), 90.0);
}

TEST_F(Test06_DistributionFunctions, CalculateVACF_and_VDOS) {
  Cell c({10, 10, 10, 90, 90, 90});
  c.addAtom("Ar", {0, 0, 0});
  Trajectory t;
  t.addFrame(c);
  t.addFrame(c); // Static
  t.calculateVelocities();

  DistributionFunctions df(c, 0.0, {{0.0}});

  df.calculateVACF(t, 1);
  EXPECT_NO_THROW(df.getHistogram("VACF"));
  const auto &vacf = df.getHistogram("VACF").partials.at("Total");
  // Should be 1.0 normalized (if constant 0 velocity? Wait 0 velocity -> ??)
  // If static, position constant -> velocity 0. Correlation of 0 with 0 is 0.
  // Let's give it velocity.

  Trajectory tMoving;
  Cell c1 = c;

  Cell c2({10, 10, 10, 90, 90, 90});
  // Atom 1 moves +1.0 in x
  c2.addAtom("Ar", {1.0, 0.0, 0.0});
  // Atom 2 moves -1.0 in x (balancing COM)
  c2.addAtom("Ar", {-1.0, 0.0, 0.0});

  Cell c3({10, 10, 10, 90, 90, 90});
  // Atom 1 moves to +2.0
  c3.addAtom("Ar", {2.0, 0.0, 0.0});
  // Atom 2 moves to -2.0
  c3.addAtom("Ar", {-2.0, 0.0, 0.0});

  // Need to update c1 (frame 0) to have 2 atoms at 0
  c1 = Cell({10, 10, 10, 90, 90, 90});
  c1.addAtom("Ar", {0.0, 0.0, 0.0});
  c1.addAtom("Ar", {0.0, 0.0, 0.0});

  tMoving.addFrame(c1);
  tMoving.addFrame(c2);
  tMoving.addFrame(c3);
  tMoving.setTimeStep(1.0);
  tMoving.calculateVelocities();

  df.calculateVACF(tMoving, 1);
  const auto &vacf2 = df.getHistogram("Normalized VACF").partials.at("Total");
  EXPECT_NEAR(vacf2[0], 1.0, 1e-6); // t=0
  EXPECT_NEAR(vacf2[1], 1.0, 1e-6); // t=1, const velocity

  // VDOS
  df.calculateVDOS();
  EXPECT_NO_THROW(df.getHistogram("VDOS"));
  const auto &vdos_hist = df.getHistogram("VDOS");
  EXPECT_TRUE(vdos_hist.partials.count("Frequency (cm-1)"));
}

TEST_F(Test06_DistributionFunctions, Smoothing) {
  updateTrajectory();
  DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  df.calculateRDF(5.0, 0.1);

  // Checks "smooth" single
  ASSERT_NO_THROW(df.smooth("g(r)", 0.2));
  const auto &hist = df.getHistogram("g(r)");
  EXPECT_FALSE(hist.smoothed_partials.empty());

  // Check "smoothAll"
  // Add another histogram
  df.calculateCoordinationNumber();
  df.smoothAll(0.2);
  const auto &cn_hist = df.getHistogram("CN");
  EXPECT_FALSE(cn_hist.smoothed_partials.empty());
}

TEST_F(Test06_DistributionFunctions, SetStructureAnalyzer) {
  updateTrajectory();
  // Create an external analyzer
  StructureAnalyzer analyzer(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  DistributionFunctions df(cell_, 0.0, {});
  // Should depend on analyzer for RDF
  df.setStructureAnalyzer(&analyzer);

  df.calculateRDF(5.0, 0.1);
  EXPECT_NO_THROW(df.getHistogram("g(r)"));
  EXPECT_FALSE(df.getHistogram("g(r)").partials.empty());
}

//---------------------------------------------------------------------------//
//----------------------------- Accumulation --------------------------------//
//---------------------------------------------------------------------------//

TEST_F(Test06_DistributionFunctions, AddAndScale) {
  updateTrajectory();
  DistributionFunctions df1(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  df1.calculateRDF(5.0, 0.1);

  DistributionFunctions df2(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  df2.calculateRDF(5.0, 0.1);

  // Add
  df1.add(df2);
  const auto &h1 = df1.getHistogram("g(r)").partials.at("Ar-Ar");

  // Peak should be doubled roughly (since they are identical)
  // Actually add() sums the bins.
  // If both have 1 count at peak, sum is 2.

  // Scale
  df1.scale(0.5);
  const auto &h1_scaled = df1.getHistogram("g(r)").partials.at("Ar-Ar");

  // Should be back to original magnitude
  // We check peak value
  auto max_it = std::max_element(h1_scaled.begin(), h1_scaled.end());
  double peak = *max_it;

  // Single frame RDF peak value depends on volume and density, but it's
  // consistent. Let's compare with a fresh one
  DistributionFunctions dfRef(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  dfRef.calculateRDF(5.0, 0.1);
  double refPeak =
      *std::max_element(dfRef.getHistogram("g(r)").partials.at("Ar-Ar").begin(),
                        dfRef.getHistogram("g(r)").partials.at("Ar-Ar").end());

  EXPECT_NEAR(peak, refPeak, 1e-6);
}

TEST_F(Test06_DistributionFunctions, ComputeMean) {
  updateTrajectory();
  TrajectoryAnalyzer ta(trajectory_, 5.0, trajectory_.getBondCutoffsSQ());

  AnalysisSettings settings;
  settings.r_max = 5.0;
  settings.r_bin_width = 0.1;
  settings.smoothing = false;

  auto dfMean =
      DistributionFunctions::computeMean(trajectory_, ta, 0, settings);
  ASSERT_TRUE(dfMean != nullptr);
  EXPECT_NO_THROW(dfMean->getHistogram("g(r)"));
}

TEST_F(Test06_DistributionFunctions, HandlesMissingPartialInAdd) {
  // Edge case: adding DFs with different partials (e.g. from different systems
  // or logic) Though usually DFs added should be compatible. If df2 has a key
  // that df1 doesn't, df1 should acquire it.

  updateTrajectory();
  DistributionFunctions df1(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  DistributionFunctions df2(cell_, 5.0, trajectory_.getBondCutoffsSQ());

  // Hack: manually inject a histogram if possible, or use different calculation
  // Hard to simulate without different system composition.
  // Let's try different calculation. df1 has RDF, df2 has CN.
  // Add df2 to df1. df1 should have CN now?
  // add() implementation usually iterates over other's histograms.

  df1.calculateRDF(5.0, 0.1);
  df2.calculateCoordinationNumber();

  df1.add(df2);

  EXPECT_NO_THROW(df1.getHistogram("g(r)"));
  EXPECT_NO_THROW(df1.getHistogram("CN"));
}
