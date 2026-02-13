// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <algorithm>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <iterator>

#include "../include/Cell.hpp"
#include "../include/DistributionFunctions.hpp"
#include "../include/PhysicalData.hpp"
#include "../include/Trajectory.hpp"
#include "../include/TrajectoryAnalyzer.hpp"
#include <vector>

namespace {
// Helper function to print a histogram's contents for debugging purposes.
void print_histogram(const std::string &title, const std::vector<double> &bins,
                     const std::vector<double> &values) {
  std::cout << "\n--- Histogram: " << title << " ---" << std::endl;
  std::cout << std::fixed << std::setprecision(4);
  for (size_t i = 0; i < bins.size(); ++i) {
    // Only print bins with non-zero values to keep the output concise.
    if (values.size() > i && std::abs(values[i]) > 1e-9) {
      std::cout << "Bin: " << std::setw(8) << bins[i]
                << " | Value: " << values[i] << std::endl;
    }
  }
  std::cout << "--------------------------------------\n" << std::endl;
}
} // namespace

// Test fixture for DistributionFunctions tests.
// Provides a pre-configured Cell object to reduce boilerplate in test cases.
class DistributionFunctionsTest : public ::testing::Test {
protected:
  void SetUp() override {
    // A simple cubic cell containing two atoms, 1.5 Angstroms apart.
    // This is a common setup for testing pair-based calculations.
    cell_ = Cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
    cell_.addAtom("Ar", {5.0, 5.0, 5.0});
    cell_.addAtom("Ar", {6.5, 5.0, 5.0});
  }

  void updateTrajectory() {
    trajectory_ = Trajectory();
    trajectory_.addFrame(cell_);
    trajectory_.precomputeBondCutoffs();
  }

  void updateTrajectory(const Cell &cell) {
    trajectory_ = Trajectory();
    trajectory_.addFrame(cell);
    trajectory_.precomputeBondCutoffs();
  }

  Cell cell_{};
  Trajectory trajectory_;
};

TEST_F(DistributionFunctionsTest, CalculateRDFThrowsOnInvalidParameters) {
  // Arrange
  updateTrajectory();
  DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffs());

  // Act & Assert
  EXPECT_THROW(df.calculateRDF(5.0, 0.0), std::invalid_argument);
  EXPECT_THROW(df.calculateRDF(0.0, 0.1), std::invalid_argument);
}

TEST_F(DistributionFunctionsTest, RDFPeakPositionIsCorrect) {
  // Arrange
  updateTrajectory();
  DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffs());
  const double bin_width = 0.1;
  const double expected_distance = 1.5;

  // Act
  df.calculateRDF(5.0, bin_width);
  const auto &rdf_hist = df.getHistogram("g(r)");
  const auto &total_rdf = rdf_hist.partials.at("Ar-Ar");

  // Assert: Find the peak of the RDF and verify its position.
  auto max_it = std::max_element(total_rdf.begin(), total_rdf.end());
  size_t peak_index = std::distance(total_rdf.begin(), max_it);

  double peak_position = rdf_hist.bins[peak_index];

  EXPECT_NEAR(peak_position, expected_distance, bin_width);
}

TEST_F(DistributionFunctionsTest, PADPeakPositionIsCorrectForWater) {
  // Arrange: A water-like structure with a known ~109.5 degree angle.
  Cell water_cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  water_cell.addAtom("O", {0.0, 0.0, 0.0});
  water_cell.addAtom("H", {1.0, 0.0, 0.0});
  water_cell.addAtom("H",
                     {std::cos(1.916), std::sin(1.916), 0.0}); // ~109.5 deg

  updateTrajectory(water_cell);
  DistributionFunctions df(water_cell, 5.0, trajectory_.getBondCutoffs());
  const double bin_width = 1.0; // 1-degree bins

  // Act
  df.calculatePAD(bin_width);
  const auto &pad_hist = df.getHistogram("f(theta)");
  const auto &hoh_pad = pad_hist.partials.at("H-O-H");

  // Assert
  auto max_it = std::max_element(hoh_pad.begin(), hoh_pad.end());
  size_t peak_index = std::distance(hoh_pad.begin(), max_it);
  double peak_angle = pad_hist.bins[peak_index];

  EXPECT_NEAR(peak_angle, 109.5, bin_width * 2.0); // Allow for binning error
}

TEST_F(DistributionFunctionsTest, SmoothAllUpdatesSmoothedPartials) {
  // Arrange
  updateTrajectory();
  DistributionFunctions df(cell_, 5.0, trajectory_.getBondCutoffs());
  df.calculateRDF(5.0, 0.1);

  // Act
  df.smoothAll(0.2);
  const auto &rdf_hist = df.getHistogram("g(r)");

  // Assert
  ASSERT_FALSE(rdf_hist.smoothed_partials.empty());
  ASSERT_TRUE(rdf_hist.smoothed_partials.count("Ar-Ar"));

  const auto &raw_data = rdf_hist.partials.at("Ar-Ar");
  const auto &smoothed_data = rdf_hist.smoothed_partials.at("Ar-Ar");

  ASSERT_EQ(raw_data.size(), smoothed_data.size());

  // A simple check: the peak of the smoothed data should be lower than the raw
  // data's peak.
  double raw_max = *std::max_element(raw_data.begin(), raw_data.end());
  double smoothed_max =
      *std::max_element(smoothed_data.begin(), smoothed_data.end());

  EXPECT_LE(smoothed_max, raw_max);
}

TEST_F(DistributionFunctionsTest, CoordinationNumberDistributionIsCorrect) {
  // Arrange: Create a structure with a known, simple coordination environment.
  Cell test_cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  // A central "C" atom
  test_cell.addAtom("C", {10.0, 10.0, 10.0});
  // Four "H" atoms tetrahedrally coordinated around the "C"
  test_cell.addAtom("H", {11.0, 10.0, 10.0});
  test_cell.addAtom("H", {10.0, 11.0, 10.0});
  test_cell.addAtom("H", {10.0, 10.0, 11.0});
  test_cell.addAtom("H", {10.5, 10.5, 10.5}); // Not a perfect tetrahedron
  // A lone "O" atom, not bonded to anything
  test_cell.addAtom("O", {15.0, 15.0, 15.0});

  // A cutoff that includes the C-H bonds but excludes everything else.
  // The C-H distance is 1.0, C-O is ~8.6
  updateTrajectory(test_cell);
  DistributionFunctions df(test_cell, 3.0, trajectory_.getBondCutoffs());

  // Act
  df.calculateCoordinationNumber();
  const auto &cn_hist = df.getHistogram("CN");

  // Assert C-H coordination
  const auto &c_h_cn = cn_hist.partials.at("C-H");
  // There is 1 Carbon atom, and it has 4 Hydrogen neighbors.
  // So, the histogram should have a value of 1 at bin 4.
  EXPECT_EQ(c_h_cn.size(), 7); // Bins for CN=0,1,2,3,4,5,6
  EXPECT_EQ(c_h_cn[4], 1);
  EXPECT_EQ(c_h_cn[0] + c_h_cn[1] + c_h_cn[2] + c_h_cn[3], 0);

  // Assert H-C coordination
  const auto &h_c_cn = cn_hist.partials.at("H-C");
  // There are 4 Hydrogen atoms, and each has 1 Carbon neighbor.
  // So, the histogram should have a value of 4 at bin 1.
  EXPECT_EQ(h_c_cn.size(), 7); // Must be padded to the max CN+2
  EXPECT_EQ(h_c_cn[1], 4);
  EXPECT_EQ(h_c_cn[0] + h_c_cn[2] + h_c_cn[3] + h_c_cn[4], 0);

  // Assert that non-bonded pairs do not appear or are empty.
  EXPECT_EQ(cn_hist.partials.count("C-O"), 0);

  // Assert Element-Any columns
  // C-Any: equivalent to C-H since C only has H neighbors
  ASSERT_TRUE(cn_hist.partials.count("C-Any"));
  const auto &c_any = cn_hist.partials.at("C-Any");
  EXPECT_EQ(c_any[4], 1);

  // H-Any: equivalent to H-C since H only has C neighbors
  ASSERT_TRUE(cn_hist.partials.count("H-Any"));
  const auto &h_any = cn_hist.partials.at("H-Any");
  EXPECT_EQ(h_any[1], 4);

  // O-Any: O has no neighbors, so no partials started with O-.
  // Thus O-Any should not be created.
  EXPECT_EQ(cn_hist.partials.count("O-Any"), 0);

  // Assert Any-Any column: Sum of all Element-Any columns
  // Bin 4: from C-Any (1)
  // Bin 1: from H-Any (4)
  ASSERT_TRUE(cn_hist.partials.count("Any-Any"));
  const auto &any_any = cn_hist.partials.at("Any-Any");
  EXPECT_EQ(any_any[4], 1);
  EXPECT_EQ(any_any[1], 4);
}

TEST_F(DistributionFunctionsTest, SQPeakPositionIsCorrect) {
  // Arrange
  // Using the fixture's cell with two Ar atoms 1.5 Å apart.
  updateTrajectory();
  DistributionFunctions df(cell_, 10.0, trajectory_.getBondCutoffs());
  const double q_bin_width = 0.01;
  const double expected_peak_q = 0.7; // 2.0 * constants::pi / 1.5; // ~4.18

  // Act
  // S(Q) depends on g(r), so we must calculate it first.
  df.calculateRDF(10.0, 0.01);
  // df.smooth("g(r)", 0.1, KernelType::Gaussian);
  df.calculateSQ(10.0, q_bin_width, 8.0);
  const auto &g_hist = df.getHistogram("g(r)");
  const auto &total_g = g_hist.partials.at("Total");

  const auto &sq_hist = df.getHistogram("S(Q)");
  const auto &total_sq = sq_hist.partials.at("Total");

  // Assert: Find the peak of S(Q) and verify its position.
  // We'll ignore the very first few bins as S(Q) can be noisy at Q->0.
  auto search_start = total_sq.begin() + 5;
  auto max_it = std::max_element(search_start, total_sq.end());
  size_t peak_index = std::distance(total_sq.begin(), max_it);

  double peak_position = sq_hist.bins[peak_index];

  // The peak should be near 2*pi/r
  EXPECT_NEAR(peak_position, expected_peak_q,
              q_bin_width * 5.0); // Allow a generous tolerance
}

TEST_F(DistributionFunctionsTest, AccumulationAndScalingWorks) {
  // Arrange: Create two different cells
  Cell cell1 = cell_; // Default cell from fixture (Ar pair at 1.5)

  Cell cell2({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell2.addAtom("Ar", {5.0, 5.0, 5.0});
  cell2.addAtom("Ar", {7.0, 5.0, 5.0}); // Distance 2.0

  Trajectory traj1;
  traj1.addFrame(cell1);
  traj1.precomputeBondCutoffs();
  Trajectory traj2;
  traj2.addFrame(cell2);
  traj2.precomputeBondCutoffs();

  DistributionFunctions df1(cell1, 5.0, traj1.getBondCutoffs());
  DistributionFunctions df2(cell2, 5.0, traj2.getBondCutoffs());

  // Calculate RDF for both
  df1.calculateRDF(5.0, 0.1);
  df2.calculateRDF(5.0, 0.1);

  const auto &h1 = df1.getHistogram("g(r)").partials.at("Ar-Ar");
  const auto &h2 = df2.getHistogram("g(r)").partials.at("Ar-Ar");

  // Expected average
  std::vector<double> expected(h1.size());
  for (size_t i = 0; i < h1.size(); ++i) {
    expected[i] = (h1[i] + h2[i]) / 2.0;
  }

  // Act: Add df2 to df1 and scale
  df1.add(df2);
  df1.scale(0.5);

  const auto &h_avg = df1.getHistogram("g(r)").partials.at("Ar-Ar");

  // Assert
  ASSERT_EQ(h_avg.size(), expected.size());
  for (size_t i = 0; i < h_avg.size(); ++i) {
    EXPECT_NEAR(h_avg[i], expected[i], 1e-9);
  }
}

TEST_F(DistributionFunctionsTest, ComputeMeanMatchesSequential) {
  // Arrange
  Cell cell1 = cell_;
  Cell cell2({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell2.addAtom("Ar", {5.0, 5.0, 5.0});
  cell2.addAtom("Ar", {7.0, 5.0, 5.0});

  Trajectory traj;
  traj.addFrame(cell1);
  traj.addFrame(cell2);
  traj.precomputeBondCutoffs();

  std::vector<std::vector<double>> cutoffs = traj.getBondCutoffs();

  TrajectoryAnalyzer analyzer(traj, 10.0, cutoffs);

  AnalysisSettings settings;
  settings.r_max = 5.0;
  settings.r_bin_width = 0.1;
  settings.smoothing = false; // Disable smoothing to compare raw values easily
  settings.angle_bin_width = 1.0;
  settings.q_max = 10.0;

  // Act
  auto result = DistributionFunctions::computeMean(traj, analyzer, 0, settings);

  // Manual calculation for comparison
  DistributionFunctions df1(cell1, 0.0, cutoffs);
  df1.setStructureAnalyzer(analyzer.getAnalyzers()[0].get());
  df1.calculateRDF(5.0, 0.1);

  DistributionFunctions df2(cell2, 0.0, cutoffs);
  df2.setStructureAnalyzer(analyzer.getAnalyzers()[1].get());
  df2.calculateRDF(5.0, 0.1);

  df1.add(df2);
  df1.scale(0.5);

  // Assert
  const auto &h_auto = result->getHistogram("g(r)").partials.at("Ar-Ar");
  const auto &h_manual = df1.getHistogram("g(r)").partials.at("Ar-Ar");

  ASSERT_EQ(h_auto.size(), h_manual.size());
  for (size_t i = 0; i < h_auto.size(); ++i) {
    EXPECT_NEAR(h_auto[i], h_manual[i], 1e-9);
  }
}

TEST_F(DistributionFunctionsTest, CalculateVACFReturnsCorrectValues) {
  // Arrange
  Cell c1({10, 10, 10, 90, 90, 90});
  c1.addAtom("Ar", {0, 0, 0});
  Cell c2({10, 10, 10, 90, 90, 90});
  c2.addAtom("Ar", {1, 0, 0});
  Cell c3({10, 10, 10, 90, 90, 90});
  c3.addAtom("Ar", {2, 0, 0});

  Trajectory traj;
  traj.addFrame(c1);
  traj.addFrame(c2);
  traj.addFrame(c3);
  traj.setTimeStep(1.0);
  traj.calculateVelocities();

  // Dummy cutoffs
  std::vector<std::vector<double>> cutoffs = {{0.0}};
  DistributionFunctions df(c1, 0.0, cutoffs);

  // Act
  df.calculateVACF(traj, 2);

  // Assert
  ASSERT_NO_THROW(df.getHistogram("VACF"));
  const auto &hist = df.getHistogram("VACF");
  const auto &total = hist.partials.at("Total");

  // With constant velocity (1,0,0), VACF should be constant 1.0
  EXPECT_EQ(total.size(), 3);
  EXPECT_NEAR(total[0], 1.0, 1e-6);
  EXPECT_NEAR(total[1], 1.0, 1e-6);

  // Check Normalized VACF
  ASSERT_NO_THROW(df.getHistogram("Normalized VACF"));
  const auto &norm_hist = df.getHistogram("Normalized VACF");
  const auto &norm_total = norm_hist.partials.at("Total");
  EXPECT_NEAR(norm_total[0], 1.0, 1e-6);
  EXPECT_NEAR(norm_total[1], 1.0, 1e-6);
}

TEST_F(DistributionFunctionsTest, CalculateVDOSGeneratesExtraUnits) {
  // Arrange
  Cell c1({10, 10, 10, 90, 90, 90});
  c1.addAtom("Ar", {0, 0, 0});
  Trajectory traj;
  traj.addFrame(c1);
  traj.addFrame(c1); // 2 frames needed minimum
  traj.setTimeStep(1.0);
  traj.calculateVelocities();

  // Dummy cutoffs
  std::vector<std::vector<double>> cutoffs = {{0.0}};
  DistributionFunctions df(c1, 0.0, cutoffs);

  df.calculateVACF(traj, 1);

  // Act
  df.calculateVDOS();

  // Assert
  ASSERT_NO_THROW(df.getHistogram("VDOS"));
  const auto &hist = df.getHistogram("VDOS");

  EXPECT_TRUE(hist.partials.count("Frequency (cm-1)"));
  EXPECT_TRUE(hist.partials.count("Frequency (meV)"));

  const auto &freq_thz = hist.bins; // Bins are Frequency (THz)
  const auto &freq_cm = hist.partials.at("Frequency (cm-1)");
  const auto &freq_mev = hist.partials.at("Frequency (meV)");

  ASSERT_EQ(freq_thz.size(), freq_cm.size());
  ASSERT_EQ(freq_thz.size(), freq_mev.size());

  // Check conversion factors for a non-zero frequency point
  // Index 0 is usually 0 frequency or negative max depending on implementation
  // Let's check a point
  if (freq_thz.size() > 1) {
    // Find a non-zero frequency
    for (size_t i = 0; i < freq_thz.size(); ++i) {
      if (std::abs(freq_thz[i]) > 1e-6) {
        EXPECT_NEAR(freq_cm[i], freq_thz[i] * constants::THz_to_cmInv, 1e-5);
        EXPECT_NEAR(freq_mev[i], freq_thz[i] * constants::THz_to_meV, 1e-5);
        break;
      }
    }
  }
}
