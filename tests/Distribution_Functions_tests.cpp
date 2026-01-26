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
#include <vector>
#include "../include/PhysicalData.hpp"
#include "../include/Trajectory.hpp"

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

  void updateTrajectory(const Cell& cell) {
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
  df.calculatePAD(180.0, bin_width);
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
