// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <algorithm>
#include <gtest/gtest.h>
#include <iomanip>  // Required for std::setw and std::fixed
#include <iostream> // Required for std::cout
#include <iterator>

#include "../include/Cell.hpp"
#include "../include/DistributionFunctions.hpp"

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

  Cell cell_{};
};

TEST_F(DistributionFunctionsTest, CalculateRDFThrowsOnInvalidParameters) {
  // Arrange
  auto neighbors_ = NeighborList(cell_);
  DistributionFunctions df(cell_, neighbors_);

  // Act & Assert
  EXPECT_THROW(df.calculateRDF(5.0, 0.0), std::invalid_argument);
  EXPECT_THROW(df.calculateRDF(0.0, 0.1), std::invalid_argument);
}

TEST_F(DistributionFunctionsTest, RDFPeakPositionIsCorrect) {
  // Arrange
  auto neighbors_ = NeighborList(cell_);
  DistributionFunctions df(cell_, neighbors_);
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
  auto neighbors_ = NeighborList(water_cell);
  DistributionFunctions df(water_cell, neighbors_);
  const double bin_width = 1.0; // 1-degree bins

  // Act
  df.calculatePAD(180.0, bin_width);
  const auto &pad_hist = df.getHistogram("f(theta)");
  const auto &hoh_pad = pad_hist.partials.at("H-O-H");

  // --- DEBUGGING STEP ---
  // This will print the contents of the H-O-H partial to the console.
  print_histogram("H-O-H Partial", pad_hist.bins, hoh_pad);
  // --- END DEBUGGING STEP ---

  // Assert
  auto max_it = std::max_element(hoh_pad.begin(), hoh_pad.end());
  size_t peak_index = std::distance(hoh_pad.begin(), max_it);
  double peak_angle = pad_hist.bins[peak_index];

  EXPECT_NEAR(peak_angle, 109.5, bin_width * 2.0); // Allow for binning error
}

TEST_F(DistributionFunctionsTest, SmoothAllUpdatesSmoothedPartials) {
  // Arrange
  auto neighbors_ = NeighborList(cell_);
  DistributionFunctions df(cell_, neighbors_);
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

  EXPECT_LT(smoothed_max, raw_max);
}
