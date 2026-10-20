/**
 * @file FirstMinimumDetectorTests.cpp
 * @brief Unit tests for FirstMinimumDetector.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/FirstMinimumDetector.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <vector>

using namespace correlation::analysis;

class FirstMinimumDetectorTests : public ::testing::Test {
protected:
  // Generates synthetic g(r) representing a covalent or amorphous pair with primary peak & trough
  static auto generateSyntheticRDF(real_t peak_pos = 1.62, real_t peak_width = 0.12,
                                   real_t trough_pos = 2.15, real_t r_max_grid = 5.0,
                                   real_t dr = 0.02)
      -> std::pair<std::vector<real_t>, std::vector<real_t>> {
    std::vector<real_t> bins;
    std::vector<real_t> gr;
    const size_t n_bins = static_cast<size_t>(r_max_grid / dr);
    bins.reserve(n_bins);
    gr.reserve(n_bins);

    for (size_t i = 0; i < n_bins; ++i) {
      const real_t r = static_cast<real_t>(i) * dr;
      bins.push_back(r);

      real_t val = 0.0;
      if (r > 1.3) {
        // First coordination peak
        const real_t diff1 = (r - peak_pos) / peak_width;
        const real_t peak1 =
            static_cast<real_t>(4.5) * std::exp(-static_cast<real_t>(0.5) * diff1 * diff1);

        // Second coordination peak
        const real_t diff2 = (r - static_cast<real_t>(2.8)) / static_cast<real_t>(0.35);
        const real_t peak2 =
            static_cast<real_t>(1.8) * std::exp(-static_cast<real_t>(0.5) * diff2 * diff2);

        // Trough baseline between peaks
        val = peak1 + peak2 + static_cast<real_t>(0.1);
      }
      gr.push_back(val);
    }
    return {bins, gr};
  }
};

TEST_F(FirstMinimumDetectorTests, ReturnsErrorOnInsufficientData) {
  const std::vector<real_t> bins = {1.0, 1.1, 1.2};
  const std::vector<real_t> gr = {0.0, 1.0, 0.0};

  const auto result = FirstMinimumDetector::detectCutoffs(bins, gr);
  ASSERT_FALSE(result.has_value());
  EXPECT_EQ(result.error(), DetectionError::InsufficientData);
}

TEST_F(FirstMinimumDetectorTests, ReturnsErrorOnMismatchedSizes) {
  const std::vector<real_t> bins = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
  const std::vector<real_t> gr = {0.0, 1.0, 2.0, 3.0, 2.0, 1.0, 0.5};

  const auto result = FirstMinimumDetector::detectCutoffs(bins, gr);
  ASSERT_FALSE(result.has_value());
  EXPECT_EQ(result.error(), DetectionError::MismatchedSizes);
}

TEST_F(FirstMinimumDetectorTests, ReturnsErrorWhenNoPeakFound) {
  const std::vector<real_t> bins = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  const std::vector<real_t> gr = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

  const auto result = FirstMinimumDetector::detectCutoffs(bins, gr);
  ASSERT_FALSE(result.has_value());
  EXPECT_EQ(result.error(), DetectionError::NoPeakFound);
}

TEST_F(FirstMinimumDetectorTests, DetectsPeakOnsetAndMinimumAccurately) {
  const auto [bins, gr] = generateSyntheticRDF(1.62, 0.12, 2.15);

  const auto result = FirstMinimumDetector::detectCutoffs(bins, gr);
  ASSERT_TRUE(result.has_value()) << FirstMinimumDetector::formatError(result.error());

  const auto cutoffs = *result;
  // Primary peak should be located at approximately 1.62 Å
  EXPECT_NEAR(cutoffs.r_peak, 1.62, 0.04);
  EXPECT_GT(cutoffs.peak_height, 3.5);

  // Peak onset r_min should be between 1.3 Å and 1.55 Å (left side of peak)
  EXPECT_GE(cutoffs.r_min, 1.30);
  EXPECT_LE(cutoffs.r_min, 1.55);

  // First minimum r_max should be in the trough between peaks (~2.0 - 2.3 Å)
  EXPECT_GT(cutoffs.r_max, cutoffs.r_peak);
  EXPECT_NEAR(cutoffs.r_max, 2.15, 0.15);
  EXPECT_LT(cutoffs.min_height, cutoffs.peak_height);
}

TEST_F(FirstMinimumDetectorTests, HistogramOverloadDetectsPartial) {
  const auto [bins, gr] = generateSyntheticRDF(1.60, 0.10, 2.10);

  Histogram hist;
  hist.bins = bins;
  hist.partials["Si-O"] = gr;

  const auto result = FirstMinimumDetector::detectCutoffs(hist, "Si-O");
  ASSERT_TRUE(result.has_value());
  EXPECT_NEAR(result->r_peak, 1.60, 0.04);
  EXPECT_GT(result->r_max, 1.80);

  // Partial not found
  const auto missing = FirstMinimumDetector::detectCutoffs(hist, "O-O");
  ASSERT_FALSE(missing.has_value());
  EXPECT_EQ(missing.error(), DetectionError::PartialNotFound);
}

TEST_F(FirstMinimumDetectorTests, FormatErrorReturnsNonEmptyDescriptions) {
  EXPECT_FALSE(FirstMinimumDetector::formatError(DetectionError::InsufficientData).empty());
  EXPECT_FALSE(FirstMinimumDetector::formatError(DetectionError::MismatchedSizes).empty());
  EXPECT_FALSE(FirstMinimumDetector::formatError(DetectionError::NoPeakFound).empty());
  EXPECT_FALSE(FirstMinimumDetector::formatError(DetectionError::NoMinimumFound).empty());
  EXPECT_FALSE(FirstMinimumDetector::formatError(DetectionError::InvalidRange).empty());
  EXPECT_FALSE(FirstMinimumDetector::formatError(DetectionError::PartialNotFound).empty());
}
