// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "calculators/HyperuniformityCalculator.hpp"
#include "core/Cell.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <random>
#include <vector>

namespace correlation::calculators {

namespace {
/**
 * @brief Struct to hold the parameters for creating a random cell.
 */

struct CreateRandomCellParams {
  double box_length;
  size_t num_atoms;
};

/**
 * @brief Helper to create a simple cubic (SC) lattice cell.
 *
 * Creates a cubic box of side L with atoms placed on a regular grid
 * with spacing `a` (lattice constant). The number of atoms per axis
 * is n = L / a.
 */
correlation::core::Cell createSCLattice(real_t box_length, int atoms_per_axis) {
  real_t const lattice_constant = box_length / static_cast<real_t>(atoms_per_axis);
  correlation::core::Cell cell({box_length, box_length, box_length, static_cast<real_t>(90.0),
                                static_cast<real_t>(90.0), static_cast<real_t>(90.0)});
  for (int ix = 0; ix < atoms_per_axis; ++ix) {
    for (int iy = 0; iy < atoms_per_axis; ++iy) {
      for (int iz = 0; iz < atoms_per_axis; ++iz) {
        real_t x_pos = (static_cast<real_t>(ix) + static_cast<real_t>(0.5)) * lattice_constant;
        real_t y_pos = (static_cast<real_t>(iy) + static_cast<real_t>(0.5)) * lattice_constant;
        real_t z_pos = (static_cast<real_t>(iz) + static_cast<real_t>(0.5)) * lattice_constant;
        cell.addAtom("Ar", {x_pos, y_pos, z_pos});
      }
    }
  }
  return cell;
}

/**
 * @brief Helper to create a cell with randomly placed atoms (Poisson process).
 *
 * Uses a fixed seed for reproducibility.
 */
correlation::core::Cell createRandomCell(CreateRandomCellParams params) {
  correlation::core::Cell cell({static_cast<real_t>(params.box_length), static_cast<real_t>(params.box_length),
                                static_cast<real_t>(params.box_length), static_cast<real_t>(90.0),
                                static_cast<real_t>(90.0), static_cast<real_t>(90.0)});
  std::mt19937_64 rng(12345); // NOLINT(cert-msc51-cpp, cert-msc32-c)
  std::uniform_real_distribution<double> dist(0.0, params.box_length);
  for (size_t i = 0; i < params.num_atoms; ++i) {
    cell.addAtom("Ar",
                 {static_cast<real_t>(dist(rng)), static_cast<real_t>(dist(rng)), static_cast<real_t>(dist(rng))});
  }
  return cell;
}

class HyperuniformityCalculatorTests : public ::testing::Test {};

} // namespace

// =============================================================================
// Basic functionality tests
// =============================================================================

TEST_F(HyperuniformityCalculatorTests, ReturnsEmptyForEmptyCell) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  auto results = HyperuniformityCalculator::calculate(cell, {.num_samples = 100, .r_bin_width = 0.5});
  EXPECT_TRUE(results.empty());
}

TEST_F(HyperuniformityCalculatorTests, ReturnsEmptyForTinyBox) {
  // Box too small: L/2 = 1.5, which is < r_min=2.0
  correlation::core::Cell cell({3.0, 3.0, 3.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {1.5, 1.5, 1.5});
  auto results = HyperuniformityCalculator::calculate(cell, {.num_samples = 100, .r_bin_width = 0.5});
  EXPECT_TRUE(results.empty());
}

TEST_F(HyperuniformityCalculatorTests, ProducesExpectedHistograms) {
  auto cell = createSCLattice(20.0, 4);
  auto results = HyperuniformityCalculator::calculate(cell, {.num_samples = 500, .r_bin_width = 1.0});

  ASSERT_TRUE(results.contains("sigma2_N"));
  ASSERT_TRUE(results.contains("chi_H"));

  const auto &sigma2 = results.at("sigma2_N");
  const auto &chi_h = results.at("chi_H");

  // Bins should be populated
  EXPECT_FALSE(sigma2.bins.empty());
  EXPECT_FALSE(chi_h.bins.empty());
  EXPECT_EQ(sigma2.bins.size(), chi_h.bins.size());

  // Should have "Total" partial
  EXPECT_TRUE(sigma2.partials.contains("Total"));
  EXPECT_TRUE(chi_h.partials.contains("Total"));

  // Variance should be non-negative everywhere
  for (double variance : sigma2.partials.at("Total")) {
    EXPECT_GE(variance, -1e-9); // Allow tiny numerical errors
  }
}

TEST_F(HyperuniformityCalculatorTests, ThrowsOnInvalidBinWidth) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});
  EXPECT_THROW(HyperuniformityCalculator::calculate(cell, {100, 0.0}), std::invalid_argument);
  EXPECT_THROW(HyperuniformityCalculator::calculate(cell, {100, -1.0}), std::invalid_argument);
}

TEST_F(HyperuniformityCalculatorTests, ThrowsOnZeroSamples) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});
  EXPECT_THROW(HyperuniformityCalculator::calculate(cell, {0, 0.5}), std::invalid_argument);
}

// =============================================================================
// Variance scaling tests
// =============================================================================

/**
 * For a Poisson (random) point process, the number variance scales as:
 *   σ²_N(R) ∝ R³  (volume-dominated)
 *
 * We verify this by fitting log(σ²_N) vs log(R) and checking the slope
 * is close to 3.
 */
TEST_F(HyperuniformityCalculatorTests, RandomPointsVarianceScalesAsR3) {
  auto cell = createRandomCell({.box_length = 30.0, .num_atoms = 500});
  // Use more samples for statistical stability
  auto results = HyperuniformityCalculator::calculate(cell, {.num_samples = 5000, .r_bin_width = 0.5});

  ASSERT_TRUE(results.contains("sigma2_N"));
  const auto &sigma2 = results.at("sigma2_N");
  const auto &bins = sigma2.bins;
  const auto &total = sigma2.partials.at("Total");

  // Fit log(σ²) = slope * log(R) + intercept using least squares
  // Only use bins where variance > 0 and R is in a reasonable range
  double sum_x = 0.0;
  double sum_y = 0.0;
  double sum_xx = 0.0;
  double sum_xy = 0.0;
  int count = 0;
  for (size_t k = 0; k < bins.size(); ++k) {
    if (total[k] > 0.0 && bins[k] > 2.5 && bins[k] < 12.0) {
      double log_r = std::log(bins[k]);
      double log_var = std::log(total[k]);
      sum_x += log_r;
      sum_y += log_var;
      sum_xx += log_r * log_r;
      sum_xy += log_r * log_var;
      ++count;
    }
  }

  ASSERT_GT(count, 3) << "Not enough data points for fitting";

  double const slope = (count * sum_xy - sum_x * sum_y) / (count * sum_xx - sum_x * sum_x);

  // For Poisson points, slope should be ~3.0
  // Allow generous tolerance due to finite-size effects and random sampling
  EXPECT_GT(slope, 2.0) << "Slope too low for Poisson points (expected ~3.0, got " << slope << ")";
  EXPECT_LT(slope, 4.5) << "Slope too high for Poisson points (expected ~3.0, got " << slope << ")";
}

/**
 * For a perfect lattice (hyperuniform system), the number variance scales as:
 *   σ²_N(R) ∝ R²  (surface-dominated)
 *
 * We verify this by fitting log(σ²_N) vs log(R) and checking the slope
 * is close to 2.
 */
TEST_F(HyperuniformityCalculatorTests, LatticeVarianceScalesAsR2) {
  // Large SC lattice: 8³ = 512 atoms in 40 Å box
  auto cell = createSCLattice(40.0, 8);
  auto results = HyperuniformityCalculator::calculate(cell, {5000, 0.5});

  ASSERT_TRUE(results.contains("sigma2_N"));
  const auto &sigma2 = results.at("sigma2_N");
  const auto &bins = sigma2.bins;
  const auto &total = sigma2.partials.at("Total");

  // Fit log(σ²) = slope * log(R) + intercept
  // Use intermediate R range to avoid boundary/commensurability artifacts
  double sum_x = 0;
  double sum_y = 0;
  double sum_xx = 0;
  double sum_xy = 0;
  int count = 0;
  for (size_t k = 0; k < bins.size(); ++k) {
    if (total[k] > 0.0 && bins[k] > 4.0 && bins[k] < 15.0) {
      double log_r = std::log(bins[k]);
      double log_var = std::log(total[k]);
      sum_x += log_r;
      sum_y += log_var;
      sum_xx += log_r * log_r;
      sum_xy += log_r * log_var;
      ++count;
    }
  }

  ASSERT_GT(count, 3) << "Not enough data points for fitting";

  double const slope = (count * sum_xy - sum_x * sum_y) / (count * sum_xx - sum_x * sum_x);

  // For a hyperuniform system, slope should be ~2.0
  // Allow generous tolerance for finite-size effects
  EXPECT_GT(slope, 1.0) << "Slope too low for lattice (expected ~2.0, got " << slope << ")";
  EXPECT_LT(slope, 3.5) << "Slope too high for lattice (expected ~2.0, got " << slope << ")";
}

/**
 * The lattice should have a lower slope than random points, confirming
 * that the hyperuniformity metric can distinguish ordered vs disordered systems.
 */
TEST_F(HyperuniformityCalculatorTests, LatticeHasLowerSlopeThanRandom) {
  // Create both systems with same density
  auto lattice_cell = createSCLattice(30.0, 6); // 216 atoms
  auto random_cell = createRandomCell({.box_length = 30.0, .num_atoms = 216});

  auto lattice_results = HyperuniformityCalculator::calculate(lattice_cell, {.num_samples = 3000, .r_bin_width = 0.5});
  auto random_results = HyperuniformityCalculator::calculate(random_cell, {.num_samples = 3000, .r_bin_width = 0.5});

  ASSERT_TRUE(lattice_results.contains("sigma2_N"));
  ASSERT_TRUE(random_results.contains("sigma2_N"));

  // Compute slopes for both
  auto compute_slope = [](const analysis::Histogram &sigma2) -> double {
    double sum_x = 0;
    double sum_y = 0;
    double sum_xx = 0;
    double sum_xy = 0;
    int count = 0;
    for (size_t k = 0; k < sigma2.bins.size(); ++k) {
      double var = sigma2.partials.at("Total")[k];
      double r_bin = sigma2.bins[k];
      if (var > 0.0 && r_bin > 3.0 && r_bin < 12.0) {
        double log_r_bin = std::log(r_bin);
        double log_var = std::log(var);
        sum_x += log_r_bin;
        sum_y += log_var;
        sum_xx += log_r_bin * log_r_bin;
        sum_xy += log_r_bin * log_var;
        ++count;
      }
    }
    if (count < 3) {
      return -1.0;
    }
    return (count * sum_xy - sum_x * sum_y) / (count * sum_xx - sum_x * sum_x);
  };

  double const lattice_slope = compute_slope(lattice_results.at("sigma2_N"));
  double const random_slope = compute_slope(random_results.at("sigma2_N"));

  ASSERT_GT(lattice_slope, 0.0) << "Lattice slope computation failed";
  ASSERT_GT(random_slope, 0.0) << "Random slope computation failed";

  EXPECT_LT(lattice_slope, random_slope) << "Lattice (slope=" << lattice_slope
                                         << ") should have lower variance scaling "
                                         << "than random (slope=" << random_slope << ")";
}

// =============================================================================
// Histogram metadata tests
// =============================================================================

TEST_F(HyperuniformityCalculatorTests, HistogramMetadataIsCorrect) {
  auto cell = createSCLattice(20.0, 4);
  auto results = HyperuniformityCalculator::calculate(cell, {100, 1.0});

  ASSERT_TRUE(results.contains("sigma2_N"));
  ASSERT_TRUE(results.contains("chi_H"));

  const auto &sigma2 = results.at("sigma2_N");
  EXPECT_EQ(sigma2.x_label, "R");
  EXPECT_EQ(sigma2.y_label, "σ²_N(R)");
  EXPECT_EQ(sigma2.x_unit, "Å");
  EXPECT_EQ(sigma2.file_suffix, "_sigma2_N");
  EXPECT_FALSE(sigma2.title.empty());
  EXPECT_FALSE(sigma2.description.empty());

  const auto &chi_h = results.at("chi_H");
  EXPECT_EQ(chi_h.x_label, "R");
  EXPECT_EQ(chi_h.y_label, "χ_H(R)");
  EXPECT_EQ(chi_h.x_unit, "Å");
  EXPECT_EQ(chi_h.file_suffix, "_chi_H");
  EXPECT_FALSE(chi_h.title.empty());
  EXPECT_FALSE(chi_h.description.empty());
}

// =============================================================================
// BaseCalculator interface tests
// =============================================================================

TEST_F(HyperuniformityCalculatorTests, CalculatorInterfaceIsCorrect) {
  HyperuniformityCalculator calc;
  EXPECT_FALSE(calc.getName().empty());
  EXPECT_FALSE(calc.getShortName().empty());
  EXPECT_EQ(calc.getGroup(), "Advanced");
  EXPECT_TRUE(calc.isFrameCalculator());
  EXPECT_FALSE(calc.isTrajectoryCalculator());
  EXPECT_FALSE(calc.getDescription().empty());
}

} // namespace correlation::calculators
