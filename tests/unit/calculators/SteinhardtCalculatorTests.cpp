// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/StructureAnalyzer.hpp"
#include "calculators/SteinhardtCalculator.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>
#include <numbers>

namespace correlation::analysis {
namespace {

class SteinhardtCalculatorTests : public ::testing::Test {
protected:
  static void checkOutputs(const std::map<std::string, Histogram> &hists, double expected_Q4, double expected_Q6,
                           double expected_W6_hat) {
    const auto &hist_Q4 = hists.at("Q4").partials.at("Total");
    const auto &hist_Q6 = hists.at("Q6").partials.at("Total");
    const auto &hist_W6 = hists.at("W6_hat").partials.at("Total");

    double q4_val = 0;
    double q6_val = 0;
    double w6_val = 0;

    // Find the non-zero bins
    size_t const q4_bins = 100;
    double const d_Q = 1.0 / q4_bins;
    size_t const w6_bins = 100;
    double const W_min = -0.2;
    double const W_max = 0.2;
    double const d_W = (W_max - W_min) / w6_bins;

    for (size_t bin = 0; bin < q4_bins; ++bin) {
      if (hist_Q4[bin] > 0) {
        q4_val = (static_cast<double>(bin) + 0.5) * d_Q;
      }
      if (hist_Q6[bin] > 0) {
        q6_val = (static_cast<double>(bin) + 0.5) * d_Q;
      }
    }
    for (size_t bin = 0; bin < w6_bins; ++bin) {
      if (hist_W6[bin] > 0) {
        w6_val = W_min + (static_cast<double>(bin) + 0.5) * d_W;
      }
    }

    EXPECT_NEAR(q4_val, expected_Q4, 0.015);
    EXPECT_NEAR(q6_val, expected_Q6, 0.015);
    EXPECT_NEAR(w6_val, expected_W6_hat, 0.015);
  }
};

TEST_F(SteinhardtCalculatorTests, SimpleCubic) {
  correlation::core::Cell cell({1.0, 1.0, 1.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {0.5, 0.5, 0.5});

  // ignore_periodic_self_interactions = false
  StructureAnalyzer const analyzer(cell, 1.1, {{1.1 * 1.1}}, false);
  auto hists = correlation::calculators::SteinhardtCalculator::calculate(cell, &analyzer);

  checkOutputs(hists, 0.764, 0.354, 0.013);
}

TEST_F(SteinhardtCalculatorTests, BCC) {
  correlation::core::Cell cell({1.0, 1.0, 1.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {0.0, 0.0, 0.0});
  cell.addAtom("Ar", {0.5, 0.5, 0.5});

  StructureAnalyzer const analyzer(cell, 1.1, {{1.1 * 1.1}},
                                   false); // dist = sqrt(0.75) ~ 0.866 and 1.0
  auto hists = correlation::calculators::SteinhardtCalculator::calculate(cell, &analyzer);

  checkOutputs(hists, 0.036, 0.511, 0.013);
}

TEST_F(SteinhardtCalculatorTests, FCC) {
  correlation::core::Cell cell({1.0, 1.0, 1.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {0.0, 0.0, 0.0});
  cell.addAtom("Ar", {0.5, 0.5, 0.0});
  cell.addAtom("Ar", {0.5, 0.0, 0.5});
  cell.addAtom("Ar", {0.0, 0.5, 0.5});

  StructureAnalyzer const analyzer(cell, 0.8, {{0.8 * 0.8}},
                                   false); // dist = sqrt(0.5) ~ 0.707
  auto hists = correlation::calculators::SteinhardtCalculator::calculate(cell, &analyzer);

  checkOutputs(hists, 0.191, 0.575, -0.013); // W6 for FCC is approx -0.013
}

TEST_F(SteinhardtCalculatorTests, Icosahedral) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});

  double const phi = std::numbers::phi;
  double const scale = 1.0 / std::sqrt(1.0 + phi * phi); // Normalize dist to 1.0
  double const phi_scale = phi * scale;

  cell.addAtom("Ar", {5.0, 5.0 + scale, 5.0 + phi_scale});
  cell.addAtom("Ar", {5.0, 5.0 + scale, 5.0 - phi_scale});
  cell.addAtom("Ar", {5.0, 5.0 - scale, 5.0 + phi_scale});
  cell.addAtom("Ar", {5.0, 5.0 - scale, 5.0 - phi_scale});

  cell.addAtom("Ar", {5.0 + scale, 5.0 + phi_scale, 5.0});
  cell.addAtom("Ar", {5.0 + scale, 5.0 - phi_scale, 5.0});
  cell.addAtom("Ar", {5.0 - scale, 5.0 + phi_scale, 5.0});
  cell.addAtom("Ar", {5.0 - scale, 5.0 - phi_scale, 5.0});

  cell.addAtom("Ar", {5.0 + phi_scale, 5.0, 5.0 + scale});
  cell.addAtom("Ar", {5.0 + phi_scale, 5.0, 5.0 - scale});
  cell.addAtom("Ar", {5.0 - phi_scale, 5.0, 5.0 + scale});
  cell.addAtom("Ar", {5.0 - phi_scale, 5.0, 5.0 - scale});

  // Edge length is ~1.05. Using cutoff 1.02 ensures surface atoms only see
  // center. Thus they will have 1 neighbor, Ql=1.0, and be excluded from
  // histogram!
  StructureAnalyzer const analyzer(cell, 1.02, {{1.02 * 1.02}}, true);
  auto hists = correlation::calculators::SteinhardtCalculator::calculate(cell, &analyzer);

  checkOutputs(hists, 0.000, 0.663,
               -0.169); // W6_hat for Icosahedral is approx -0.1697
}

TEST_F(SteinhardtCalculatorTests, HandlesAcosNumericalNoiseSafely) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});
  cell.addAtom("Ar", {5.0, 5.0, 6.000000000000001});
  cell.addAtom("Ar", {5.0, 5.0, 4.0});

  StructureAnalyzer const analyzer(cell, 1.1, {{1.1 * 1.1}}, false);
  ASSERT_NO_THROW({
    auto hists = correlation::calculators::SteinhardtCalculator::calculate(cell, &analyzer);
    EXPECT_FALSE(hists.empty());
  });
}

TEST_F(SteinhardtCalculatorTests, EmptySystemOrNoNeighborsFillsPartialsWithZeros) {
  // Cell with 1 atom (no neighbors)
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});

  StructureAnalyzer const analyzer(cell, 1.1, {{1.1 * 1.1}}, false);
  auto hists = correlation::calculators::SteinhardtCalculator::calculate(cell, &analyzer);

  EXPECT_FALSE(hists.empty());
  for (const auto &name : {"Q4", "Q6", "W6_hat"}) {
    ASSERT_TRUE(hists.count(name));
    const auto &hist = hists.at(name);
    ASSERT_TRUE(hist.partials.count("Total"));

    // Check that all bins in Total are exactly 0
    const auto &total = hist.partials.at("Total");
    for (double const val : total) {
      EXPECT_DOUBLE_EQ(val, 0.0);
    }
  }
}

TEST_F(SteinhardtCalculatorTests, SphericalHarmonics) {
  using correlation::calculators::SteinhardtCalculator;

  // L = 0, M = 0: Y_0^0 = 0.5 * sqrt(1/pi) ~ 0.28209479
  {
    auto val = SteinhardtCalculator::sphericalHarmonic(0, 0, {.theta = 0.5, .phi = 0.2});
    EXPECT_NEAR(val.real(), 0.28209479177, 1e-8);
    EXPECT_NEAR(val.imag(), 0.0, 1e-8);
  }

  // L = 1, M = 0: Y_1^0 = 0.5 * sqrt(3/pi) * cos(theta) ~ 0.4886025 * cos(theta)
  {
    double const theta = 1.0;
    auto val = SteinhardtCalculator::sphericalHarmonic(1, 0, {.theta = theta, .phi = 0.5});
    EXPECT_NEAR(val.real(), 0.4886025119 * std::cos(theta), 1e-8);
    EXPECT_NEAR(val.imag(), 0.0, 1e-8);
  }

  // L = 1, M = 1: Y_1^1 = 0.5 * sqrt(3/(2*pi)) * sin(theta) * e^(i*phi) ~ 0.345494149 * sin(theta) * e^(i*phi)
  // (Condon-Shortley phase is cancelled)
  {
    double const theta = 0.8;
    double const phi = 0.6;
    auto val = SteinhardtCalculator::sphericalHarmonic(1, 1, {.theta = theta, .phi = phi});
    std::complex<double> expected = 0.345494149 * std::sin(theta) * std::polar(1.0, phi);
    EXPECT_NEAR(val.real(), expected.real(), 1e-8);
    EXPECT_NEAR(val.imag(), expected.imag(), 1e-8);
  }

  // L = 1, M = -1: Y_1^-1 = - (Y_1^1)*
  {
    double const theta = 0.8;
    double const phi = 0.6;
    auto val = SteinhardtCalculator::sphericalHarmonic(1, -1, {.theta = theta, .phi = phi});
    std::complex<double> expected = -std::conj(0.345494149 * std::sin(theta) * std::polar(1.0, phi));
    EXPECT_NEAR(val.real(), expected.real(), 1e-8);
    EXPECT_NEAR(val.imag(), expected.imag(), 1e-8);
  }
}

TEST_F(SteinhardtCalculatorTests, Wigner3j) {
  using correlation::calculators::SteinhardtCalculator;

  // Invalid selection where magnetic projections do not sum to 0
  EXPECT_DOUBLE_EQ(SteinhardtCalculator::wigner3j(1, 1, 1, 0, 0, 1), 0.0);

  // Invalid selection violating triangle inequality
  EXPECT_DOUBLE_EQ(SteinhardtCalculator::wigner3j(1, 1, 3, 0, 0, 0), 0.0);

  // Invalid selection with magnetic projection larger than angular momentum
  EXPECT_DOUBLE_EQ(SteinhardtCalculator::wigner3j(1, 1, 1, 2, 0, -2), 0.0);

  // Known analytical values
  // 3j(1, 1, 0, 0, 0, 0) = -1/sqrt(3) ~ -0.57735027
  EXPECT_NEAR(SteinhardtCalculator::wigner3j(1, 1, 0, 0, 0, 0), -1.0 / std::sqrt(3.0), 1e-8);

  // 3j(2, 2, 2, 0, 0, 0) = -sqrt(2/35) ~ -0.23904572
  EXPECT_NEAR(SteinhardtCalculator::wigner3j(2, 2, 2, 0, 0, 0), -std::sqrt(2.0 / 35.0), 1e-8);

  // 3j(2, 2, 1, 1, -1, 0) = -1/sqrt(30) ~ -0.18257419
  EXPECT_NEAR(SteinhardtCalculator::wigner3j(2, 2, 1, 1, -1, 0), -1.0 / std::sqrt(30.0), 1e-8);
}

} // namespace
} // namespace correlation::analysis
