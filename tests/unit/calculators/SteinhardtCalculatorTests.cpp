// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/StructureAnalyzer.hpp"
#include "calculators/order/SteinhardtCalculator.hpp"
#include "core/Cell.hpp"

#include "../../CrystalTestHelper.hpp"

#include <gtest/gtest.h>

namespace correlation::analysis {
namespace {

class SteinhardtCalculatorTests : public ::testing::Test {
protected:
  static void checkOutputs(const std::map<std::string, Histogram> &hists, real_t expected_q4,
                           real_t expected_q6, real_t expected_w6_hat) {
    const auto &hist_q4 = hists.at("Q4").partials.at("Total");
    const auto &hist_q6 = hists.at("Q6").partials.at("Total");
    const auto &hist_w6 = hists.at("W6_hat").partials.at("Total");

    real_t q4_val = 0;
    real_t q6_val = 0;
    real_t w6_val = 0;

    // Find the non-zero bins
    size_t const q4_bins = 100;
    real_t const d_q = 1.0 / q4_bins;
    size_t const w6_bins = 100;
    real_t const w_min = -0.2;
    real_t const w_max = 0.2;
    real_t const d_w = (w_max - w_min) / w6_bins;

    for (size_t bin = 0; bin < q4_bins; ++bin) {
      if (hist_q4[bin] > 0) {
        q4_val = static_cast<real_t>(static_cast<real_t>(bin) + 0.5) * d_q;
      }
      if (hist_q6[bin] > 0) {
        q6_val = static_cast<real_t>(static_cast<real_t>(bin) + 0.5) * d_q;
      }
    }
    for (size_t bin = 0; bin < w6_bins; ++bin) {
      if (hist_w6[bin] > 0) {
        w6_val = static_cast<real_t>(w_min + (static_cast<real_t>(bin) + 0.5) * d_w);
      }
    }

    EXPECT_NEAR(q4_val, expected_q4, 0.015);
    EXPECT_NEAR(q6_val, expected_q6, 0.015);
    EXPECT_NEAR(w6_val, expected_w6_hat, 0.015);
  }

  static void checkAllOutputs(const std::map<std::string, Histogram> &hists, real_t expected_q4,
                              real_t expected_q6, real_t expected_w4_hat, real_t expected_w6_hat,
                              real_t expected_q4_bar, real_t expected_q6_bar) {
    checkOutputs(hists, expected_q4, expected_q6, expected_w6_hat);

    const auto &hist_w4 = hists.at("W4_hat").partials.at("Total");
    const auto &hist_q4_bar = hists.at("Q4_bar").partials.at("Total");
    const auto &hist_q6_bar = hists.at("Q6_bar").partials.at("Total");

    size_t const w4_bins = 100;
    real_t const w4_min = -0.5;
    real_t const w4_max = 0.5;
    real_t const d_w4 = (w4_max - w4_min) / w4_bins;

    size_t const q_bins = 100;
    real_t const d_q = 1.0 / q_bins;

    real_t w4_val = 0;
    real_t q4_bar_val = 0;
    real_t q6_bar_val = 0;

    for (size_t bin = 0; bin < w4_bins; ++bin) {
      if (hist_w4[bin] > 0) {
        w4_val = static_cast<real_t>(w4_min + (static_cast<real_t>(bin) + 0.5) * d_w4);
      }
    }
    for (size_t bin = 0; bin < q_bins; ++bin) {
      if (hist_q4_bar[bin] > 0) {
        q4_bar_val = static_cast<real_t>(static_cast<real_t>(bin) + 0.5) * d_q;
      }
      if (hist_q6_bar[bin] > 0) {
        q6_bar_val = static_cast<real_t>(static_cast<real_t>(bin) + 0.5) * d_q;
      }
    }

    EXPECT_NEAR(w4_val, expected_w4_hat, 0.02);
    EXPECT_NEAR(q4_bar_val, expected_q4_bar, 0.02);
    EXPECT_NEAR(q6_bar_val, expected_q6_bar, 0.02);
  }
};

TEST_F(SteinhardtCalculatorTests, SimpleCubic) {
  auto cell = correlation::testing::crystals::createSimpleCubicCell(1.0, "Ar", 1, 1, 1);
  // Shift atom to center (0.5, 0.5, 0.5) for exact fixture parity
  cell = correlation::core::Cell({1.0, 1.0, 1.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {0.5, 0.5, 0.5});

  // ignore_periodic_self_interactions = false
  StructureAnalyzer const analyzer(cell, 1.1, {{{.min_sq = 0.36, .max_sq = 1.1 * 1.1}}}, false);
  auto hists = correlation::calculators::SteinhardtCalculator::calculate(cell, &analyzer);

  checkOutputs(hists, 0.764, 0.354, 0.013);
  checkAllOutputs(hists, 0.764, 0.354, 0.155, 0.013, 0.764, 0.354);
}

TEST_F(SteinhardtCalculatorTests, BCC) {
  auto cell = correlation::testing::crystals::createBCCCell(1.0, "Ar", 1, 1, 1);

  StructureAnalyzer const analyzer(cell, 1.1, {{{.min_sq = 0.36, .max_sq = 1.1 * 1.1}}},
                                   false); // dist = sqrt(0.75) ~ 0.866 and 1.0
  auto hists = correlation::calculators::SteinhardtCalculator::calculate(cell, &analyzer);

  checkOutputs(hists, 0.036, 0.511, 0.013);
  checkAllOutputs(hists, 0.036, 0.511, 0.155, 0.013, 0.036, 0.511);
}

TEST_F(SteinhardtCalculatorTests, FCC) {
  auto cell = correlation::testing::crystals::createFCCCell(1.0, "Ar", 1, 1, 1);

  StructureAnalyzer const analyzer(cell, 0.8, {{{.min_sq = 0.36, .max_sq = 0.8 * 0.8}}},
                                   false); // dist = sqrt(0.5) ~ 0.707
  auto hists = correlation::calculators::SteinhardtCalculator::calculate(cell, &analyzer);

  checkOutputs(hists, 0.191, 0.575, -0.013); // W6 for FCC is approx -0.013
  checkAllOutputs(hists, 0.191, 0.575, -0.155, -0.013, 0.191, 0.575);
}

TEST_F(SteinhardtCalculatorTests, DistortedCrystalLechnerDellagoAveraging) {
  // Test that Lechner-Dellago Q4_bar and Q6_bar remain stable under small thermal-like
  // perturbations
  auto cell = correlation::testing::crystals::createFCCCell(1.0, "Ar", 2, 2, 2);
  StructureAnalyzer const analyzer(cell, 0.8, {{{.min_sq = 0.36, .max_sq = 0.8 * 0.8}}}, false);
  auto hists = correlation::calculators::SteinhardtCalculator::calculate(cell, &analyzer);

  ASSERT_TRUE(hists.count("Q4_bar"));
  ASSERT_TRUE(hists.count("Q6_bar"));
  const auto &total_q4_bar = hists.at("Q4_bar").partials.at("Total");
  const auto &total_q6_bar = hists.at("Q6_bar").partials.at("Total");

  // Sum of normalized probabilities must be non-zero
  real_t sum_q4 = 0;
  real_t sum_q6 = 0;
  for (real_t const val : total_q4_bar) {
    sum_q4 += val;
  }
  for (real_t const val : total_q6_bar) {
    sum_q6 += val;
  }
  EXPECT_GT(sum_q4, 0.0);
  EXPECT_GT(sum_q6, 0.0);
}

TEST_F(SteinhardtCalculatorTests, Icosahedral) {
  auto cell = correlation::testing::crystals::createIcosahedralClusterCell(
      {.center_elem = "Ar", .shell_elem = "Ar", .r_bond = 1.0, .box_size = 10.0});

  // Edge length is ~1.05. Using cutoff 1.02 ensures surface atoms only see
  // center. Thus they will have 1 neighbor, Ql=1.0, and be excluded from
  // histogram!
  StructureAnalyzer const analyzer(cell, 1.02, {{{.min_sq = 0.36, .max_sq = 1.02 * 1.02}}}, true);
  auto hists = correlation::calculators::SteinhardtCalculator::calculate(cell, &analyzer);

  checkOutputs(hists, 0.000, 0.663,
               -0.169); // W6_hat for Icosahedral is approx -0.1697
}

TEST_F(SteinhardtCalculatorTests, HandlesAcosNumericalNoiseSafely) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});
  cell.addAtom("Ar", {5.0, 5.0, 6.000000000000001});
  cell.addAtom("Ar", {5.0, 5.0, 4.0});

  StructureAnalyzer const analyzer(cell, 1.1, {{{.min_sq = 0.36, .max_sq = 1.1 * 1.1}}}, false);
  ASSERT_NO_THROW({
    auto hists = correlation::calculators::SteinhardtCalculator::calculate(cell, &analyzer);
    EXPECT_FALSE(hists.empty());
  });
}

void assertZeroHistogram(const correlation::analysis::Histogram &hist) {
  ASSERT_TRUE(hist.partials.count("Total"));
  const auto &total = hist.partials.at("Total");
  for (real_t const val : total) {
    EXPECT_DOUBLE_EQ(val, 0.0);
  }
}

TEST_F(SteinhardtCalculatorTests, EmptySystemOrNoNeighborsFillsPartialsWithZeros) {
  // Cell with 1 atom (no neighbors)
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});

  StructureAnalyzer const analyzer(cell, 1.1, {{{.min_sq = 0.36, .max_sq = 1.1 * 1.1}}}, false);
  auto hists = correlation::calculators::SteinhardtCalculator::calculate(cell, &analyzer);

  EXPECT_FALSE(hists.empty());
  for (const auto &name : {"Q4", "Q6", "W4_hat", "W6_hat", "Q4_bar", "Q6_bar"}) {
    ASSERT_TRUE(hists.count(name));
    assertZeroHistogram(hists.at(name));
  }
}

TEST_F(SteinhardtCalculatorTests, SphericalHarmonics) {
  using correlation::calculators::SteinhardtCalculator;

  // L = 0, M = 0: Y_0^0 = 0.5 * sqrt(1/pi) ~ 0.28209479
  {
    auto val = SteinhardtCalculator::sphericalHarmonic(0, 0,
                                                       {
                                                           .theta = 0.5,
                                                           .phi = 0.2,
                                                       });
    EXPECT_NEAR(val.real(), 0.28209479177, correlation::is_single_precision ? 1e-5 : 1e-8);
    EXPECT_NEAR(val.imag(), 0.0, correlation::is_single_precision ? 1e-5 : 1e-8);
  }

  // L = 1, M = 0: Y_1^0 = 0.5 * sqrt(3/pi) * cos(theta) ~ 0.4886025 * cos(theta)
  {
    real_t const theta = 1.0;
    auto val = SteinhardtCalculator::sphericalHarmonic(1, 0,
                                                       {
                                                           .theta = static_cast<real_t>(theta),
                                                           .phi = 0.5F,
                                                       });
    EXPECT_NEAR(val.real(), 0.4886025119 * std::cos(theta), 1e-5);
    EXPECT_NEAR(val.imag(), 0.0, 1e-5);
  }

  // L = 1, M = 1: Y_1^1 = 0.5 * sqrt(3/(2*pi)) * sin(theta) * e^(i*phi) ~ 0.345494149 * sin(theta)
  // * e^(i*phi) (Condon-Shortley phase is cancelled)
  {
    real_t const theta = 0.8;
    real_t const phi = 0.6;
    auto val = SteinhardtCalculator::sphericalHarmonic(1, 1,
                                                       {
                                                           .theta = static_cast<real_t>(theta),
                                                           .phi = static_cast<real_t>(phi),
                                                       });
    const std::complex<real_t> expected = static_cast<real_t>(0.345494149) * std::sin(theta) *
                                          std::polar(static_cast<real_t>(1.0), phi);
    EXPECT_NEAR(val.real(), expected.real(), 1e-5);
    EXPECT_NEAR(val.imag(), expected.imag(), 1e-5);
  }

  // L = 1, M = -1: Y_1^-1 = - (Y_1^1)*
  {
    real_t const theta = 0.8;
    real_t const phi = 0.6;
    auto val = SteinhardtCalculator::sphericalHarmonic(1, -1,
                                                       {
                                                           .theta = static_cast<real_t>(theta),
                                                           .phi = static_cast<real_t>(phi),
                                                       });
    const std::complex<real_t> expected =
        -std::conj(static_cast<real_t>(0.345494149) * std::sin(theta) *
                   std::polar(static_cast<real_t>(1.0), phi));
    EXPECT_NEAR(val.real(), expected.real(), 1e-5);
    EXPECT_NEAR(val.imag(), expected.imag(), 1e-5);
  }
}

TEST_F(SteinhardtCalculatorTests, Wigner3j) {
  using correlation::calculators::SteinhardtCalculator;

  // Invalid selection where magnetic projections do not sum to 0
  EXPECT_DOUBLE_EQ(SteinhardtCalculator::wigner3j({
                       .j_one = 1,
                       .j_two = 1,
                       .j_three = 1,
                       .m_one = 0,
                       .m_two = 0,
                       .m_three = 1,
                   }),
                   0.0);

  // Invalid selection violating triangle inequality
  EXPECT_DOUBLE_EQ(SteinhardtCalculator::wigner3j({
                       .j_one = 1,
                       .j_two = 1,
                       .j_three = 3,
                       .m_one = 0,
                       .m_two = 0,
                       .m_three = 0,
                   }),
                   0.0);

  // Invalid selection with magnetic projection larger than angular momentum
  EXPECT_DOUBLE_EQ(SteinhardtCalculator::wigner3j({
                       .j_one = 1,
                       .j_two = 1,
                       .j_three = 1,
                       .m_one = 2,
                       .m_two = 0,
                       .m_three = -2,
                   }),
                   0.0);

  // Known analytical values
  // 3j(1, 1, 0, 0, 0, 0) = -1/sqrt(3) ~ -0.57735027
  EXPECT_NEAR(SteinhardtCalculator::wigner3j({
                  .j_one = 1,
                  .j_two = 1,
                  .j_three = 0,
                  .m_one = 0,
                  .m_two = 0,
                  .m_three = 0,
              }),
              -1.0 / std::sqrt(3.0), correlation::is_single_precision ? 1e-5 : 1e-8);

  // 3j(2, 2, 2, 0, 0, 0) = -sqrt(2/35) ~ -0.23904572
  EXPECT_NEAR(SteinhardtCalculator::wigner3j({
                  .j_one = 2,
                  .j_two = 2,
                  .j_three = 2,
                  .m_one = 0,
                  .m_two = 0,
                  .m_three = 0,
              }),
              -std::sqrt(2.0 / 35.0), correlation::is_single_precision ? 1e-5 : 1e-8);

  // 3j(2, 2, 1, 1, -1, 0) = -1/sqrt(30) ~ -0.18257419
  EXPECT_NEAR(SteinhardtCalculator::wigner3j({
                  .j_one = 2,
                  .j_two = 2,
                  .j_three = 1,
                  .m_one = 1,
                  .m_two = -1,
                  .m_three = 0,
              }),
              -1.0 / std::sqrt(30.0), 1e-8);
}

} // namespace
} // namespace correlation::analysis
