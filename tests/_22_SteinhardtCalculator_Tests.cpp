// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "Cell.hpp"
#include "StructureAnalyzer.hpp"
#include "calculators/SteinhardtCalculator.hpp"
#include <gtest/gtest.h>

class _22_SteinhardtCalculator_Tests : public ::testing::Test {
protected:
  void checkOutputs(const std::map<std::string, Histogram> &hists,
                    double expected_Q4, double expected_Q6,
                    double expected_W6_hat) {
    const auto &hist_Q4 = hists.at("Q4").partials.at("Total");
    const auto &hist_Q6 = hists.at("Q6").partials.at("Total");
    const auto &hist_W6 = hists.at("W6_hat").partials.at("Total");

    double q4_val = 0, q6_val = 0, w6_val = 0;

    // Find the non-zero bins
    size_t q4_bins = 100;
    double dQ = 1.0 / q4_bins;
    size_t w6_bins = 100;
    double W_min = -0.2;
    double W_max = 0.2;
    double dW = (W_max - W_min) / w6_bins;

    for (size_t b = 0; b < q4_bins; ++b) {
      if (hist_Q4[b] > 0)
        q4_val = (b + 0.5) * dQ;
      if (hist_Q6[b] > 0)
        q6_val = (b + 0.5) * dQ;
    }
    for (size_t b = 0; b < w6_bins; ++b) {
      if (hist_W6[b] > 0)
        w6_val = W_min + (b + 0.5) * dW;
    }

    EXPECT_NEAR(q4_val, expected_Q4, 0.015);
    EXPECT_NEAR(q6_val, expected_Q6, 0.015);
    EXPECT_NEAR(w6_val, expected_W6_hat, 0.015);
  }
};

TEST_F(_22_SteinhardtCalculator_Tests, SimpleCubic) {
  Cell cell({1.0, 1.0, 1.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {0.5, 0.5, 0.5});

  // ignore_periodic_self_interactions = false
  StructureAnalyzer analyzer(cell, 1.1, {{1.1 * 1.1}}, false);
  auto hists = SteinhardtCalculator::calculate(cell, &analyzer);

  checkOutputs(hists, 0.764, 0.354, 0.013);
}

TEST_F(_22_SteinhardtCalculator_Tests, BCC) {
  Cell cell({1.0, 1.0, 1.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {0.0, 0.0, 0.0});
  cell.addAtom("Ar", {0.5, 0.5, 0.5});

  StructureAnalyzer analyzer(cell, 1.1, {{1.1 * 1.1}},
                             false); // dist = sqrt(0.75) ~ 0.866 and 1.0
  auto hists = SteinhardtCalculator::calculate(cell, &analyzer);

  checkOutputs(hists, 0.036, 0.511, 0.013);
}

TEST_F(_22_SteinhardtCalculator_Tests, FCC) {
  Cell cell({1.0, 1.0, 1.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {0.0, 0.0, 0.0});
  cell.addAtom("Ar", {0.5, 0.5, 0.0});
  cell.addAtom("Ar", {0.5, 0.0, 0.5});
  cell.addAtom("Ar", {0.0, 0.5, 0.5});

  StructureAnalyzer analyzer(cell, 0.8, {{0.8 * 0.8}},
                             false); // dist = sqrt(0.5) ~ 0.707
  auto hists = SteinhardtCalculator::calculate(cell, &analyzer);

  checkOutputs(hists, 0.191, 0.575, -0.013); // W6 for FCC is approx -0.013
}

TEST_F(_22_SteinhardtCalculator_Tests, Icosahedral) {
  Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});

  double phi = (1.0 + std::sqrt(5.0)) / 2.0;
  double L = 1.0 / std::sqrt(1.0 + phi * phi); // Normalize dist to 1.0
  double pL = phi * L;

  cell.addAtom("Ar", {5.0, 5.0 + L, 5.0 + pL});
  cell.addAtom("Ar", {5.0, 5.0 + L, 5.0 - pL});
  cell.addAtom("Ar", {5.0, 5.0 - L, 5.0 + pL});
  cell.addAtom("Ar", {5.0, 5.0 - L, 5.0 - pL});

  cell.addAtom("Ar", {5.0 + L, 5.0 + pL, 5.0});
  cell.addAtom("Ar", {5.0 + L, 5.0 - pL, 5.0});
  cell.addAtom("Ar", {5.0 - L, 5.0 + pL, 5.0});
  cell.addAtom("Ar", {5.0 - L, 5.0 - pL, 5.0});

  cell.addAtom("Ar", {5.0 + pL, 5.0, 5.0 + L});
  cell.addAtom("Ar", {5.0 + pL, 5.0, 5.0 - L});
  cell.addAtom("Ar", {5.0 - pL, 5.0, 5.0 + L});
  cell.addAtom("Ar", {5.0 - pL, 5.0, 5.0 - L});

  // Edge length is ~1.05. Using cutoff 1.02 ensures surface atoms only see
  // center. Thus they will have 1 neighbor, Ql=1.0, and be excluded from
  // histogram!
  StructureAnalyzer analyzer(cell, 1.02, {{1.02 * 1.02}}, true);
  auto hists = SteinhardtCalculator::calculate(cell, &analyzer);

  checkOutputs(hists, 0.000, 0.663,
               -0.169); // W6_hat for Icosahedral is approx -0.1697
}
