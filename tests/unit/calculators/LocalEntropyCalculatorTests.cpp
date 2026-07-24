// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/StructureAnalyzer.hpp"
#include "calculators/LocalEntropyCalculator.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>
#include <random>
#include <vector>

namespace correlation::analysis {
namespace {

class LocalEntropyCalculatorTests : public ::testing::Test {
protected:
  static real_t getPeakEntropy(const Histogram &hist) {
    const auto &vals = hist.partials.at("Total");
    real_t max_val = -1.0;
    size_t max_bin = 0;
    for (size_t i = 0; i < vals.size(); ++i) {
      if (vals[i] > max_val) {
        max_val = vals[i];
        max_bin = i;
      }
    }
    return hist.bins[max_bin];
  }
};

TEST_F(LocalEntropyCalculatorTests, SimpleCubic) {
  // 1. Simple Cubic (SC) Structure
  // Box length 12.0 with SC lattice spacing 4.0 (27 atoms)
  correlation::core::Cell cell_sc({12.0, 12.0, 12.0, 90.0, 90.0, 90.0});
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        cell_sc.addAtom("Ar",
                        {static_cast<real_t>(i * 4.0), static_cast<real_t>(j * 4.0), static_cast<real_t>(k * 4.0)});
      }
    }
  }
  StructureAnalyzer const analyzer_sc(cell_sc, 6.5, {{6.5 * 6.5}}, false);
  auto hist_sc =
      correlation::calculators::LocalEntropyCalculator::calculate(cell_sc, &analyzer_sc, {.cutoff = 6.0, .sigma = 0.2});
  real_t const entropy_sc = getPeakEntropy(hist_sc);
  EXPECT_NEAR(entropy_sc, -6.15, 1e-4);
}

TEST_F(LocalEntropyCalculatorTests, BodyCenteredCubic) {
  // 2. Body-Centered Cubic (BCC) Structure
  // Box length 10.0 with BCC lattice spacing 5.0 (16 atoms)
  correlation::core::Cell cell_bcc({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        cell_bcc.addAtom("Ar",
                         {static_cast<real_t>(i * 5.0), static_cast<real_t>(j * 5.0), static_cast<real_t>(k * 5.0)});
        cell_bcc.addAtom("Ar", {i * 5.0 + 2.5, j * 5.0 + 2.5, k * 5.0 + 2.5});
      }
    }
  }
  StructureAnalyzer const analyzer_bcc(cell_bcc, 6.5, {{6.5 * 6.5}}, false);
  auto hist_bcc = correlation::calculators::LocalEntropyCalculator::calculate(cell_bcc, &analyzer_bcc,
                                                                              {.cutoff = 6.0, .sigma = 0.2});
  real_t const entropy_bcc = getPeakEntropy(hist_bcc);
  EXPECT_NEAR(entropy_bcc, -5.95, 1e-4);
}

TEST_F(LocalEntropyCalculatorTests, FaceCenteredCubic) {
  // 3. Face-Centered Cubic (FCC) Structure
  // Box length 10.0 with FCC lattice spacing 5.0 (32 atoms)
  correlation::core::Cell cell_fcc({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        cell_fcc.addAtom("Ar", {i * 5.0, j * 5.0, k * 5.0});
        cell_fcc.addAtom("Ar", {i * 5.0 + 2.5, j * 5.0 + 2.5, k * 5.0});
        cell_fcc.addAtom("Ar", {i * 5.0 + 2.5, j * 5.0, k * 5.0 + 2.5});
        cell_fcc.addAtom("Ar", {i * 5.0, j * 5.0 + 2.5, k * 5.0 + 2.5});
      }
    }
  }
  StructureAnalyzer const analyzer_fcc(cell_fcc, 6.5, {{6.5 * 6.5}}, false);
  auto hist_fcc = correlation::calculators::LocalEntropyCalculator::calculate(cell_fcc, &analyzer_fcc,
                                                                              {.cutoff = 6.0, .sigma = 0.2});
  real_t const entropy_fcc = getPeakEntropy(hist_fcc);
  EXPECT_NEAR(entropy_fcc, -8.95, 1e-4);
}

TEST_F(LocalEntropyCalculatorTests, Random) {
  // 4. Random / Gas Structure (Poisson distribution)
  // Box length 12.0 with 30 randomly placed atoms
  correlation::core::Cell cell_rand({12.0, 12.0, 12.0, 90.0, 90.0, 90.0});
  std::seed_seq seed{42};
  std::mt19937 gen{seed};

  std::uniform_real_distribution<double> dis(0.0, 12.0);
  for (int i = 0; i < 30; ++i) {
    cell_rand.addAtom("Ar",
                      {static_cast<real_t>(dis(gen)), static_cast<real_t>(dis(gen)), static_cast<real_t>(dis(gen))});
  }
  StructureAnalyzer const analyzer_rand(cell_rand, 6.5, {{6.5 * 6.5}}, false);
  auto hist_rand = correlation::calculators::LocalEntropyCalculator::calculate(cell_rand, &analyzer_rand,
                                                                               {.cutoff = 6.0, .sigma = 0.2});
  real_t const entropy_rand = getPeakEntropy(hist_rand);
  EXPECT_NEAR(entropy_rand, -2.05, correlation::is_single_precision ? 0.6 : 1e-4);
}

} // namespace
} // namespace correlation::analysis
