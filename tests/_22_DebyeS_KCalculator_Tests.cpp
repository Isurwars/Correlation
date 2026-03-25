// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>

#include "../include/Cell.hpp"
#include "../include/DistributionFunctions.hpp"
#include "../include/calculators/DebyeS_KCalculator.hpp"

// Test fixture for Debye S(k) Calculator
class _22_DebyeS_KCalculator_Tests : public ::testing::Test {
protected:
  void SetUp() override {}
};

TEST_F(_22_DebyeS_KCalculator_Tests, CalculatesDiatomicMolecule) {
  // Arrange
  Cell cell;
  // A diatomic molecule CO with bond length 1.2 Å
  cell.addAtom("C", {0.0, 0.0, 0.0});
  cell.addAtom("O", {1.2, 0.0, 0.0});

  DistributionFunctions df(cell);

  AnalysisSettings settings;
  settings.q_max = 10.0;
  settings.q_bin_width = 1.0;

  DebyeS_KCalculator calc;

  // Act
  calc.calculateFrame(df, settings);

  // Assert
  const auto &hist = df.getHistogram("debye_S_q");
  
  EXPECT_EQ(hist.bins.size(), 10);
  EXPECT_EQ(hist.bin_label, "Q (Å⁻¹)");

  // Check partials
  EXPECT_TRUE(hist.partials.count("C-C"));
  EXPECT_TRUE(hist.partials.count("C-O") || hist.partials.count("O-C"));
  EXPECT_TRUE(hist.partials.count("O-O"));

  std::string cross_key = "C-O";

  for (size_t i = 0; i < hist.bins.size(); ++i) {
    double Q = hist.bins[i];
    
    // For single atoms, S(Q) = 1.0 (self term) + 0 (distinct terms because N_A=1)
    EXPECT_NEAR(hist.partials.at("C-C")[i], 1.0, 1e-6);
    EXPECT_NEAR(hist.partials.at("O-O")[i], 1.0, 1e-6);
    
    // For cross term, S(Q) = (1/sqrt(Na*Nb)) * sum_sinc = 1 * sinc = sinc
    double expected_cross = std::sin(Q * 1.2) / (Q * 1.2);
    EXPECT_NEAR(hist.partials.at(cross_key)[i], expected_cross, 1e-5);
  }
}

TEST_F(_22_DebyeS_KCalculator_Tests, VerifiesSZeroLimit) {
  // Arrange
  Cell cell;
  // 4 identical atoms
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {2.0, 0.0, 0.0});
  cell.addAtom("Si", {0.0, 2.0, 0.0});
  cell.addAtom("Si", {0.0, 0.0, 2.0});

  DistributionFunctions df(cell);
  DebyeS_KCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 0.02;
  settings.q_bin_width = 0.02; // first bin at Q=0.01

  // Act
  calc.calculateFrame(df, settings);
  const auto &hist = df.getHistogram("debye_S_q");

  // Assert
  // S(0) should be N = 4.0. At Q=0.01, it should be very close.
  // sinc(0.01 * 2.0) = sin(0.02)/0.02 = 0.999933
  // S(Q) = 1 + (2/4) * 6 * 0.999933 = 1 + 3 * 0.9999 = 3.9998
  EXPECT_NEAR(hist.partials.at("Total")[0], 4.0, 0.01);
}

TEST_F(_22_DebyeS_KCalculator_Tests, VerifiesSZeroLimitMultiComponent) {
  // Arrange
  Cell cell;
  // 1 C and 1 O very close
  cell.addAtom("C", {0.0, 0.0, 0.0});
  cell.addAtom("O", {0.001, 0.0, 0.0});

  DistributionFunctions df(cell);
  DebyeS_KCalculator calc;
  AnalysisSettings settings;
  settings.q_max = 0.01;
  settings.q_bin_width = 0.01; // first bin at Q=0.005

  // Act
  calc.calculateFrame(df, settings);
  const auto &hist = df.getHistogram("debye_S_q");

  // Assert
  // S(0) should be N = 2.0. 
  // wCC=0.25, wOO=0.25, wCO=0.5
  // SCC=1, SOO=1, SCO=sinc(0.005*0.001) approx 1
  // Stotal = 0.25*1 + 0.25*1 + 0.5*1 = 1.0? 
  // Wait, my manual calculation before said 2. Let's re-verify.
  // cA=0.5, cB=0.5. 
  // S_total = cA*SAA + cB*SBB + 2*sqrt(cA*cB)*SAB
  // S_total = 0.5*1 + 0.5*1 + 2*0.5*1 = 2.0. YES!
  EXPECT_NEAR(hist.partials.at("Total")[0], 2.0, 0.01);
}
