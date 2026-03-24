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
  const auto &hist = df.getHistogram("Debye_S(Q)");
  
  EXPECT_EQ(hist.bins.size(), 10);
  EXPECT_EQ(hist.bin_label, "Q (Å⁻¹)");

  // Check partials
  EXPECT_TRUE(hist.partials.count("C-C"));
  EXPECT_TRUE(hist.partials.count("C-O") || hist.partials.count("O-C"));
  EXPECT_TRUE(hist.partials.count("O-O"));

  std::string cross_key = hist.partials.count("C-O") ? "C-O" : "O-C";

  for (size_t i = 0; i < hist.bins.size(); ++i) {
    double Q = hist.bins[i];
    
    // For single atoms, S(Q) = 1.0
    EXPECT_NEAR(hist.partials.at("C-C")[i], 1.0, 1e-6);
    EXPECT_NEAR(hist.partials.at("O-O")[i], 1.0, 1e-6);
    
    // For cross term, S(Q) = sin(Q * r) / (Q * r)
    double expected_cross = std::sin(Q * 1.2) / (Q * 1.2);
    EXPECT_NEAR(hist.partials.at(cross_key)[i], expected_cross, 1e-5);
  }
}
