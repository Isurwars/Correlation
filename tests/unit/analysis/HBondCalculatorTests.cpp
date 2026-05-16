// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/StructureAnalyzer.hpp"
#include "calculators/HBondCalculator.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>

namespace correlation::calculators {

// ============================================================================
// Null / Empty Input
// ============================================================================

/// Passing a null neighbor pointer should return an empty histogram.
TEST(HBondCalculatorTests, NullNeighborsReturnsEmptyHistogram) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("O", {5.0, 5.0, 5.0});

  auto hist = HBondCalculator::calculate(cell, nullptr);

  EXPECT_TRUE(hist.bins.empty());
  EXPECT_TRUE(hist.partials.empty());
}

// ============================================================================
// No Electronegative Atoms
// ============================================================================

/// A system with no electronegative atoms (O, N, F, S) produces no H-bonds.
/// Note: the calculator still produces a histogram with bin {0} when there are
/// no donors. We verify that the entire distribution is at N_HB=0.
TEST(HBondCalculatorTests, NoElectronegativeAtomsNoHBonds) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("C", {5.0, 5.0, 5.0});
  cell.addAtom("H", {5.0, 5.0, 6.0});
  cell.addAtom("H", {5.0, 6.0, 5.0});

  double cutoff = 5.0;
  std::vector<std::vector<double>> bcsq(2, std::vector<double>(2, cutoff * cutoff));
  analysis::StructureAnalyzer analyzer(cell, cutoff, bcsq, true);
  auto hist = HBondCalculator::calculate(cell, &analyzer);

  // No electronegative atoms means no donor/acceptor analysis.
  // The histogram may be empty or may have bin {0} with zero frequency.
  // Either way, no non-zero H-bond counts should exist.
  for (size_t i = 0; i < hist.bins.size(); ++i) {
    if (hist.bins[i] > 0) {
      EXPECT_NEAR(hist.partials.at("Total")[i], 0.0, 1e-12)
          << "Non-zero H-bond count should have zero probability.";
    }
  }
}

/// A system with only hydrogens — no donors or acceptors.
TEST(HBondCalculatorTests, OnlyHydrogensNoHBonds) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("H", {5.0, 5.0, 5.0});
  cell.addAtom("H", {5.0, 5.0, 6.0});

  double cutoff = 5.0;
  std::vector<std::vector<double>> bcsq(1, std::vector<double>(1, cutoff * cutoff));
  analysis::StructureAnalyzer analyzer(cell, cutoff, bcsq, true);
  auto hist = HBondCalculator::calculate(cell, &analyzer);

  // No electronegative atoms → no H-bond counting, but histogram may
  // still contain bin {0}. Verify no non-zero H-bond counts.
  for (size_t i = 0; i < hist.bins.size(); ++i) {
    if (hist.bins[i] > 0) {
      EXPECT_NEAR(hist.partials.at("Total")[i], 0.0, 1e-12);
    }
  }
}

// ============================================================================
// Histogram Consistency
// ============================================================================

/// Bins and all partial vectors must have matching sizes.
TEST(HBondCalculatorTests, HistogramDimensionsAreConsistent) {
  // Simple water-like: O-H...O
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  // Donor water
  cell.addAtom("O", {10.0, 10.0, 10.0});
  cell.addAtom("H", {10.0, 10.0, 10.96}); // O-H ≈ 0.96 Å
  // Acceptor oxygen
  cell.addAtom("O", {10.0, 10.0, 12.5}); // D-A distance ≈ 2.5 Å < 3.5 Å

  double cutoff = 5.0;
  // 2 element types: O (0), H (1)
  std::vector<std::vector<double>> bcsq(2, std::vector<double>(2, cutoff * cutoff));
  analysis::StructureAnalyzer analyzer(cell, cutoff, bcsq, true);
  auto hist = HBondCalculator::calculate(cell, &analyzer);

  for (const auto &[key, values] : hist.partials) {
    EXPECT_EQ(values.size(), hist.bins.size())
        << "Partial '" << key << "' size mismatch with bins.";
  }
}

/// The distribution must sum to 1.0 (normalized over electronegative atoms).
TEST(HBondCalculatorTests, TotalDistributionSumsToOne) {
  // 3 water molecules forming H-bonds
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  // Water 1: donor
  cell.addAtom("O", {10.0, 10.0, 10.0});
  cell.addAtom("H", {10.0, 10.0, 10.96});
  cell.addAtom("H", {10.0, 10.96, 10.0});
  // Water 2: acceptor
  cell.addAtom("O", {10.0, 10.0, 12.5});
  cell.addAtom("H", {10.0, 10.0, 13.46});
  cell.addAtom("H", {10.0, 10.96, 12.5});
  // Water 3: far away, no H-bond
  cell.addAtom("O", {18.0, 18.0, 18.0});
  cell.addAtom("H", {18.0, 18.0, 18.96});
  cell.addAtom("H", {18.0, 18.96, 18.0});

  double cutoff = 5.0;
  // 2 elements: O(0), H(1)
  std::vector<std::vector<double>> bcsq(2, std::vector<double>(2, cutoff * cutoff));
  analysis::StructureAnalyzer analyzer(cell, cutoff, bcsq, true);
  auto hist = HBondCalculator::calculate(cell, &analyzer);

  ASSERT_TRUE(hist.partials.count("Total"));
  double sum = 0.0;
  for (double v : hist.partials.at("Total")) {
    sum += v;
  }
  EXPECT_NEAR(sum, 1.0, 1e-12);
}

/// Histogram metadata should be correctly set.
TEST(HBondCalculatorTests, HistogramMetadataIsPopulated) {
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  cell.addAtom("O", {10.0, 10.0, 10.0});
  cell.addAtom("H", {10.0, 10.0, 10.96});
  cell.addAtom("O", {10.0, 10.0, 12.5});

  double cutoff = 5.0;
  std::vector<std::vector<double>> bcsq(2, std::vector<double>(2, cutoff * cutoff));
  analysis::StructureAnalyzer analyzer(cell, cutoff, bcsq, true);
  auto hist = HBondCalculator::calculate(cell, &analyzer);

  EXPECT_EQ(hist.title, "Hydrogen Bond Distribution");
  EXPECT_EQ(hist.x_label, "N_HB");
  EXPECT_EQ(hist.y_label, "P(N_HB)");
  EXPECT_EQ(hist.file_suffix, "_HBond");
  EXPECT_FALSE(hist.description.empty());
}

// ============================================================================
// Known H-Bond Geometry
// ============================================================================

/// A textbook O-H...O hydrogen bond: D-A < 3.5 Å, angle H-D...A < 30°.
/// Place atoms in a straight line: O(donor) — H — O(acceptor).
/// Angle H-D...A = 0° (perfectly aligned), distance D-A = 2.5 Å.
TEST(HBondCalculatorTests, LinearHBond_IsDetected) {
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  // O donor at origin
  cell.addAtom("O", {10.0, 10.0, 10.0});
  // H along z-axis, O-H = 0.96 Å
  cell.addAtom("H", {10.0, 10.0, 10.96});
  // O acceptor along z-axis, D-A = 2.5 Å
  cell.addAtom("O", {10.0, 10.0, 12.5});

  double cutoff = 5.0;
  std::vector<std::vector<double>> bcsq(2, std::vector<double>(2, cutoff * cutoff));
  analysis::StructureAnalyzer analyzer(cell, cutoff, bcsq, true);
  auto hist = HBondCalculator::calculate(cell, &analyzer);

  // There should be at least one H-bond detected
  ASSERT_FALSE(hist.bins.empty()) << "Expected at least one H-bond bin.";

  // Check that some atoms have non-zero H-bond counts
  bool found_nonzero = false;
  for (const auto &v : hist.partials.at("Total")) {
    if (v > 0 && hist.bins[&v - &hist.partials.at("Total")[0]] > 0) {
      found_nonzero = true;
      break;
    }
  }
  EXPECT_TRUE(found_nonzero)
      << "Linear O-H...O bond should produce non-zero H-bond counts.";
}

/// An O-H...O geometry where the H points away from the acceptor.
/// By placing H on the opposite side of the donor from the acceptor AND
/// ensuring H is far from the acceptor (> bond cutoff), no H-bond forms.
TEST(HBondCalculatorTests, OppositeDirection_NotDetected) {
  correlation::core::Cell cell({40.0, 40.0, 40.0, 90.0, 90.0, 90.0});
  // O donor at center
  cell.addAtom("O", {20.0, 20.0, 20.0});
  // H pointing AWAY from the acceptor, O-H = 0.96 Å along -z
  cell.addAtom("H", {20.0, 20.0, 19.04});
  // O acceptor along +z direction, D-A = 3.0 Å (within 3.5 Å H-bond cutoff)
  // H-to-acceptor distance = 3.96 Å (> typical bond cutoff of ~1.5 Å)
  cell.addAtom("O", {20.0, 20.0, 23.0});

  // Use a small bond cutoff (1.5 Å) so that H is only bonded to the donor O,
  // not to the acceptor O.
  double cutoff = 1.5;
  // 2 elements: O(0), H(1)
  std::vector<std::vector<double>> bcsq(2, std::vector<double>(2, cutoff * cutoff));
  analysis::StructureAnalyzer analyzer(cell, cutoff, bcsq, true);
  auto hist = HBondCalculator::calculate(cell, &analyzer);

  // H points opposite to acceptor → angle ≈ 180° >> 30°, no H-bond.
  // All electronegative atoms should have 0 H-bonds.
  for (size_t i = 0; i < hist.bins.size(); ++i) {
    if (hist.bins[i] > 0) {
      EXPECT_NEAR(hist.partials.at("Total")[i], 0.0, 1e-12)
          << "Opposite-direction geometry should yield zero H-bonds.";
    }
  }
}

/// D-A distance > 3.5 Å: even with perfect alignment, no H-bond.
TEST(HBondCalculatorTests, TooFarApart_NotDetected) {
  correlation::core::Cell cell({30.0, 30.0, 30.0, 90.0, 90.0, 90.0});
  // O donor
  cell.addAtom("O", {10.0, 10.0, 10.0});
  // H along z
  cell.addAtom("H", {10.0, 10.0, 10.96});
  // O acceptor far along z (D-A = 5.0 Å > 3.5 Å cutoff)
  cell.addAtom("O", {10.0, 10.0, 15.0});

  double cutoff = 6.0;
  std::vector<std::vector<double>> bcsq(2, std::vector<double>(2, cutoff * cutoff));
  analysis::StructureAnalyzer analyzer(cell, cutoff, bcsq, true);
  auto hist = HBondCalculator::calculate(cell, &analyzer);

  // All electronegative atoms should have 0 H-bonds
  if (!hist.bins.empty()) {
    EXPECT_NEAR(hist.partials.at("Total")[0], 1.0, 1e-12)
        << "D-A > 3.5 Å should yield zero H-bonds.";
  }
}

// ============================================================================
// Electronegative Element Coverage
// ============================================================================

/// N-H...O hydrogen bond should also be detected.
TEST(HBondCalculatorTests, NitrogenDonor_IsDetected) {
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  // N donor
  cell.addAtom("N", {10.0, 10.0, 10.0});
  // H along z-axis, N-H ≈ 1.01 Å
  cell.addAtom("H", {10.0, 10.0, 11.01});
  // O acceptor, D-A ≈ 2.8 Å
  cell.addAtom("O", {10.0, 10.0, 12.8});

  double cutoff = 5.0;
  // 3 elements: N(0), H(1), O(2)
  std::vector<std::vector<double>> bcsq(3, std::vector<double>(3, cutoff * cutoff));
  analysis::StructureAnalyzer analyzer(cell, cutoff, bcsq, true);
  auto hist = HBondCalculator::calculate(cell, &analyzer);

  ASSERT_FALSE(hist.bins.empty());

  // N and O should each have at least 1 H-bond, so bins > 0 should be present
  bool found_hbond = false;
  for (size_t i = 0; i < hist.bins.size(); ++i) {
    if (hist.bins[i] > 0 && hist.partials.at("Total")[i] > 0) {
      found_hbond = true;
      break;
    }
  }
  EXPECT_TRUE(found_hbond) << "N-H...O bond should be detected.";
}

/// F-H...F hydrogen bond (fluorine donor and acceptor).
TEST(HBondCalculatorTests, FluorineDonorAndAcceptor) {
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  // F donor
  cell.addAtom("F", {10.0, 10.0, 10.0});
  // H along z-axis, F-H ≈ 0.92 Å
  cell.addAtom("H", {10.0, 10.0, 10.92});
  // F acceptor, D-A ≈ 2.5 Å
  cell.addAtom("F", {10.0, 10.0, 12.5});

  double cutoff = 5.0;
  std::vector<std::vector<double>> bcsq(2, std::vector<double>(2, cutoff * cutoff));
  analysis::StructureAnalyzer analyzer(cell, cutoff, bcsq, true);
  auto hist = HBondCalculator::calculate(cell, &analyzer);

  ASSERT_FALSE(hist.bins.empty());

  bool found_hbond = false;
  for (size_t i = 0; i < hist.bins.size(); ++i) {
    if (hist.bins[i] > 0 && hist.partials.at("Total")[i] > 0) {
      found_hbond = true;
      break;
    }
  }
  EXPECT_TRUE(found_hbond) << "F-H...F bond should be detected.";
}

// ============================================================================
// Extreme / Edge Cases
// ============================================================================

/// H and donor at the same position (degenerate). Should not crash (NaN guard).
TEST(HBondCalculatorTests, CoincidentHAndDonor_DoesNotCrash) {
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  // O and H at exactly the same position
  cell.addAtom("O", {10.0, 10.0, 10.0});
  cell.addAtom("H", {10.0, 10.0, 10.0}); // Degenerate!
  cell.addAtom("O", {10.0, 10.0, 12.5});

  double cutoff = 5.0;
  std::vector<std::vector<double>> bcsq(2, std::vector<double>(2, cutoff * cutoff));
  analysis::StructureAnalyzer analyzer(cell, cutoff, bcsq, true);

  // Should not throw or produce NaN
  ASSERT_NO_THROW({
    auto hist = HBondCalculator::calculate(cell, &analyzer);

    // Verify no NaN values in the output
    for (const auto &[key, values] : hist.partials) {
      for (double v : values) {
        EXPECT_FALSE(std::isnan(v))
            << "NaN detected in partial '" << key << "'";
      }
    }
  });
}

/// Electronegative atom with no bonded hydrogens is just an acceptor (0 H-bonds
/// as donor). It may still appear with 0 or more H-bonds if it's an acceptor.
TEST(HBondCalculatorTests, AcceptorOnly_NoHydrogens) {
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  // O with no bonded H — can only be acceptor, not donor
  cell.addAtom("O", {10.0, 10.0, 10.0});
  // Another O with bonded H (donor)
  cell.addAtom("O", {10.0, 10.0, 12.5});
  cell.addAtom("H", {10.0, 10.0, 11.54}); // H between the two O atoms

  double cutoff = 5.0;
  std::vector<std::vector<double>> bcsq(2, std::vector<double>(2, cutoff * cutoff));
  analysis::StructureAnalyzer analyzer(cell, cutoff, bcsq, true);
  auto hist = HBondCalculator::calculate(cell, &analyzer);

  // Should produce a valid histogram without crashes
  for (const auto &[key, values] : hist.partials) {
    EXPECT_EQ(values.size(), hist.bins.size());
  }
}

/// Large system with many atoms but no valid H-bond geometry.
/// Tests performance and correctness under bulk non-bonding conditions.
TEST(HBondCalculatorTests, BulkMetalNoHBonds) {
  // 4-atom FCC copper — no H, O, N, F, S atoms at all
  double a = 3.6;
  correlation::core::Cell cell({a, a, a, 90.0, 90.0, 90.0});
  cell.addAtom("Cu", {0.0, 0.0, 0.0});
  cell.addAtom("Cu", {a / 2, a / 2, 0.0});
  cell.addAtom("Cu", {a / 2, 0.0, a / 2});
  cell.addAtom("Cu", {0.0, a / 2, a / 2});

  double cutoff = 3.0;
  std::vector<std::vector<double>> bcsq(1, std::vector<double>(1, cutoff * cutoff));
  analysis::StructureAnalyzer analyzer(cell, cutoff, bcsq, false);
  auto hist = HBondCalculator::calculate(cell, &analyzer);

  // Pure metal: no electronegative atoms, so no H-bond counting.
  // The histogram may still exist but should have no non-zero H-bond counts.
  for (size_t i = 0; i < hist.bins.size(); ++i) {
    if (hist.bins[i] > 0) {
      EXPECT_NEAR(hist.partials.at("Total")[i], 0.0, 1e-12)
          << "Pure metal should have no H-bonds.";
    }
  }
}

/// Deterministic: running twice gives the same result.
TEST(HBondCalculatorTests, DeterministicResults) {
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  cell.addAtom("O", {10.0, 10.0, 10.0});
  cell.addAtom("H", {10.0, 10.0, 10.96});
  cell.addAtom("O", {10.0, 10.0, 12.5});

  double cutoff = 5.0;
  std::vector<std::vector<double>> bcsq(2, std::vector<double>(2, cutoff * cutoff));
  analysis::StructureAnalyzer analyzer(cell, cutoff, bcsq, true);

  auto hist1 = HBondCalculator::calculate(cell, &analyzer);
  auto hist2 = HBondCalculator::calculate(cell, &analyzer);

  ASSERT_EQ(hist1.bins.size(), hist2.bins.size());
  for (size_t i = 0; i < hist1.bins.size(); ++i) {
    EXPECT_DOUBLE_EQ(hist1.bins[i], hist2.bins[i]);
  }
  for (const auto &[key, vals1] : hist1.partials) {
    ASSERT_TRUE(hist2.partials.count(key));
    const auto &vals2 = hist2.partials.at(key);
    ASSERT_EQ(vals1.size(), vals2.size());
    for (size_t i = 0; i < vals1.size(); ++i) {
      EXPECT_DOUBLE_EQ(vals1[i], vals2[i]);
    }
  }
}

} // namespace correlation::calculators
