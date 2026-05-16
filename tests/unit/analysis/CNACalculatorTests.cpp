// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/StructureAnalyzer.hpp"
#include "calculators/CNACalculator.hpp"
#include "core/Cell.hpp"

#include <gtest/gtest.h>

namespace correlation::calculators {

// ============================================================================
// Null / Empty Input
// ============================================================================

/// Passing a null neighbor pointer should return an empty histogram.
TEST(CNACalculatorTests, NullNeighborsReturnsEmptyHistogram) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});

  auto hist = CNACalculator::calculate(cell, nullptr);

  EXPECT_TRUE(hist.bins.empty());
  EXPECT_TRUE(hist.partials.empty());
}

/// A single isolated atom with no neighbors produces no CNA pairs.
TEST(CNACalculatorTests, SingleIsolatedAtomReturnsEmptyHistogram) {
  correlation::core::Cell cell({100.0, 100.0, 100.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {50.0, 50.0, 50.0});

  // cutoff is tiny: no neighbors found
  analysis::StructureAnalyzer analyzer(cell, 0.5, {{0.25}}, true);
  auto hist = CNACalculator::calculate(cell, &analyzer);

  EXPECT_TRUE(hist.bins.empty());
  EXPECT_TRUE(hist.partials.empty());
}

// ============================================================================
// Histogram Consistency
// ============================================================================

/// All partial vectors must have the same size as the bins vector.
/// This was the crash bug fixed earlier — ensure it never regresses.
TEST(CNACalculatorTests, HistogramDimensionsAreConsistent) {
  // FCC conventional cell (4 atoms in a cube)
  correlation::core::Cell cell({4.0, 4.0, 4.0, 90.0, 90.0, 90.0});
  cell.addAtom("Al", {0.0, 0.0, 0.0});
  cell.addAtom("Al", {2.0, 2.0, 0.0});
  cell.addAtom("Al", {2.0, 0.0, 2.0});
  cell.addAtom("Al", {0.0, 2.0, 2.0});

  double cutoff = 3.0;
  analysis::StructureAnalyzer analyzer(cell, cutoff,
                                       {{cutoff * cutoff}}, false);
  auto hist = CNACalculator::calculate(cell, &analyzer);

  ASSERT_FALSE(hist.bins.empty());

  for (const auto &[key, values] : hist.partials) {
    EXPECT_EQ(values.size(), hist.bins.size())
        << "Partial '" << key << "' has size " << values.size()
        << " but bins has size " << hist.bins.size();
  }
}

/// The "Total" partial must sum to exactly 1.0 (within fp tolerance).
TEST(CNACalculatorTests, TotalPartialSumsToOne) {
  correlation::core::Cell cell({4.0, 4.0, 4.0, 90.0, 90.0, 90.0});
  cell.addAtom("Cu", {0.0, 0.0, 0.0});
  cell.addAtom("Cu", {2.0, 2.0, 0.0});
  cell.addAtom("Cu", {2.0, 0.0, 2.0});
  cell.addAtom("Cu", {0.0, 2.0, 2.0});

  double cutoff = 3.0;
  analysis::StructureAnalyzer analyzer(cell, cutoff,
                                       {{cutoff * cutoff}}, false);
  auto hist = CNACalculator::calculate(cell, &analyzer);

  ASSERT_TRUE(hist.partials.count("Total"));
  double sum = 0.0;
  for (double v : hist.partials.at("Total")) {
    sum += v;
  }
  EXPECT_NEAR(sum, 1.0, 1e-12);
}

/// Histogram metadata is correctly populated.
TEST(CNACalculatorTests, HistogramMetadataIsPopulated) {
  correlation::core::Cell cell({4.0, 4.0, 4.0, 90.0, 90.0, 90.0});
  cell.addAtom("Cu", {0.0, 0.0, 0.0});
  cell.addAtom("Cu", {2.0, 2.0, 0.0});

  double cutoff = 3.0;
  analysis::StructureAnalyzer analyzer(cell, cutoff,
                                       {{cutoff * cutoff}}, false);
  auto hist = CNACalculator::calculate(cell, &analyzer);

  EXPECT_EQ(hist.title, "Common Neighbor Analysis");
  EXPECT_EQ(hist.x_label, "CNA Index");
  EXPECT_EQ(hist.y_label, "Frequency");
  EXPECT_EQ(hist.file_suffix, "_CNA");
  EXPECT_FALSE(hist.description.empty());
}

// ============================================================================
// Crystal Structure Identification
// ============================================================================

/// FCC 2×2×2 supercell (32 atoms) — standard CNA should produce
/// non-empty results. This is a structural integration test verifying that
/// the calculator runs correctly on a realistic periodic FCC system.
TEST(CNACalculatorTests, FCC_Supercell_ProducesNonEmptyResult) {
  double a = 4.0;
  double L = 2.0 * a; // 2×2×2 supercell
  correlation::core::Cell cell({L, L, L, 90.0, 90.0, 90.0});

  // Generate 2×2×2 FCC supercell (4 basis atoms × 8 unit cells = 32 atoms)
  for (int ix = 0; ix < 2; ++ix)
    for (int iy = 0; iy < 2; ++iy)
      for (int iz = 0; iz < 2; ++iz) {
        double ox = ix * a, oy = iy * a, oz = iz * a;
        cell.addAtom("Al", {ox, oy, oz});
        cell.addAtom("Al", {ox + a / 2, oy + a / 2, oz});
        cell.addAtom("Al", {ox + a / 2, oy, oz + a / 2});
        cell.addAtom("Al", {ox, oy + a / 2, oz + a / 2});
      }

  // 1st NN distance = a/sqrt(2) ≈ 2.83; cutoff between 1st and 2nd NN
  double cutoff = 3.2;
  analysis::StructureAnalyzer analyzer(cell, cutoff,
                                       {{cutoff * cutoff}}, false);
  auto hist = CNACalculator::calculate(cell, &analyzer);

  // Should produce non-empty results
  ASSERT_FALSE(hist.bins.empty()) << "FCC supercell should produce CNA output.";

  // 1421 should be the dominant index for FCC
  EXPECT_TRUE(hist.partials.count("1421"))
      << "FCC should produce CNA index 1421.";

  // The Total partial should sum to 1.0
  double sum = 0;
  for (double v : hist.partials.at("Total"))
    sum += v;
  EXPECT_NEAR(sum, 1.0, 1e-12);
}

/// BCC 3×3×3 supercell (54 atoms) including both 1st and 2nd NN shells.
/// BCC CNA requires 2nd shell neighbors as common neighbors for 1st-shell pairs.
/// With a cutoff that captures both shells (14 neighbors per atom), the
/// resulting CNA indices depend on the full neighbor topology.
TEST(CNACalculatorTests, BCC_Supercell_ProducesOutput) {
  double a = 3.0;
  double L = 3.0 * a;
  correlation::core::Cell cell({L, L, L, 90.0, 90.0, 90.0});

  for (int ix = 0; ix < 3; ++ix)
    for (int iy = 0; iy < 3; ++iy)
      for (int iz = 0; iz < 3; ++iz) {
        double ox = ix * a, oy = iy * a, oz = iz * a;
        cell.addAtom("Fe", {ox, oy, oz});
        cell.addAtom("Fe", {ox + a / 2, oy + a / 2, oz + a / 2});
      }

  // Include both 1st NN (a*sqrt(3)/2 ≈ 2.598) and 2nd NN (a = 3.0)
  double cutoff = 3.1;
  analysis::StructureAnalyzer analyzer(cell, cutoff,
                                       {{cutoff * cutoff}}, false);
  auto hist = CNACalculator::calculate(cell, &analyzer);

  ASSERT_FALSE(hist.bins.empty()) << "BCC supercell should produce CNA output.";

  // Validate histogram consistency
  for (const auto &[key, vals] : hist.partials) {
    EXPECT_EQ(vals.size(), hist.bins.size())
        << "Partial '" << key << "' size mismatch.";
  }

  // Total should sum to 1.0
  double sum = 0;
  for (double v : hist.partials.at("Total"))
    sum += v;
  EXPECT_NEAR(sum, 1.0, 1e-12);
}

// ============================================================================
// Two Atoms – Bonded But No Common Neighbors
// ============================================================================

/// Two bonded atoms with no common neighbors should produce no CNA output.
TEST(CNACalculatorTests, TwoAtomsNoCommonNeighbors) {
  correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {10.0, 10.0, 10.0});
  cell.addAtom("Ar", {12.0, 10.0, 10.0}); // 2 Å apart

  double cutoff = 3.0;
  analysis::StructureAnalyzer analyzer(cell, cutoff,
                                       {{cutoff * cutoff}}, true);
  auto hist = CNACalculator::calculate(cell, &analyzer);

  // Two atoms that are neighbors but share no common neighbors ⟹ no CNA pairs
  EXPECT_TRUE(hist.bins.empty());
}

// ============================================================================
// Deterministic Ordering
// ============================================================================

/// Running CNA twice on the same structure should produce identical results.
TEST(CNACalculatorTests, DeterministicResults) {
  double a = 4.0;
  double L = 2.0 * a;
  correlation::core::Cell cell({L, L, L, 90.0, 90.0, 90.0});
  for (int ix = 0; ix < 2; ++ix)
    for (int iy = 0; iy < 2; ++iy)
      for (int iz = 0; iz < 2; ++iz) {
        double ox = ix * a, oy = iy * a, oz = iz * a;
        cell.addAtom("Cu", {ox, oy, oz});
        cell.addAtom("Cu", {ox + a / 2, oy + a / 2, oz});
        cell.addAtom("Cu", {ox + a / 2, oy, oz + a / 2});
        cell.addAtom("Cu", {ox, oy + a / 2, oz + a / 2});
      }

  double cutoff = 3.2;
  analysis::StructureAnalyzer analyzer(cell, cutoff,
                                       {{cutoff * cutoff}}, false);

  auto hist1 = CNACalculator::calculate(cell, &analyzer);
  auto hist2 = CNACalculator::calculate(cell, &analyzer);

  ASSERT_EQ(hist1.bins.size(), hist2.bins.size());
  ASSERT_EQ(hist1.partials.size(), hist2.partials.size());

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

/// Tiny periodic cell: 4-atom FCC unit cell where every atom neighbors every
/// other through multiple images. Verifies CNA doesn't crash on this edge case
/// and still produces consistent histograms.
TEST(CNACalculatorTests, TinyPeriodicCell_DoesNotCrash) {
  double a = 4.0;
  correlation::core::Cell cell({a, a, a, 90.0, 90.0, 90.0});
  cell.addAtom("Cu", {0.0, 0.0, 0.0});
  cell.addAtom("Cu", {a / 2, a / 2, 0.0});
  cell.addAtom("Cu", {a / 2, 0.0, a / 2});
  cell.addAtom("Cu", {0.0, a / 2, a / 2});

  double cutoff = 3.2;
  analysis::StructureAnalyzer analyzer(cell, cutoff,
                                       {{cutoff * cutoff}}, false);

  ASSERT_NO_THROW({
    auto hist = CNACalculator::calculate(cell, &analyzer);
    // Whatever indices are produced, they must be dimensionally consistent
    for (const auto &[key, vals] : hist.partials) {
      EXPECT_EQ(vals.size(), hist.bins.size())
          << "Partial '" << key << "' has inconsistent size.";
    }
  });
}

} // namespace correlation::calculators
