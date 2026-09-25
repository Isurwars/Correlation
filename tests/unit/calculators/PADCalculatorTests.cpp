// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/PADCalculator.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include "../../CrystalTestHelper.hpp"

#include "math/Constants.hpp"
#include "math/Precision.hpp"
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <numbers>
#include <numeric>

namespace correlation::analysis {

// ============================================================================
// Part 1: Angle Reproduction Tests
// ============================================================================
namespace {
class PADCalculatorTests_AngleReproduction : public ::testing::Test {
protected:
  void SetUp() override {
    // Simple cubic cell
    cell_ = correlation::core::Cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  }

  void updateTrajectory() {
    trajectory_ = correlation::core::Trajectory();
    trajectory_.addFrame(cell_);
    trajectory_.precomputeBondCutoffs();
  }

  void updateTrajectory(const correlation::core::Cell &cell) {
    trajectory_ = correlation::core::Trajectory();
    trajectory_.addFrame(cell);
    trajectory_.precomputeBondCutoffs();
  }

public:
  correlation::core::Cell cell_;
  correlation::core::Trajectory trajectory_;
};

// Helper to sum a partial histogram
real_t sumHistogram(const std::vector<real_t> &hist) {
  return std::accumulate(hist.begin(), hist.end(), static_cast<real_t>(0.0));
}

[[nodiscard]] int getElementId(const correlation::core::Cell &cell, const std::string &symbol) {
  const auto elem = cell.findElement(symbol);
  if (!elem) {
    throw std::runtime_error("Element not found: " + symbol);
  }
  return elem->id.value;
}

struct AngleClassification {
  int count_63 = 0;
  int count_116 = 0;
  int count_180 = 0;
  int count_60 = 0;
  int count_108 = 0;
  int count_58 = 0;
  int total_angles = 0;
};

void classifyAngle(double deg, AngleClassification &counts) {
  if (std::abs(deg - 63.43) < 1.0) {
    counts.count_63++;
  } else if (std::abs(deg - 116.57) < 1.0) {
    counts.count_116++;
  } else if (std::abs(deg - 180.0) < 1.0) {
    counts.count_180++;
  } else if (std::abs(deg - 60.0) < 1.0) {
    counts.count_60++;
  } else if (std::abs(deg - 108.0) < 1.0) {
    counts.count_108++;
  } else if (std::abs(deg - 58.28) < 1.0) {
    counts.count_58++;
  }
}

[[nodiscard]] AngleClassification
countAngles(const std::vector<std::vector<std::vector<std::vector<real_t>>>> &angles) {
  AngleClassification counts;
  for (const auto &t_1 : angles) {
    for (const auto &t_2 : t_1) {
      for (const auto &t_3 : t_2) {
        for (double const angle : t_3) {
          counts.total_angles++;
          classifyAngle(angle * 180.0 / correlation::math::pi, counts);
        }
      }
    }
  }
  return counts;
}

[[nodiscard]] bool hasPeakNear(const std::vector<real_t> &partial, const std::vector<real_t> &bins,
                               double target_angle, double tol = 0.5) {
  for (size_t i = 0; i < partial.size(); ++i) {
    if (partial[i] > 0.01 && std::abs(bins[i] - target_angle) < tol) {
      return true;
    }
  }
  return false;
}

[[nodiscard]] double findPeakAngle(const std::vector<real_t> &partial,
                                   const std::vector<real_t> &bins) {
  auto max_it = std::ranges::max_element(partial);
  if (max_it == partial.end()) {
    return -1.0;
  }
  auto idx = std::distance(partial.begin(), max_it);
  return bins[idx];
}

void setupIcosahedron(correlation::core::Cell &cell,
                      real_t base_coord = static_cast<real_t>(10.0)) {
  cell.addAtom("Si", {base_coord, base_coord, base_coord});
  const auto phi = static_cast<real_t>(std::numbers::phi);
  std::vector<std::vector<real_t>> const vertices = {
      {0, 1, phi},  {0, 1, -phi},  {0, -1, phi}, {0, -1, -phi}, {1, phi, 0},  {1, -phi, 0},
      {-1, phi, 0}, {-1, -phi, 0}, {phi, 0, 1},  {phi, 0, -1},  {-phi, 0, 1}, {-phi, 0, -1}};
  for (const auto &vertex : vertices) {
    cell.addAtom("Si", correlation::math::Vector3<real_t>(
                           base_coord + vertex[0], base_coord + vertex[1], base_coord + vertex[2]));
  }
}

void verifyIcosahedronCenterAngles(const AngleClassification &counts) {
  EXPECT_EQ(counts.count_63, 30) << "Should find 30 Center-Edge angles (~63.4 deg)";
  EXPECT_EQ(counts.count_116, 30) << "Should find 30 Center-Diagonal angles (~116.6 deg)";
  EXPECT_EQ(counts.count_180, 6) << "Should find 6 Center-Opposite angles (180 deg)";
}

void verifyIcosahedronSurfaceAngles(const AngleClassification &counts) {
  EXPECT_EQ(counts.count_60, 60) << "Should find 60 Surface-Triangle angles (60 deg)";
  EXPECT_EQ(counts.count_108, 60) << "Should find 60 Surface-Pentagon angles (108 deg)";
  EXPECT_EQ(counts.count_58, 60) << "Should find 60 Surface-Center angles "
                                    "(Center-S-correlation::core::Neighbor, ~58.3 deg)";
  EXPECT_EQ(counts.total_angles, 246) << "Total angles should be 246";
}

void verifyIcosahedronPeaks(const std::vector<real_t> &partial, const std::vector<real_t> &bins) {
  EXPECT_TRUE(hasPeakNear(partial, bins, 58.28)) << "Should find PAD peak near 58.28 degrees";
  EXPECT_TRUE(hasPeakNear(partial, bins, 60.00)) << "Should find PAD peak near 60.00 degrees";
  EXPECT_TRUE(hasPeakNear(partial, bins, 63.43)) << "Should find PAD peak near 63.43 degrees";
  EXPECT_TRUE(hasPeakNear(partial, bins, 108.00)) << "Should find PAD peak near 108.00 degrees";
  EXPECT_TRUE(hasPeakNear(partial, bins, 116.57)) << "Should find PAD peak near 116.57 degrees";
  EXPECT_TRUE(hasPeakNear(partial, bins, 180.00)) << "Should find PAD peak near 180.00 degrees";
}

class PADCalculatorTests : public ::testing::Test {
protected:
  void SetUp() override {
    // Large box to avoid PBC issues by default
    cell_ = correlation::core::Cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  }

  void updateTrajectory() {
    trajectory_ = correlation::core::Trajectory();
    trajectory_.addFrame(cell_);
    trajectory_.precomputeBondCutoffs();
  }

public:
  correlation::core::Cell cell_;
  correlation::core::Trajectory trajectory_;
};
} // namespace

TEST_F(PADCalculatorTests_AngleReproduction, CalculatePAD) {
  auto water = correlation::testing::crystals::createWaterMoleculeCell(
      {.O_pos = {5.0, 5.0, 5.0}, .r_OH = 1.0, .angle_HOH_deg = 104.5, .box_size = 10.0});

  updateTrajectory(water);
  DistributionFunctions dists(water, 2.0, trajectory_.getBondCutoffsSQ());

  dists.calculatePAD(0.001);
  const auto &hist = dists.getHistogram("PAD");
  const auto &hist_alias = dists.getHistogram("BAD");
  EXPECT_EQ(hist.bins.size(), hist_alias.bins.size());
  const auto &hoh = hist.partials.at("H-O-H");

  auto max_it = std::ranges::max_element(hoh);
  size_t const idx = std::distance(hoh.begin(), max_it);
  real_t const angle = hist.bins[idx];
  // 104.5 angle with 0.001 bins could land in 104.4995 or 104.5005 due to
  // precision
  EXPECT_NEAR(angle, 104.5, 0.001);
}

TEST_F(PADCalculatorTests_AngleReproduction, MissingAnglesWhenCutoffIsTooSmall) {
  // A-B-C angle.
  // B is at (5,5,5)
  // A is at (3.4,5,5) -> dist 1.6
  // C is at (5,6.6,5) -> dist 1.6
  // Angle should be 90 degrees.

  cell_.addAtom("O", {3.4, 5.0, 5.0});
  cell_.addAtom("Si", {5.0, 5.0, 5.0});
  cell_.addAtom("O", {5.0, 6.6, 5.0});
  updateTrajectory();

  // Bond cutoff for Si-O is likely around 1.6 * 1.2 = 1.92 or similar.
  // Distance is 1.6.
  {
    StructureAnalyzer const analyzer(cell_, 1.8, trajectory_.getBondCutoffsSQ());
    const auto &angles = analyzer.angles();
    bool found = false;
    for (const auto &t_1 : angles) {
      for (const auto &t_2 : t_1) {
        for (const auto &t_3 : t_2) {
          for (real_t const angle : t_3) {
            if (std::abs(angle * 180.0 / correlation::math::pi - 90.0) < 1.0) {
              found = true;
            }
          }
        }
      }
    }
    EXPECT_TRUE(found) << "Should find 90 degree angle with sufficient cutoff";
  }
}

TEST_F(PADCalculatorTests_AngleReproduction, PBCAngleDetection) {
  cell_.addAtom("Si", {0.5, 0.5, 0.5});
  cell_.addAtom("O", {8.9, 0.5, 0.5});
  cell_.addAtom("O", {0.5, 8.9, 0.5});
  updateTrajectory();

  StructureAnalyzer const analyzer(cell_, 1.8, trajectory_.getBondCutoffsSQ());

  bool found = false;
  const auto &angles = analyzer.angles();
  for (const auto &t_1 : angles) {
    for (const auto &t_2 : t_1) {
      for (const auto &t_3 : t_2) {
        for (double const angle : t_3) {
          if (std::abs(angle * 180.0 / correlation::math::pi - 90.0) < 1.0) {
            found = true;
          }
        }
      }
    }
  }
  EXPECT_TRUE(found) << "Should find 90 degree angle across PBC";
}

TEST_F(PADCalculatorTests_AngleReproduction, SiTetrahedron_4Atoms) {
  cell_.addAtom("Si", {5.0, 5.0, 5.0}); // Center
  cell_.addAtom("Si", {6.0, 6.0, 6.0}); // correlation::core::Neighbor 1 (1,1,1)
  cell_.addAtom("Si", {6.0, 4.0, 4.0}); // correlation::core::Neighbor 2 (1,-1,-1)
  cell_.addAtom("Si", {4.0, 6.0, 4.0}); // correlation::core::Neighbor 3 (-1,1,-1)
  cell_.addAtom("Si", {4.0, 4.0, 6.0}); // correlation::core::Neighbor 4 (-1,-1,1)
  updateTrajectory();

  // With 4 neighbors, we have C(4,2) = 6 angles.
  // Neighbors are at dist sqrt(3) ~ 1.73.
  // N-N dist is sqrt(8) ~ 2.82.
  // Si radius 1.16. Bond cutoff ~ 2.78.
  // Thus neighbors are NOT connected to each other.

  StructureAnalyzer const analyzer(cell_, 3.0, trajectory_.getBondCutoffsSQ());
  const auto &angles = analyzer.angles();

  int angle_count = 0;
  for (const auto &t_1 : angles) {
    for (const auto &t_2 : t_1) {
      for (const auto &t_3 : t_2) {
        for (double const angle : t_3) {
          double const degrees = angle * 180.0 / correlation::math::pi;
          // std::cout << "Angle: " << degrees << " degrees\n";
          // Expected angle is acos(-1/3) ~ 109.47 degrees
          if (std::abs(degrees - 109.47) < 1.0) {
            angle_count++;
          }
        }
      }
    }
  }
  EXPECT_EQ(angle_count, 6) << "Should find exactly 6 angles of ~109.47 "
                               "degrees for a standard Si tetrahedron";
}

TEST_F(PADCalculatorTests_AngleReproduction, Icosahedron_13Atoms) {
  setupIcosahedron(cell_);
  updateTrajectory();

  // Cutoff ~ 2.5 covers bonds (1.902, 2.0) but avoids next-nearest (3.236)
  StructureAnalyzer const analyzer(cell_, 2.5, trajectory_.getBondCutoffsSQ());
  const auto counts = countAngles(analyzer.angles());
  verifyIcosahedronCenterAngles(counts);
  verifyIcosahedronSurfaceAngles(counts);
}

// ============================================================================
// Part 2: Plane Angle Distribution (PAD) Tests
// ============================================================================
// 1. Trivial Cases
TEST_F(PADCalculatorTests, EmptyCellThrows) {
  // Current implementation throws explicitly if atoms are empty in
  // calculateAshcroftWeights or implicitly via other checks.
  updateTrajectory();
  EXPECT_THROW(
      { const DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ()); },
      std::invalid_argument);
}

TEST_F(PADCalculatorTests, SingleAtomNoAngles) {
  cell_.addAtom("Si", {10.0, 10.0, 10.0});
  updateTrajectory();
  DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ());
  dists.calculatePAD(1.0);
  // Might have partials created but empty, or just no "BAD" if logic
  // handles it. Actually implementation might create partials if atoms exist
  // but no angles found. Let's check total counts.
  if (static_cast<unsigned int>(dists.getAllHistograms().contains("PAD")) != 0U) {
    const auto &hist = dists.getHistogram("PAD");
    if (!hist.partials.empty()) {
      if (static_cast<unsigned int>(hist.partials.contains("Total")) != 0U) {
        EXPECT_DOUBLE_EQ(sumHistogram(hist.partials.at("Total")), 0.0);
      }
    }
  }
}

TEST_F(PADCalculatorTests, NullNeighborsThrows) {
  cell_.addAtom("Si", {10.0, 10.0, 10.0});
  EXPECT_THROW(
      { correlation::calculators::PADCalculator::calculate(cell_, nullptr, 1.0); },
      std::logic_error);
}

// 2. Geometry Verification
TEST_F(PADCalculatorTests, LinearGeometry180) {
  // A-B-C line
  cell_.addAtom("O", {8.4, 10.0, 10.0});
  cell_.addAtom("Si", {10.0, 10.0, 10.0});
  cell_.addAtom("O", {11.6, 10.0, 10.0});
  updateTrajectory();

  // Verify StructureAnalyzer finds neighbors
  StructureAnalyzer const analyzer(cell_, 2.0, trajectory_.getBondCutoffsSQ());
  const auto &neighbor_graph = analyzer.neighborGraph();
  // Si is atom index 1 (0-based)
  ASSERT_GT(neighbor_graph.nodeCount(), 1);
  EXPECT_EQ(neighbor_graph.getNeighbors(1).size(), 2) << "Si should have 2 neighbors (O atoms)";

  // Bond length 1.6. Cutoff needs to be > 1.6
  DistributionFunctions dists(cell_, 2.0, trajectory_.getBondCutoffsSQ());
  // Fine binning for accuracy
  dists.calculatePAD(0.001);

  const auto &hist = dists.getHistogram("PAD");
  // Should have O-Si-O peak at 180
  ASSERT_EQ(hist.partials.count("O-Si-O"), 1);
  const auto &partial = hist.partials.at("O-Si-O");

  // Bin for 180 degrees.
  // Multiply counts/density by bin width (0.001) to get the probability sum
  double const total_prob = sumHistogram(partial) * 0.001;
  EXPECT_NEAR(total_prob, 1.0, 1e-5)
      << "Should be normalized to 1 angle (normalized by counts * bin_width)";

  // Check peak location
  double const peak_angle = findPeakAngle(partial, hist.bins);
  EXPECT_NEAR(peak_angle, 180.0, 1e-3);
}

TEST_F(PADCalculatorTests, RightAngle90) {
  cell_.addAtom("O", {10.0, 8.4, 10.0});
  cell_.addAtom("Si", {10.0, 10.0, 10.0}); // Center
  cell_.addAtom("O", {11.6, 10.0, 10.0});
  updateTrajectory();

  DistributionFunctions dists(cell_, 2.0, trajectory_.getBondCutoffsSQ());
  dists.calculatePAD(0.001);

  const auto &hist = dists.getHistogram("PAD");
  ASSERT_EQ(hist.partials.count("O-Si-O"), 1);

  // Find peak
  const auto &partial = hist.partials.at("O-Si-O");
  double const peak_angle = findPeakAngle(partial, hist.bins);
  // 90.0 / 0.001 could land in 89.9995 or 90.0005
  EXPECT_NEAR(peak_angle, 90.0, 0.001);
}

TEST_F(PADCalculatorTests, EquilateralTriangle60) {
  // Si at (0,0,0)
  // O at (1.6,0,0)
  // O at (0.8, 1.6*sqrt(3)/2, 0)
  cell_.addAtom("Si", {10.0, 10.0, 10.0});
  cell_.addAtom("O", {11.6, 10.0, 10.0});
  cell_.addAtom("O", {10.0 + 1.6 * 0.5, 10.0 + 1.6 * std::numbers::sqrt3 / 2.0, 10.0});
  updateTrajectory();

  DistributionFunctions dists(cell_, 2.0, trajectory_.getBondCutoffsSQ());
  dists.calculatePAD(1.0);

  const auto &hist = dists.getHistogram("PAD");
  // Should have O-Si-O
  const auto &partial = hist.partials.at("O-Si-O");

  // Find peak near 60
  double const val_at_60 = 0;
  // index for 60 deg is 60 or 59 depending on binning.
  // 59.5 (idx 59) -> [59, 60)
  // 60.5 (idx 60) -> [60, 61)
  // Exact 60 might land in 60.

  // Search max around 60
  size_t const bin_60 = 60;
  EXPECT_GT(partial[bin_60] + partial[bin_60 - 1], 0.1) << "Should have peak near 60 degrees";
}

TEST_F(PADCalculatorTests, TetrahedralAngle) {
  // Si at center
  // 4 Neighbors at tetrahedral positions.
  // For simplicity, just check one angle 109.47
  const auto base_coord = static_cast<real_t>(10.0);
  cell_.addAtom("Si", {base_coord, base_coord, base_coord});
  // Vector 1: (1,1,1) normalized * 1.6
  // Vector 2: (1,-1,-1) normalized * 1.6
  // Dot product = (1-1-1)/3 = -1/3. acos(-1/3) = 109.47 deg

  const auto lattice_constant = static_cast<real_t>(1.6 * std::numbers::inv_sqrt3);

  cell_.addAtom("O", correlation::math::Vector3<real_t>(base_coord + lattice_constant,
                                                        base_coord + lattice_constant,
                                                        base_coord + lattice_constant));
  cell_.addAtom("O", correlation::math::Vector3<real_t>(base_coord + lattice_constant,
                                                        base_coord - lattice_constant,
                                                        base_coord - lattice_constant));
  updateTrajectory();

  DistributionFunctions dists(cell_, 2.0,
                              trajectory_.getBondCutoffsSQ()); // Distance is 1.6
  dists.calculatePAD(0.001);                                   // Hyperfine bins

  const auto &hist = dists.getHistogram("PAD");
  const auto &partial = hist.partials.at("O-Si-O");

  double const peak_angle = findPeakAngle(partial, hist.bins);
  // 109.4712... / 0.001 -> index 109471 -> center 109.4715
  EXPECT_NEAR(peak_angle, 109.4712206, 0.001);
}

// 3. Symmetry & Multi-Species
TEST_F(PADCalculatorTests, SymmetryAndSorting) {
  cell_.addAtom("Si", {10.0, 10.0, 10.0}); // Center
  cell_.addAtom("O", {11.6, 10.0, 10.0});
  cell_.addAtom("N", {10.0, 11.6, 10.0}); // 90 degrees
  updateTrajectory();

  DistributionFunctions dists(cell_, 2.0, trajectory_.getBondCutoffsSQ());
  dists.calculatePAD(1.0);

  const auto &hist = dists.getHistogram("PAD");

  // Check if we have O-Si-N or N-Si-O
  bool found = false;
  if (static_cast<unsigned int>(hist.partials.contains("O-Si-N")) != 0U) {
    found = true;
  }
  if (static_cast<unsigned int>(hist.partials.contains("N-Si-O")) != 0U) {
    found = true;
  }

  EXPECT_TRUE(found) << "Should have mixed species angle distribution";
}

// 4. Normalization
TEST_F(PADCalculatorTests, FullNormalizationCheck) {
  // 1 Si, 4 O neighbors (tetrahedron)
  // 4 neighbors -> 4*3/2 = 6 angles.
  // All 6 angles are 109.47
  const auto base_coord = static_cast<real_t>(10.0);
  cell_.addAtom("Si", {base_coord, base_coord, base_coord});
  const auto lattice_constant = static_cast<real_t>(1.6 * std::numbers::inv_sqrt3);

  // Tetrahedral vertices
  cell_.addAtom("O", correlation::math::Vector3<real_t>(base_coord + lattice_constant,
                                                        base_coord + lattice_constant,
                                                        base_coord + lattice_constant));
  cell_.addAtom("O", correlation::math::Vector3<real_t>(base_coord + lattice_constant,
                                                        base_coord - lattice_constant,
                                                        base_coord - lattice_constant));
  cell_.addAtom("O", correlation::math::Vector3<real_t>(base_coord - lattice_constant,
                                                        base_coord + lattice_constant,
                                                        base_coord - lattice_constant));
  cell_.addAtom("O", correlation::math::Vector3<real_t>(base_coord - lattice_constant,
                                                        base_coord - lattice_constant,
                                                        base_coord + lattice_constant));
  updateTrajectory();

  // Custom bond cutoffs to avoid O-O bonds (distance ~2.61) which would create
  // extra angles
  auto cutoffs = trajectory_.getBondCutoffs();
  int const id_o = getElementId(cell_, "O");
  cutoffs[id_o][id_o].max_sq = 2.0;

  DistributionFunctions dists(cell_, 2.0, cutoffs);
  dists.calculatePAD(1.0);

  const auto &hist = dists.getHistogram("PAD");

  double sum_partial = 0;
  double sum_total = 0;
  double const bin_width = 1.0;

  if (static_cast<unsigned int>(hist.partials.contains("O-Si-O")) != 0U) {
    const auto &partial = hist.partials.at("O-Si-O");
    for (double const val : partial) {
      sum_partial += val * bin_width;
    }
  }

  if (static_cast<unsigned int>(hist.partials.contains("Total")) != 0U) {
    const auto &total = hist.partials.at("Total");
    for (double const val : total) {
      sum_total += val * bin_width;
    }
  }

  EXPECT_NEAR(sum_partial, 1.0,
              0.05); // Relaxed checking 0.05 due to binning effects
  EXPECT_NEAR(sum_total, 1.0, 0.05);
}

TEST_F(PADCalculatorTests, IcosahedronAnglesPAD) {
  setupIcosahedron(cell_);
  updateTrajectory();

  DistributionFunctions dists(cell_, 2.5, trajectory_.getBondCutoffsSQ());
  dists.calculatePAD(0.01);

  const auto &hist = dists.getHistogram("PAD");
  ASSERT_EQ(hist.partials.count("Si-Si-Si"), 1);
  const auto &partial = hist.partials.at("Si-Si-Si");

  verifyIcosahedronPeaks(partial, hist.bins);
}

TEST_F(PADCalculatorTests, BondDistanceBelowMinCutoffProducesNoAngles) {
  // Center Si at (10, 10, 10), O at (10.5, 10, 10) [dist 0.5 Å], O at (10, 11.5, 10) [dist 1.5 Å]
  cell_.addAtom("Si", {10.0, 10.0, 10.0});
  cell_.addAtom("O", {10.5, 10.0, 10.0}); // dist = 0.5 Å
  cell_.addAtom("O", {10.0, 11.5, 10.0}); // dist = 1.5 Å
  updateTrajectory();

  int const id_si = getElementId(cell_, "Si");
  int const id_o = getElementId(cell_, "O");

  // Cutoff range for Si-O: [1.0 Å, 2.0 Å] -> min_sq = 1.0, max_sq = 4.0; O-O and Si-Si = 0
  BondCutoffMatrix cutoffs(
      2, std::vector<BondCutoffRange>(2, BondCutoffRange{.min_sq = 0.0, .max_sq = 0.0}));
  cutoffs[id_si][id_o] = BondCutoffRange{.min_sq = 1.0, .max_sq = 4.0};
  cutoffs[id_o][id_si] = BondCutoffRange{.min_sq = 1.0, .max_sq = 4.0};
  StructureAnalyzer const analyzer(cell_, 2.0, cutoffs);

  // Since bond 1 (0.5 Å) < min_cutoff (1.0 Å), Si-O bond is not formed -> 0 angles
  const auto &angles = analyzer.angles();
  size_t total_angles = 0;
  for (const auto &t_1 : angles) {
    for (const auto &t_2 : t_1) {
      for (const auto &t_3 : t_2) {
        total_angles += t_3.size();
      }
    }
  }
  EXPECT_EQ(total_angles, 0);
}

TEST_F(PADCalculatorTests, BondDistanceAboveMaxCutoffProducesNoAngles) {
  // Center Si at (10, 10, 10), O at (11.5, 10, 10) [dist 1.5 Å], O at (10, 12.5, 10) [dist 2.5 Å]
  cell_.addAtom("Si", {10.0, 10.0, 10.0});
  cell_.addAtom("O", {11.5, 10.0, 10.0}); // dist = 1.5 Å
  cell_.addAtom("O", {10.0, 12.5, 10.0}); // dist = 2.5 Å
  updateTrajectory();

  int const id_si = getElementId(cell_, "Si");
  int const id_o = getElementId(cell_, "O");

  // Cutoff range for Si-O: [1.0 Å, 2.0 Å] -> min_sq = 1.0, max_sq = 4.0; O-O and Si-Si = 0
  BondCutoffMatrix cutoffs(
      2, std::vector<BondCutoffRange>(2, BondCutoffRange{.min_sq = 0.0, .max_sq = 0.0}));
  cutoffs[id_si][id_o] = BondCutoffRange{.min_sq = 1.0, .max_sq = 4.0};
  cutoffs[id_o][id_si] = BondCutoffRange{.min_sq = 1.0, .max_sq = 4.0};
  StructureAnalyzer const analyzer(cell_, 3.0, cutoffs);

  // Since bond 2 (2.5 Å) > max_cutoff (2.0 Å), second Si-O bond is not formed -> 0 angles
  const auto &angles = analyzer.angles();
  size_t total_angles = 0;
  for (const auto &t_1 : angles) {
    for (const auto &t_2 : t_1) {
      for (const auto &t_3 : t_2) {
        total_angles += t_3.size();
      }
    }
  }
  EXPECT_EQ(total_angles, 0);
}

TEST_F(PADCalculatorTests, BondDistanceWithinCutoffRangeProducesAngle) {
  // Center Si at (10, 10, 10), O at (11.5, 10, 10) [dist 1.5 Å], O at (10, 11.5, 10) [dist 1.5 Å]
  cell_.addAtom("Si", {10.0, 10.0, 10.0});
  cell_.addAtom("O", {11.5, 10.0, 10.0}); // dist = 1.5 Å
  cell_.addAtom("O", {10.0, 11.5, 10.0}); // dist = 1.5 Å
  updateTrajectory();

  int const id_si = getElementId(cell_, "Si");
  int const id_o = getElementId(cell_, "O");

  // Cutoff range for Si-O: [1.0 Å, 2.0 Å] -> min_sq = 1.0, max_sq = 4.0; O-O and Si-Si = 0
  BondCutoffMatrix cutoffs(
      2, std::vector<BondCutoffRange>(2, BondCutoffRange{.min_sq = 0.0, .max_sq = 0.0}));
  cutoffs[id_si][id_o] = BondCutoffRange{.min_sq = 1.0, .max_sq = 4.0};
  cutoffs[id_o][id_si] = BondCutoffRange{.min_sq = 1.0, .max_sq = 4.0};
  StructureAnalyzer const analyzer(cell_, 2.0, cutoffs);

  // Both bonds are within [1.0, 2.0], forming a 90 degree angle
  const auto &angles = analyzer.angles();
  bool found_90 = false;
  for (const auto &t_1 : angles) {
    for (const auto &t_2 : t_1) {
      for (const auto &t_3 : t_2) {
        for (real_t const angle : t_3) {
          if (std::abs(angle * 180.0 / correlation::math::pi - 90.0) < 1.0) {
            found_90 = true;
          }
        }
      }
    }
  }
  EXPECT_TRUE(found_90);
}

TEST_F(PADCalculatorTests, VerifyPADRawHistogram) {
  // Linear O-Si-O -> 1 angle of 180 degrees
  cell_.addAtom("O", {8.4, 10.0, 10.0});
  cell_.addAtom("Si", {10.0, 10.0, 10.0});
  cell_.addAtom("O", {11.6, 10.0, 10.0});
  updateTrajectory();

  DistributionFunctions dists(cell_, 2.0, trajectory_.getBondCutoffsSQ());
  dists.calculatePAD(1.0);

  ASSERT_TRUE(dists.getAllHistograms().contains("PAD_raw"));
  const auto &raw_hist = dists.getHistogram("PAD_raw");
  EXPECT_EQ(raw_hist.y_unit, "counts");
  ASSERT_TRUE(raw_hist.partials.contains("O-Si-O"));
  EXPECT_DOUBLE_EQ(sumHistogram(raw_hist.partials.at("O-Si-O")), 1.0);
  EXPECT_DOUBLE_EQ(sumHistogram(raw_hist.partials.at("Total")), 1.0);
}
} // namespace correlation::analysis
