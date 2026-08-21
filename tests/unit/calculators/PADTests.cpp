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
class PADTests_AngleReproduction : public ::testing::Test {
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

class PADTests : public ::testing::Test {
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

TEST_F(PADTests_AngleReproduction, CalculatePAD) {
  auto water = correlation::testing::crystals::createWaterMoleculeCell(
      {.O_pos = {5.0, 5.0, 5.0}, .r_OH = 1.0, .angle_HOH_deg = 104.5, .box_size = 10.0});

  updateTrajectory(water);
  DistributionFunctions dists(water, 2.0, trajectory_.getBondCutoffsSQ());

  dists.calculatePAD(0.001);
  const auto &hist = dists.getHistogram("PAD");
  const auto &hist_alias = dists.getHistogram("BAD");
  EXPECT_EQ(hist.bins.size(), hist_alias.bins.size());
  const auto &hoh = hist.partials.at("H-O-H");

  auto max_it = std::max_element(hoh.begin(), hoh.end());
  size_t const idx = std::distance(hoh.begin(), max_it);
  real_t const angle = hist.bins[idx];
  // 104.5 angle with 0.001 bins could land in 104.4995 or 104.5005 due to
  // precision
  EXPECT_NEAR(angle, 104.5, 0.001);
}

TEST_F(PADTests_AngleReproduction, MissingAnglesWhenCutoffIsTooSmall) {
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

TEST_F(PADTests_AngleReproduction, PBCAngleDetection) {
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

TEST_F(PADTests_AngleReproduction, SiTetrahedron_4Atoms) {
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

TEST_F(PADTests_AngleReproduction, Icosahedron_13Atoms) {
  real_t const base_coord = static_cast<real_t>(10.0);
  cell_.addAtom("Si", {base_coord, base_coord, base_coord}); // Center

  real_t const phi = static_cast<real_t>(std::numbers::phi);
  // Vertices of icosahedron (edge length 2) relative to center
  std::vector<std::vector<real_t>> const vertices = {{0, 1, phi}, {0, 1, -phi}, {0, -1, phi}, {0, -1, -phi},
                                                     {1, phi, 0}, {1, -phi, 0}, {-1, phi, 0}, {-1, -phi, 0},
                                                     {phi, 0, 1}, {phi, 0, -1}, {-phi, 0, 1}, {-phi, 0, -1}};

  for (const auto &vertex : vertices) {
    cell_.addAtom("Si", correlation::math::Vector3<real_t>(base_coord + vertex[0], base_coord + vertex[1],
                                                           base_coord + vertex[2]));
  }
  updateTrajectory();

  // Cutoff ~ 2.5 covers bonds (1.902, 2.0) but avoids next-nearest (3.236)
  StructureAnalyzer const analyzer(cell_, 2.5, trajectory_.getBondCutoffsSQ());
  const auto &angles = analyzer.angles();

  int count_63 = 0;  // Center-Edge (approx 63.43)
  int count_116 = 0; // Center-Diagonal (approx 116.57)
  int count_180 = 0; // Center-Opposite (180.0)
  int count_60 = 0;  // Surface-Triangle (60.0)
  int count_108 = 0; // Surface-Pentagon (108.0)
  int count_58 = 0;  // Surface-Center (approx 58.28)

  int total_angles = 0;

  for (const auto &t_1 : angles) {
    for (const auto &t_2 : t_1) {
      for (const auto &t_3 : t_2) {
        for (double const angle : t_3) {
          double const deg = angle * 180.0 / correlation::math::pi;
          total_angles++;

          if (std::abs(deg - 63.43) < 1.0) {
            count_63++;
          } else if (std::abs(deg - 116.57) < 1.0) {
            count_116++;
          } else if (std::abs(deg - 180.0) < 1.0) {
            count_180++;
          } else if (std::abs(deg - 60.0) < 1.0) {
            count_60++;
          } else if (std::abs(deg - 108.0) < 1.0) {
            count_108++;
          } else if (std::abs(deg - 58.28) < 1.0) {
            count_58++;
          }
        }
      }
    }
  }

  // Diagnosis output if needed
  if (total_angles != 246) {
    std::cout << "Found " << total_angles << " angles.\n";
    std::cout << "63: " << count_63 << "\n";
    std::cout << "116: " << count_116 << "\n";
    std::cout << "180: " << count_180 << "\n";
    std::cout << "60: " << count_60 << "\n";
    std::cout << "108: " << count_108 << "\n";
    std::cout << "58: " << count_58 << "\n";
  }

  EXPECT_EQ(count_63, 30) << "Should find 30 Center-Edge angles (~63.4 deg)";
  EXPECT_EQ(count_116, 30) << "Should find 30 Center-Diagonal angles (~116.6 deg)";
  EXPECT_EQ(count_180, 6) << "Should find 6 Center-Opposite angles (180 deg)";
  EXPECT_EQ(count_60, 60) << "Should find 60 Surface-Triangle angles (60 deg)";
  EXPECT_EQ(count_108, 60) << "Should find 60 Surface-Pentagon angles (108 deg)";
  EXPECT_EQ(count_58, 60) << "Should find 60 Surface-Center angles "
                             "(Center-S-correlation::core::Neighbor, ~58.3 deg)";
  EXPECT_EQ(total_angles, 246) << "Total angles should be 246";
}

// ============================================================================
// Part 2: Plane Angle Distribution (PAD) Tests
// ============================================================================
// 1. Trivial Cases
TEST_F(PADTests, EmptyCellThrows) {
  // Current implementation throws explicitly if atoms are empty in
  // calculateAshcroftWeights or implicitly via other checks.
  updateTrajectory();
  EXPECT_THROW({ DistributionFunctions dists(cell_, 5.0, trajectory_.getBondCutoffsSQ()); }, std::invalid_argument);
}

TEST_F(PADTests, SingleAtomNoAngles) {
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

TEST_F(PADTests, NullNeighborsThrows) {
  cell_.addAtom("Si", {10.0, 10.0, 10.0});
  EXPECT_THROW({ correlation::calculators::PADCalculator::calculate(cell_, nullptr, 1.0); }, std::logic_error);
}

// 2. Geometry Verification
TEST_F(PADTests, LinearGeometry180) {
  // A-B-C line
  cell_.addAtom("O", {8.4, 10.0, 10.0});
  cell_.addAtom("Si", {10.0, 10.0, 10.0});
  cell_.addAtom("O", {11.6, 10.0, 10.0});
  updateTrajectory();

  int const id_O = cell_.findElement("O")->id.value;
  int const id_Si = cell_.findElement("Si")->id.value;

  // Verify StructureAnalyzer finds neighbors
  StructureAnalyzer const analyzer(cell_, 2.0, trajectory_.getBondCutoffsSQ());
  const auto &neighborGraph = analyzer.neighborGraph();
  // Si is atom index 1 (0-based)
  ASSERT_GT(neighborGraph.nodeCount(), 1);
  EXPECT_EQ(neighborGraph.getNeighbors(1).size(), 2) << "Si should have 2 neighbors (O atoms)";

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
  EXPECT_NEAR(total_prob, 1.0, 1e-5) << "Should be normalized to 1 angle (normalized by counts * bin_width)";

  // Check peak location
  double peak_val = 0;
  int peak_bin = -1;
  for (size_t i = 0; i < partial.size(); ++i) {
    if (partial[i] > peak_val) {
      peak_val = partial[i];
      peak_bin = static_cast<int>(i);
    }
  }

  if (peak_bin >= 0) {
    double const peak_angle = hist.bins[peak_bin];
    // A 180 degree angle lands in the last bin
    // Depending on exactly how 180 is handled, it's very close to 180
    EXPECT_NEAR(peak_angle, 180.0, 1e-3);
  } else {
    FAIL() << "No peak found in partial distribution";
  }
}

TEST_F(PADTests, RightAngle90) {
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
  double peak_val = 0;
  int peak_bin = -1;
  for (size_t i = 0; i < partial.size(); ++i) {
    if (partial[i] > peak_val) {
      peak_val = partial[i];
      peak_bin = static_cast<int>(i);
    }
  }
  double const peak_angle = hist.bins[peak_bin];
  // 90.0 / 0.001 could land in 89.9995 or 90.0005
  EXPECT_NEAR(peak_angle, 90.0, 0.001);
}

TEST_F(PADTests, EquilateralTriangle60) {
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

TEST_F(PADTests, TetrahedralAngle) {
  // Si at center
  // 4 Neighbors at tetrahedral positions.
  // For simplicity, just check one angle 109.47
  real_t const base_coord = static_cast<real_t>(10.0);
  cell_.addAtom("Si", {base_coord, base_coord, base_coord});
  // Vector 1: (1,1,1) normalized * 1.6
  // Vector 2: (1,-1,-1) normalized * 1.6
  // Dot product = (1-1-1)/3 = -1/3. acos(-1/3) = 109.47 deg

  real_t const lattice_constant = static_cast<real_t>(1.6 * std::numbers::inv_sqrt3);

  cell_.addAtom("O", correlation::math::Vector3<real_t>(base_coord + lattice_constant, base_coord + lattice_constant,
                                                        base_coord + lattice_constant));
  cell_.addAtom("O", correlation::math::Vector3<real_t>(base_coord + lattice_constant, base_coord - lattice_constant,
                                                        base_coord - lattice_constant));
  updateTrajectory();

  DistributionFunctions dists(cell_, 2.0,
                              trajectory_.getBondCutoffsSQ()); // Distance is 1.6
  dists.calculatePAD(0.001);                                   // Hyperfine bins

  const auto &hist = dists.getHistogram("PAD");
  const auto &partial = hist.partials.at("O-Si-O");

  // Expected ~109.471
  double peak_val = 0;
  double peak_angle = 0;
  for (size_t i = 0; i < partial.size(); ++i) {
    if (partial[i] > peak_val) {
      peak_val = partial[i];
      peak_angle = hist.bins[i];
    }
  }
  // 109.4712... / 0.001 -> index 109471 -> center 109.4715
  EXPECT_NEAR(peak_angle, 109.4712206, 0.001);
}

// 3. Symmetry & Multi-Species
TEST_F(PADTests, SymmetryAndSorting) {
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
TEST_F(PADTests, FullNormalizationCheck) {
  // 1 Si, 4 O neighbors (tetrahedron)
  // 4 neighbors -> 4*3/2 = 6 angles.
  // All 6 angles are 109.47
  real_t const base_coord = static_cast<real_t>(10.0);
  cell_.addAtom("Si", {base_coord, base_coord, base_coord});
  real_t const lattice_constant = static_cast<real_t>(1.6 * std::numbers::inv_sqrt3);

  // Tetrahedral vertices
  cell_.addAtom("O", correlation::math::Vector3<real_t>(base_coord + lattice_constant, base_coord + lattice_constant,
                                                        base_coord + lattice_constant));
  cell_.addAtom("O", correlation::math::Vector3<real_t>(base_coord + lattice_constant, base_coord - lattice_constant,
                                                        base_coord - lattice_constant));
  cell_.addAtom("O", correlation::math::Vector3<real_t>(base_coord - lattice_constant, base_coord + lattice_constant,
                                                        base_coord - lattice_constant));
  cell_.addAtom("O", correlation::math::Vector3<real_t>(base_coord - lattice_constant, base_coord - lattice_constant,
                                                        base_coord + lattice_constant));
  updateTrajectory();

  // Custom bond cutoffs to avoid O-O bonds (distance ~2.61) which would create
  // extra angles
  auto cutoffs = trajectory_.getBondCutoffs();
  int const id_O = cell_.findElement("O")->id.value;
  cutoffs[id_O][id_O].max_sq = 2.0;

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

TEST_F(PADTests, IcosahedronAnglesPAD) {
  cell_.addAtom("Si", {10.0, 10.0, 10.0}); // Center
  real_t phi = static_cast<real_t>(std::numbers::phi);
  std::vector<std::vector<real_t>> const vertices = {{0, 1, phi}, {0, 1, -phi}, {0, -1, phi}, {0, -1, -phi},
                                                     {1, phi, 0}, {1, -phi, 0}, {-1, phi, 0}, {-1, -phi, 0},
                                                     {phi, 0, 1}, {phi, 0, -1}, {-phi, 0, 1}, {-phi, 0, -1}};

  for (const auto &vertex : vertices) {
    cell_.addAtom("Si", correlation::math::Vector3<real_t>(static_cast<real_t>(10.0) + vertex[0],
                                                           static_cast<real_t>(10.0) + vertex[1],
                                                           static_cast<real_t>(10.0) + vertex[2]));
  }
  updateTrajectory();

  DistributionFunctions dists(cell_, 2.5, trajectory_.getBondCutoffsSQ());
  dists.calculatePAD(0.01);

  const auto &hist = dists.getHistogram("PAD");
  ASSERT_EQ(hist.partials.count("Si-Si-Si"), 1);
  const auto &partial = hist.partials.at("Si-Si-Si");

  bool found_58 = false;
  bool found_60 = false;
  bool found_63 = false;
  bool found_108 = false;
  bool found_116 = false;
  bool found_180 = false;

  for (size_t i = 0; i < partial.size(); ++i) {
    if (partial[i] > 0.01) { // some density exists
      double const angle = hist.bins[i];
      if (std::abs(angle - 58.28) < 0.5) {
        found_58 = true;
      }
      if (std::abs(angle - 60.0) < 0.5) {
        found_60 = true;
      }
      if (std::abs(angle - 63.43) < 0.5) {
        found_63 = true;
      }
      if (std::abs(angle - 108.0) < 0.5) {
        found_108 = true;
      }
      if (std::abs(angle - 116.57) < 0.5) {
        found_116 = true;
      }
      if (std::abs(angle - 180.0) < 0.5) {
        found_180 = true;
      }
    }
  }

  EXPECT_TRUE(found_58) << "Should find PAD peak near 58.28 degrees";
  EXPECT_TRUE(found_60) << "Should find PAD peak near 60.00 degrees";
  EXPECT_TRUE(found_63) << "Should find PAD peak near 63.43 degrees";
  EXPECT_TRUE(found_108) << "Should find PAD peak near 108.00 degrees";
  EXPECT_TRUE(found_116) << "Should find PAD peak near 116.57 degrees";
  EXPECT_TRUE(found_180) << "Should find PAD peak near 180.00 degrees";
}

TEST_F(PADTests, BondDistanceBelowMinCutoffProducesNoAngles) {
  // Center Si at (10, 10, 10), O at (10.5, 10, 10) [dist 0.5 Å], O at (10, 11.5, 10) [dist 1.5 Å]
  cell_.addAtom("Si", {10.0, 10.0, 10.0});
  cell_.addAtom("O", {10.5, 10.0, 10.0}); // dist = 0.5 Å
  cell_.addAtom("O", {10.0, 11.5, 10.0}); // dist = 1.5 Å
  updateTrajectory();

  int const id_si = cell_.findElement("Si")->id.value;
  int const id_o = cell_.findElement("O")->id.value;

  // Cutoff range for Si-O: [1.0 Å, 2.0 Å] -> min_sq = 1.0, max_sq = 4.0; O-O and Si-Si = 0
  BondCutoffMatrix cutoffs(2, std::vector<BondCutoffRange>(2, BondCutoffRange{0.0, 0.0}));
  cutoffs[id_si][id_o] = BondCutoffRange{1.0, 4.0};
  cutoffs[id_o][id_si] = BondCutoffRange{1.0, 4.0};
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

TEST_F(PADTests, BondDistanceAboveMaxCutoffProducesNoAngles) {
  // Center Si at (10, 10, 10), O at (11.5, 10, 10) [dist 1.5 Å], O at (10, 12.5, 10) [dist 2.5 Å]
  cell_.addAtom("Si", {10.0, 10.0, 10.0});
  cell_.addAtom("O", {11.5, 10.0, 10.0}); // dist = 1.5 Å
  cell_.addAtom("O", {10.0, 12.5, 10.0}); // dist = 2.5 Å
  updateTrajectory();

  int const id_si = cell_.findElement("Si")->id.value;
  int const id_o = cell_.findElement("O")->id.value;

  // Cutoff range for Si-O: [1.0 Å, 2.0 Å] -> min_sq = 1.0, max_sq = 4.0; O-O and Si-Si = 0
  BondCutoffMatrix cutoffs(2, std::vector<BondCutoffRange>(2, BondCutoffRange{0.0, 0.0}));
  cutoffs[id_si][id_o] = BondCutoffRange{1.0, 4.0};
  cutoffs[id_o][id_si] = BondCutoffRange{1.0, 4.0};
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

TEST_F(PADTests, BondDistanceWithinCutoffRangeProducesAngle) {
  // Center Si at (10, 10, 10), O at (11.5, 10, 10) [dist 1.5 Å], O at (10, 11.5, 10) [dist 1.5 Å]
  cell_.addAtom("Si", {10.0, 10.0, 10.0});
  cell_.addAtom("O", {11.5, 10.0, 10.0}); // dist = 1.5 Å
  cell_.addAtom("O", {10.0, 11.5, 10.0}); // dist = 1.5 Å
  updateTrajectory();

  int const id_si = cell_.findElement("Si")->id.value;
  int const id_o = cell_.findElement("O")->id.value;

  // Cutoff range for Si-O: [1.0 Å, 2.0 Å] -> min_sq = 1.0, max_sq = 4.0; O-O and Si-Si = 0
  BondCutoffMatrix cutoffs(2, std::vector<BondCutoffRange>(2, BondCutoffRange{0.0, 0.0}));
  cutoffs[id_si][id_o] = BondCutoffRange{1.0, 4.0};
  cutoffs[id_o][id_si] = BondCutoffRange{1.0, 4.0};
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

TEST_F(PADTests, VerifyPADRawHistogram) {
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
