// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/Atom.hpp"
#include "../include/Cell.hpp"
#include "../include/PhysicalData.hpp"
#include "../include/StructureAnalyzer.hpp"
#include "../include/Trajectory.hpp"
#include <gtest/gtest.h>
#include <vector>

// A test fixture for StructureAnalyzer tests.
// A test fixture for StructureAnalyzer tests.
class Test04_StructureAnalyzer : public ::testing::Test {
protected:
  Trajectory trajectory_;

  void updateTrajectory(const Cell &cell) {
    trajectory_ = Trajectory();
    trajectory_.addFrame(cell);
    trajectory_.precomputeBondCutoffs();
  }
};

TEST_F(Test04_StructureAnalyzer, FindsCorrectNeighborsForSilicon) {
  // ... (existing test content)
  // Arrange: Create an 8-atom conventional unit cell of Silicon.
  // The diamond lattice structure is a robust test for neighbor finding.
  const double lattice_const = 5.43; // Angstroms
  Cell si_cell({lattice_const, lattice_const, lattice_const, 90.0, 90.0, 90.0});

  // Fractional coordinates for the 8 atoms in a diamond cubic cell
  std::vector<linalg::Vector3<double>> fractional_coords = {
      {0.0, 0.0, 0.0},    {0.5, 0.5, 0.0},    {0.5, 0.0, 0.5},
      {0.0, 0.5, 0.5},    {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
      {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};

  for (const auto &frac_pos : fractional_coords) {
    // Convert fractional to Cartesian coordinates before adding to the cell
    si_cell.addAtom("Si", si_cell.latticeVectors() * frac_pos);
  }
  updateTrajectory(si_cell);

  // Act: Calculate neighbors with a cutoff just beyond the first neighbor
  // shell. The first nearest neighbor distance in Si is (sqrt(3)/4)*a
  // approx 2.35 Å.
  StructureAnalyzer neighbors(si_cell, 3.0, trajectory_.getBondCutoffsSQ());
  const auto &neighborMatrix = neighbors.neighbors();
  const auto &atoms = si_cell.atoms();

  // Assert
  ASSERT_EQ(neighborMatrix.size(), 8);
  const double expected_distance = 2.3512; // More precise value

  // Every Si atom in a diamond lattice must have exactly 4 nearest neighbors.
  for (const auto neighbors : neighborMatrix) {
    ASSERT_EQ(neighbors.size(), 4);
    // Check that the distance to each of these neighbors is correct.
    for (const auto neighbor : neighbors) {
      EXPECT_NEAR(neighbor.distance, expected_distance, 1e-4);
    }
  }
}

TEST_F(Test04_StructureAnalyzer, DistancesTensorIsCorrect) {
  // Arrange: Create a simple 2-atom system
  Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {0.0, 0.0, 0.0});
  cell.addAtom("Ar", {3.0, 0.0, 0.0});
  updateTrajectory(cell);

  // Act: Calculate distances
  // Cutoff 5.0 covers the 3.0 distance
  StructureAnalyzer analyzer(cell, 5.0, trajectory_.getBondCutoffsSQ());
  const auto &distances = analyzer.distances();

  // Assert Same Species
  int id_Ar = cell.findElement("Ar")->id.value;
  ASSERT_EQ(id_Ar, 0);

  // StructureAnalyzer stores unique pairs for same-species (i < j logic)
  // So for 2 Ar atoms, we expect 1 distance in [Ar][Ar]
  ASSERT_GT(distances.size(), 0);
  const auto &ar_ar_dists = distances[id_Ar][id_Ar];
  EXPECT_EQ(ar_ar_dists.size(), 1);
  EXPECT_NEAR(ar_ar_dists[0], 3.0, 1e-6);

  // Arrange Mixed Species
  Cell mixed_cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  mixed_cell.addAtom("Ar", {0.0, 0.0, 0.0});
  mixed_cell.addAtom("Xe", {0.0, 4.0, 0.0});
  updateTrajectory(mixed_cell);

  StructureAnalyzer mixed_analyzer(mixed_cell, 5.0,
                                   trajectory_.getBondCutoffsSQ());
  const auto &mixed_distances = mixed_analyzer.distances();

  int id_Ar_m = mixed_cell.findElement("Ar")->id.value;
  int id_Xe_m = mixed_cell.findElement("Xe")->id.value;

  // StructureAnalyzer stores symmetric pairs for different species
  // So we expect 1 in [Ar][Xe] and 1 in [Xe][Ar]
  const auto &ar_xe_dists = mixed_distances[id_Ar_m][id_Xe_m];
  const auto &xe_ar_dists = mixed_distances[id_Xe_m][id_Ar_m];

  EXPECT_EQ(ar_xe_dists.size(), 1);
  EXPECT_EQ(xe_ar_dists.size(), 1);
  EXPECT_NEAR(ar_xe_dists[0], 4.0, 1e-6);
  EXPECT_NEAR(xe_ar_dists[0], 4.0, 1e-6);
}

TEST_F(Test04_StructureAnalyzer, CalculatesCorrectAnglesForWater) {
  // Arrange: Create a single water molecule with a known bond angle.
  Cell water_cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  const double bond_length = 0.957;    // Angstroms
  const double bond_angle_deg = 104.5; // Degrees
  const double bond_angle_rad = bond_angle_deg * constants::deg2rad;

  water_cell.addAtom("O", {10.0, 10.0, 10.0});
  water_cell.addAtom("H", {10.0 + bond_length, 10.0, 10.0});
  water_cell.addAtom("H",
                     {10.0 + bond_length * std::cos(bond_angle_rad),
                      10.0 + bond_length * std::sin(bond_angle_rad), 10.0});
  updateTrajectory(water_cell);

  // Act: Calculate neighbors and angles.
  StructureAnalyzer neighbors(water_cell, 3.0, trajectory_.getBondCutoffsSQ());
  const auto &angles = neighbors.angles();

  // Assert: We need to get the element IDs to index the angle tensor correctly.
  const int h_id = water_cell.findElement("H")->id.value;
  const int o_id = water_cell.findElement("O")->id.value;
  ASSERT_EQ(h_id, 1);
  ASSERT_EQ(o_id, 0);

  // The angle tensor is indexed by [type1][central_type][type2].
  ASSERT_TRUE(angles.size() > h_id && angles[h_id].size() > o_id &&
              angles[h_id][o_id].size() > h_id);

  const auto &hoh_angles = angles[h_id][o_id][h_id];

  // There should be exactly one H-O-H angle in this system.
  ASSERT_EQ(hoh_angles.size(), 1);
  // The calculated angle should match the known value.
  EXPECT_NEAR(hoh_angles[0] * constants::rad2deg, bond_angle_deg, 1e-4);
}

TEST_F(Test04_StructureAnalyzer, CalculatesCorrectAngleWithPBC) {
  // Arrange: Setup a system where the angle calculation requires the Minimum
  // Image Convention. Central atom B (index 1) at (0.5, 0.5, 0.5). Neighbor A
  // (index 0) at (3.5, 0.5, 0.5) -> PBC vector B->A is (-1.0, 0.0, 0.0).
  // Neighbor C (index 2) at (0.5, 3.5, 0.5) -> PBC vector B->C is (0.0, -1.0,
  // 0.0). The resulting angle A-B-C must be 90 degrees (pi/2 radians).

  const double side_length = 10.0;
  const double cutoff = 2.0;
  // Calculate pi/2 explicitly using acos(-1.0) = pi.
  const double expected_angle_rad = std::acos(-1.0) / 2.0;

  Cell pbc_cell({side_length, side_length, side_length, 90.0, 90.0, 90.0});

  // B (Central) at (0.5, 0.5, 0.5)
  // A at (9.0, 0.5, 0.5) -> Wrapped dist to B is 1.5
  // C at (0.5, 9.0, 0.5) -> Wrapped dist to B is 1.5
  // Dist A-C is sqrt(1.5^2 + 1.5^2) = 2.12 > bond_cutoff (~1.85)

  pbc_cell.addAtom("C", {9.0, 0.5, 0.5}); // Atom A
  pbc_cell.addAtom("C", {0.5, 0.5, 0.5}); // Atom B (Central)
  pbc_cell.addAtom("O", {0.5, 9.0, 0.5}); // Atom C
  updateTrajectory(pbc_cell);

  // Act: Calculate neighbors and angles.
  StructureAnalyzer analyzer(pbc_cell, cutoff, trajectory_.getBondCutoffsSQ());
  const auto &angles = analyzer.angles();

  // Assert
  const int c_id = pbc_cell.findElement("C")->id.value;
  const int o_id = pbc_cell.findElement("O")->id.value;
  ASSERT_EQ(c_id, 0);
  ASSERT_EQ(o_id, 1);

  // The angle is stored in the [C][C][O] slot.
  const auto &cco_angles = angles[c_id][c_id][o_id];

  // There must be exactly one angle calculated (A-B-C).
  ASSERT_EQ(cco_angles.size(), 1);

  // The calculated angle should be pi/2 radians.
  EXPECT_NEAR(cco_angles[0], expected_angle_rad, 1e-6);
}

TEST_F(Test04_StructureAnalyzer, FindsNoNeighborsForIsolatedAtom) {
  // Arrange: Create a large cell with a single atom.
  Cell large_cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  large_cell.addAtom("Ar", {10.0, 10.0, 10.0});
  updateTrajectory(large_cell);
  // Use default bond factor.

  // Act: Calculate neighbors with a moderate cutoff.
  StructureAnalyzer analyzer(large_cell, 5.0, trajectory_.getBondCutoffsSQ());
  const auto &neighbors = analyzer.neighbors();

  // Assert: The single atom should have no neighbors.
  ASSERT_FALSE(neighbors.empty()); // The vector itself has size 1 (for 1 atom)
  ASSERT_TRUE(neighbors[0].empty());
}

TEST_F(Test04_StructureAnalyzer, FindsNeighborsBasedOnBondCutoff) {
  // Ar covalent radius = 0.96. Sum = 1.92.
  // Bond threshold = 1.92 * 1.2 = 2.304.

  // Case 1: Within bond threshold -> Neighbors Found
  {
    Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
    cell.addAtom("Ar", {5.0, 5.0, 5.0});
    cell.addAtom("Ar", {7.0, 5.0, 5.0}); // Distance = 2.0 < 2.304
    updateTrajectory(cell);
    StructureAnalyzer analyzer(cell, 3.1, trajectory_.getBondCutoffsSQ());
    const auto &neighbors = analyzer.neighbors();
    ASSERT_EQ(neighbors.size(), 2);
    ASSERT_EQ(neighbors[0].size(), 1);
    EXPECT_EQ(neighbors[0][0].index, 1);
    EXPECT_NEAR(neighbors[0][0].distance, 2.0, 1e-6);
  }

  // Case 2: Out of bond threshold -> No Neighbors
  {
    Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
    cell.addAtom("Ar", {5.0, 5.0, 5.0});
    cell.addAtom("Ar", {8.0, 5.0, 5.0}); // Distance = 3.0 > 2.304
    updateTrajectory(cell);
    StructureAnalyzer analyzer(cell, 3.5, trajectory_.getBondCutoffsSQ());
    const auto &neighbors = analyzer.neighbors();
    ASSERT_EQ(neighbors.size(), 2);
    EXPECT_TRUE(neighbors[0].empty());
    EXPECT_TRUE(neighbors[1].empty());
  }
}

TEST_F(Test04_StructureAnalyzer, EnforcesNeighborSymmetry) {
  // Arrange: Create a random system of atoms.
  Cell random_cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  random_cell.addAtom("Ar", {1.0, 1.0, 1.0}); // A
  random_cell.addAtom("Ar", {2.0, 2.0, 2.0}); // B (Dist ~1.73)
  random_cell.addAtom("Ar", {8.0, 8.0, 8.0}); // C
  random_cell.addAtom("Ar", {1.5, 1.5, 1.5}); // D (Dist to A ~0.866)
  updateTrajectory(random_cell);

  // Ar radius 0.96. Sum = 1.92.
  // Dist A-B = 1.73 < 1.92. Should be bonded with default bond_factor 1.2?
  // 1.92 * 1.2 = 2.304.
  // 1.92 * 1.2 = 2.304.
  // So A-B should be neighbors.
  // A-D should be neighbors.
  // Act
  StructureAnalyzer analyzer(random_cell, 3.0, trajectory_.getBondCutoffsSQ());
  const auto &neighbors = analyzer.neighbors();

  // Assert: Check symmetry for all pairs.
  for (size_t i = 0; i < neighbors.size(); ++i) {
    for (const auto &neighbor : neighbors[i]) {
      size_t j = neighbor.index;

      // If i sees j, then j must see i
      bool found_reverse = false;
      for (const auto &reverse_neighbor : neighbors[j]) {
        if (reverse_neighbor.index == i) {
          found_reverse = true;
          // Distances must match
          EXPECT_NEAR(reverse_neighbor.distance, neighbor.distance, 1e-6);
          // Vectors must be opposite (r_ij = j - i)
          EXPECT_NEAR(reverse_neighbor.r_ij.x(), -neighbor.r_ij.x(), 1e-6);
          EXPECT_NEAR(reverse_neighbor.r_ij.y(), -neighbor.r_ij.y(), 1e-6);
          EXPECT_NEAR(reverse_neighbor.r_ij.z(), -neighbor.r_ij.z(), 1e-6);
          break;
        }
      }
      EXPECT_TRUE(found_reverse) << "Symmetry broken: Atom " << i << " sees "
                                 << j << " but " << j << " does not see " << i;
    }
  }
}

TEST_F(Test04_StructureAnalyzer, HandlesPeriodicSelfInteractions) {
  // 1.92 * 1.2 = 2.304.
  // We need images to be within the bond threshold. Let's use image distance
  // = 2.0.

  Cell small_cell({2.0, 2.0, 2.0, 90.0, 90.0, 90.0});
  small_cell.addAtom("Ar", {1.0, 1.0, 1.0});
  updateTrajectory(small_cell);

  // Case 1: Ignore Periodic Self Interactions = true (Reference default)
  {
    StructureAnalyzer analyzer(small_cell, 3.0, trajectory_.getBondCutoffsSQ(),
                               true);
    const auto &neighbors = analyzer.neighbors();
    ASSERT_EQ(neighbors.size(), 1);
    // Should NOT see itself
    EXPECT_TRUE(neighbors[0].empty());
  }

  // Case 2: Ignore Periodic Self Interactions = false
  {
    StructureAnalyzer analyzer(small_cell, 3.0, trajectory_.getBondCutoffsSQ(),
                               false);
    const auto &neighbors = analyzer.neighbors();
    ASSERT_EQ(neighbors.size(), 1);

    // In a cubic cell with side 2.0, nearest images are 6 face-sharing images
    // at d=2.0. 2.0 < 2.304, so they should be found as bonded neighbors.
    EXPECT_EQ(neighbors[0].size(), 6);
    for (const auto &n : neighbors[0]) {
      EXPECT_NEAR(n.distance, 2.0, 1e-6);
    }
  }
}

TEST_F(Test04_StructureAnalyzer, TriangleMoleculeConnectivity) {
  // Arrange: Create 3 atoms in an equilateral triangle.
  Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});

  // Equilateral triangle with side length 2.0
  // Covalent radius Ar = 0.96 -> Bond cutoff = 1.92 * 1.2 = 2.304.
  // Distance 2.0 < 2.304, so they should be connected.

  const double d = 2.0;
  const double h = std::sqrt(d * d - (d / 2) * (d / 2)); // sqrt(3)

  cell.addAtom("Ar", {5.0, 5.0, 5.0});             // A
  cell.addAtom("Ar", {5.0 + d, 5.0, 5.0});         // B
  cell.addAtom("Ar", {5.0 + d / 2, 5.0 + h, 5.0}); // C
  updateTrajectory(cell);

  // Act
  StructureAnalyzer analyzer(cell, 3.0, trajectory_.getBondCutoffsSQ());
  const auto &neighbors = analyzer.neighbors();

  // Assert
  ASSERT_EQ(neighbors.size(), 3);

  // Each atom should see exactly 2 neighbors (the other two vertices)
  for (size_t i = 0; i < 3; ++i) {
    ASSERT_EQ(neighbors[i].size(), 2)
        << "Atom " << i << " does not have 2 neighbors";
    for (const auto &n : neighbors[i]) {
      EXPECT_NEAR(n.distance, d, 1e-6);
    }
  }
}
