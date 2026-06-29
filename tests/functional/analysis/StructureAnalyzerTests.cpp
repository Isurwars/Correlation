// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/StructureAnalyzer.hpp"
#include "core/Atom.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::analysis {
namespace {

// A test fixture for StructureAnalyzer tests.
class StructureAnalyzerTests : public ::testing::Test {
protected:
  void updateTrajectory(const correlation::core::Cell &cell) {
    trajectory_ = correlation::core::Trajectory();
    trajectory_.addFrame(cell);
    trajectory_.precomputeBondCutoffs();
  }

  [[nodiscard]] const correlation::core::Trajectory &trajectory() const { return trajectory_; }

private:
  correlation::core::Trajectory trajectory_;
};

TEST_F(StructureAnalyzerTests, FindsCorrectNeighborsForSilicon) {
  // Arrange: Create an 8-atom conventional unit cell of Silicon.
  // The diamond lattice structure is a robust test for neighbor finding.
  const double lattice_const = 5.43; // Angstroms
  correlation::core::Cell si_cell({lattice_const, lattice_const, lattice_const, 90.0, 90.0, 90.0});

  // Fractional coordinates for the 8 atoms in a diamond cubic cell
  std::vector<correlation::math::Vector3<double>> const fractional_coords = {
      {0.0, 0.0, 0.0},    {0.5, 0.5, 0.0},    {0.5, 0.0, 0.5},    {0.0, 0.5, 0.5},
      {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25}, {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};

  for (const auto &frac_pos : fractional_coords) {
    // Convert fractional to Cartesian coordinates before adding to the cell
    si_cell.addAtom("Si", si_cell.latticeVectors() * frac_pos);
  }
  updateTrajectory(si_cell);

  // Act: Calculate neighbors with a cutoff just beyond the first neighbor
  // shell. The first nearest neighbor distance in Si is (sqrt(3)/4)*a
  // approx 2.35 Å.
  StructureAnalyzer const neighbors(si_cell, 3.0, trajectory().getBondCutoffsSQ());
  const auto &neighborGraph = neighbors.neighborGraph();
  const auto &atoms = si_cell.atoms();

  // Assert
  ASSERT_EQ(neighborGraph.nodeCount(), 8);
  const double expected_distance = 2.3512; // More precise value

  // Every Si atom in a diamond lattice must have exactly 4 nearest neighbors.
  for (size_t i = 0; i < neighborGraph.nodeCount(); ++i) {
    const auto &node_neighbors = neighborGraph.getNeighbors(i);
    ASSERT_EQ(node_neighbors.size(), 4);
    // Check that the distance to each of these neighbors is correct.
    for (const auto neighbor : node_neighbors) {
      EXPECT_NEAR(neighbor.distance, expected_distance, 1e-4);
    }
  }
}

TEST_F(StructureAnalyzerTests, DistancesTensorIsCorrect) {
  // Arrange: Create a simple 2-atom system
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {0.0, 0.0, 0.0});
  cell.addAtom("Ar", {3.0, 0.0, 0.0});
  updateTrajectory(cell);

  // Act: Calculate distances
  // Cutoff 5.0 covers the 3.0 distance
  StructureAnalyzer const analyzer(cell, 5.0, trajectory().getBondCutoffsSQ());
  const auto &distances = analyzer.distances();

  // Assert Same Species
  int const id_Ar = cell.findElement("Ar")->id.value;
  ASSERT_EQ(id_Ar, 0);

  // StructureAnalyzer stores unique pairs for same-species (i < j logic)
  // So for 2 Ar atoms, we expect 1 distance in [Ar][Ar]
  ASSERT_GT(distances.size(), 0);
  const auto &ar_ar_dists = distances[id_Ar][id_Ar];
  EXPECT_EQ(ar_ar_dists.size(), 1);
  EXPECT_NEAR(ar_ar_dists[0], 3.0, 1e-6);

  // Arrange Mixed Species
  correlation::core::Cell mixed_cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  mixed_cell.addAtom("Ar", {0.0, 0.0, 0.0});
  mixed_cell.addAtom("Xe", {0.0, 4.0, 0.0});
  updateTrajectory(mixed_cell);

  StructureAnalyzer const mixed_analyzer(mixed_cell, 5.0, trajectory().getBondCutoffsSQ());
  const auto &mixed_distances = mixed_analyzer.distances();

  int const id_Ar_m = mixed_cell.findElement("Ar")->id.value;
  int const id_Xe_m = mixed_cell.findElement("Xe")->id.value;

  // StructureAnalyzer stores symmetric pairs for different species
  // So we expect 1 in [Ar][Xe] and 1 in [Xe][Ar]
  const auto &ar_xe_dists = mixed_distances[id_Ar_m][id_Xe_m];
  const auto &xe_ar_dists = mixed_distances[id_Xe_m][id_Ar_m];

  EXPECT_EQ(ar_xe_dists.size(), 1);
  EXPECT_EQ(xe_ar_dists.size(), 1);
  EXPECT_NEAR(ar_xe_dists[0], 4.0, 1e-6);
  EXPECT_NEAR(xe_ar_dists[0], 4.0, 1e-6);
}

TEST_F(StructureAnalyzerTests, CalculatesCorrectAnglesForWater) {
  // Arrange: Create a single water molecule with a known bond angle.
  correlation::core::Cell water_cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  const double bond_length = 0.957;    // Angstroms
  const double bond_angle_deg = 104.5; // Degrees
  const double bond_angle_rad = bond_angle_deg * correlation::math::deg_to_rad;

  water_cell.addAtom("O", {10.0, 10.0, 10.0});
  water_cell.addAtom("H", {10.0 + bond_length, 10.0, 10.0});
  water_cell.addAtom(
      "H", {10.0 + bond_length * std::cos(bond_angle_rad), 10.0 + bond_length * std::sin(bond_angle_rad), 10.0});
  updateTrajectory(water_cell);

  // Act: Calculate neighbors and angles.
  StructureAnalyzer const neighbors(water_cell, 3.0, trajectory().getBondCutoffsSQ());
  const auto &angles = neighbors.angles();

  // Assert: We need to get the element IDs to index the angle tensor correctly.
  const int h_id = water_cell.findElement("H")->id.value;
  const int o_id = water_cell.findElement("O")->id.value;
  ASSERT_EQ(h_id, 1);
  ASSERT_EQ(o_id, 0);

  // The angle tensor is indexed by [type1][central_type][type2].
  ASSERT_TRUE(angles.size() > h_id && angles[h_id].size() > o_id && angles[h_id][o_id].size() > h_id);

  const auto &hoh_angles = angles[h_id][o_id][h_id];

  // There should be exactly one H-O-H angle in this system.
  ASSERT_EQ(hoh_angles.size(), 1);
  // The calculated angle should match the known value.
  EXPECT_NEAR(hoh_angles[0] * correlation::math::rad_to_deg, bond_angle_deg, 1e-4);
}

TEST_F(StructureAnalyzerTests, CalculatesCorrectAngleWithPBC) {
  // Arrange: Setup a system where the angle calculation requires the Minimum
  // Image Convention. Central atom B (index 1) at (0.5, 0.5, 0.5).
  // correlation::core::Neighbor A (index 0) at (3.5, 0.5, 0.5) -> PBC vector
  // B->A is (-1.0, 0.0, 0.0). correlation::core::Neighbor C (index 2) at
  // (0.5, 3.5, 0.5) -> PBC vector B->C is (0.0, -1.0, 0.0). The resulting angle
  // A-B-C must be 90 degrees (pi/2 radians).

  const double side_length = 10.0;
  const double cutoff = 2.0;
  // Calculate pi/2 explicitly using acos(-1.0) = pi.
  const double expected_angle_rad = std::acos(-1.0) / 2.0;

  correlation::core::Cell pbc_cell({side_length, side_length, side_length, 90.0, 90.0, 90.0});

  // B (Central) at (0.5, 0.5, 0.5)
  // A at (9.0, 0.5, 0.5) -> Wrapped dist to B is 1.5
  // C at (0.5, 9.0, 0.5) -> Wrapped dist to B is 1.5
  // Dist A-C is sqrt(1.5^2 + 1.5^2) = 2.12 > bond_cutoff (~1.85)

  pbc_cell.addAtom("C", {9.0, 0.5, 0.5}); // correlation::core::Atom A
  pbc_cell.addAtom("C", {0.5, 0.5, 0.5}); // correlation::core::Atom B (Central)
  pbc_cell.addAtom("O", {0.5, 9.0, 0.5}); // correlation::core::Atom C
  updateTrajectory(pbc_cell);

  // Act: Calculate neighbors and angles.
  StructureAnalyzer const analyzer(pbc_cell, cutoff, trajectory().getBondCutoffsSQ());
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

TEST_F(StructureAnalyzerTests, FindsNoNeighborsForIsolatedAtom) {
  // Arrange: Create a large cell with a single atom.
  correlation::core::Cell large_cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  large_cell.addAtom("Ar", {10.0, 10.0, 10.0});
  updateTrajectory(large_cell);
  // Use default bond factor.

  // Act: Calculate neighbors with a moderate cutoff.
  StructureAnalyzer const analyzer(large_cell, 5.0, trajectory().getBondCutoffsSQ());
  const auto &neighborGraph = analyzer.neighborGraph();

  // Assert: The single atom should have no neighbors.
  ASSERT_GT(neighborGraph.nodeCount(),
            0); // The vector itself has size 1 (for 1 atom)
  ASSERT_TRUE(neighborGraph.getNeighbors(0).empty());
}

TEST_F(StructureAnalyzerTests, FindsNeighborsBasedOnBondCutoff) {
  // Ar covalent radius = 0.96. Sum = 1.92.
  // Bond threshold = 1.92 * 1.2 = 2.304.

  // Case 1: Within bond threshold -> Neighbors Found
  {
    correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
    cell.addAtom("Ar", {5.0, 5.0, 5.0});
    cell.addAtom("Ar", {7.0, 5.0, 5.0}); // Distance = 2.0 < 2.304
    updateTrajectory(cell);
    StructureAnalyzer const analyzer(cell, 3.1, trajectory().getBondCutoffsSQ());
    const auto &neighborGraph = analyzer.neighborGraph();
    ASSERT_EQ(neighborGraph.nodeCount(), 2);
    ASSERT_EQ(neighborGraph.getNeighbors(0).size(), 1);
    EXPECT_EQ(neighborGraph.getNeighbors(0)[0].index, 1);
    EXPECT_NEAR(neighborGraph.getNeighbors(0)[0].distance, 2.0, 1e-6);
  }

  // Case 2: Out of bond threshold -> No Neighbors
  {
    correlation::core::Cell cell({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
    cell.addAtom("Ar", {5.0, 5.0, 5.0});
    cell.addAtom("Ar", {8.0, 5.0, 5.0}); // Distance = 3.0 > 2.304
    updateTrajectory(cell);
    StructureAnalyzer const analyzer(cell, 3.5, trajectory().getBondCutoffsSQ());
    const auto &neighborGraph = analyzer.neighborGraph();
    ASSERT_EQ(neighborGraph.nodeCount(), 2);
    EXPECT_TRUE(neighborGraph.getNeighbors(0).empty());
    EXPECT_TRUE(neighborGraph.getNeighbors(1).empty());
  }
}

TEST_F(StructureAnalyzerTests, EnforcesNeighborSymmetry) {
  // Arrange: Create a random system of atoms.
  correlation::core::Cell random_cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
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
  StructureAnalyzer const analyzer(random_cell, 3.0, trajectory().getBondCutoffsSQ());
  const auto &neighborGraph = analyzer.neighborGraph();

  // Assert: Check symmetry for all pairs.
  for (size_t i = 0; i < neighborGraph.nodeCount(); ++i) {
    for (const auto &neighbor : neighborGraph.getNeighbors(i)) {
      size_t const neighbor_idx = neighbor.index;

      // If i sees j, then j must see i
      bool found_reverse = false;
      for (const auto &reverse_neighbor : neighborGraph.getNeighbors(neighbor_idx)) {
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
      EXPECT_TRUE(found_reverse) << "Symmetry broken: correlation::core::Atom " << i << " sees " << neighbor_idx
                                 << " but " << neighbor_idx << " does not see " << i;
    }
  }
}

TEST_F(StructureAnalyzerTests, HandlesPeriodicSelfInteractions) {
  // 1.92 * 1.2 = 2.304.
  // We need images to be within the bond threshold. Let's use image distance
  // = 2.0.

  correlation::core::Cell small_cell({2.0, 2.0, 2.0, 90.0, 90.0, 90.0});
  small_cell.addAtom("Ar", {1.0, 1.0, 1.0});
  updateTrajectory(small_cell);

  // Case 1: Ignore Periodic Self Interactions = true (Reference default)
  {
    StructureAnalyzer const analyzer(small_cell, 3.0, trajectory().getBondCutoffsSQ(), true);
    const auto &neighborGraph = analyzer.neighborGraph();
    ASSERT_EQ(neighborGraph.nodeCount(), 1);
    // Should NOT see itself
    EXPECT_TRUE(neighborGraph.getNeighbors(0).empty());
  }

  // Case 2: Ignore Periodic Self Interactions = false
  {
    StructureAnalyzer const analyzer(small_cell, 3.0, trajectory().getBondCutoffsSQ(), false);
    const auto &neighborGraph = analyzer.neighborGraph();
    ASSERT_EQ(neighborGraph.nodeCount(), 1);

    // In a cubic cell with side 2.0, nearest images are 6 face-sharing images
    // at d=2.0. 2.0 < 2.304, so they should be found as bonded neighbors.
    EXPECT_EQ(neighborGraph.getNeighbors(0).size(), 6);
    for (const auto &neighbor : neighborGraph.getNeighbors(0)) {
      EXPECT_NEAR(neighbor.distance, 2.0, 1e-6);
    }
  }
}

TEST_F(StructureAnalyzerTests, TriangleMoleculeConnectivity) {
  // Arrange: Create 3 atoms in an equilateral triangle.
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});

  // Equilateral triangle with side length 2.0
  // Covalent radius Ar = 0.96 -> Bond cutoff = 1.92 * 1.2 = 2.304.
  // Distance 2.0 < 2.304, so they should be connected.

  const double distance = 2.0;
  const double height_triangle = std::sqrt(distance * distance - (distance / 2) * (distance / 2)); // sqrt(3)

  cell.addAtom("Ar", {5.0, 5.0, 5.0});                                  // A
  cell.addAtom("Ar", {5.0 + distance, 5.0, 5.0});                       // B
  cell.addAtom("Ar", {5.0 + distance / 2, 5.0 + height_triangle, 5.0}); // C
  updateTrajectory(cell);

  // Act
  StructureAnalyzer const analyzer(cell, 3.0, trajectory().getBondCutoffsSQ());
  const auto &neighborGraph = analyzer.neighborGraph();

  // Assert
  ASSERT_EQ(neighborGraph.nodeCount(), 3);

  // Each atom should see exactly 2 neighbors (the other two vertices)
  for (size_t i = 0; i < 3; ++i) {
    ASSERT_EQ(neighborGraph.getNeighbors(i).size(), 2)
        << "correlation::core::Atom " << i << " does not have 2 neighbors";
    for (const auto &neighbor : neighborGraph.getNeighbors(i)) {
      EXPECT_NEAR(neighbor.distance, distance, 1e-6);
    }
  }
}

TEST_F(StructureAnalyzerTests, CalculatesCorrectDihedralAngles) {
  // Arrange: Create a 4-atom chain A-B-C-D with known dihedral of -pi/2
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});

  // A at (0, 0, 0), B at (1.7, 0, 0), C at (1.7, 1.7, 0), D at (1.7, 1.7, 1.7)
  // Distance A-B = 1.7, B-C = 1.7, C-D = 1.7
  cell.addAtom("C", {0.0, 0.0, 0.0}); // A (index 0)
  cell.addAtom("C", {1.7, 0.0, 0.0}); // B (index 1)
  cell.addAtom("C", {1.7, 1.7, 0.0}); // C (index 2)
  cell.addAtom("C", {1.7, 1.7, 1.7}); // D (index 3)
  updateTrajectory(cell);

  // Act: Compute structure analysis
  // Covalent radius C is ~0.77. Sum is 1.54. Bond cutoff with 1.2 factor is 1.848.
  // 1.7 < 1.848, so A-B, B-C, C-D are bonded.
  // Non-adjacent distances are 2.404 and 2.944, which are > 1.848, so no other bonds.
  StructureAnalyzer const analyzer(cell, 3.0, trajectory().getBondCutoffsSQ());
  const auto &dihedrals = analyzer.dihedrals();

  // Assert
  const int carbon_id = cell.findElement("C")->id.value;
  ASSERT_EQ(carbon_id, 0);

  // The dihedral tensor is indexed by [typeA][typeB][typeC][typeD]
  ASSERT_GT(dihedrals.size(), carbon_id);
  ASSERT_GT(dihedrals[carbon_id].size(), carbon_id);
  ASSERT_GT(dihedrals[carbon_id][carbon_id].size(), carbon_id);
  ASSERT_GT(dihedrals[carbon_id][carbon_id][carbon_id].size(), carbon_id);

  const auto &c_c_c_c_dihedrals = dihedrals[carbon_id][carbon_id][carbon_id][carbon_id];
  // Since we have a single chain, we expect exactly 1 unique dihedral angle (A-B-C-D)
  ASSERT_EQ(c_c_c_c_dihedrals.size(), 1);
  EXPECT_NEAR(c_c_c_c_dihedrals[0], -std::numbers::pi / 2.0, 1e-6);
}

TEST_F(StructureAnalyzerTests, ThrowsOnInvalidCutoff) {
  correlation::core::Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});
  updateTrajectory(cell);

  EXPECT_THROW(StructureAnalyzer(cell, 0.0, trajectory().getBondCutoffsSQ()), std::invalid_argument);
  EXPECT_THROW(StructureAnalyzer(cell, -1.0, trajectory().getBondCutoffsSQ()), std::invalid_argument);
}

TEST_F(StructureAnalyzerTests, HandlesEmptyCellGracefully) {
  correlation::core::Cell empty_cell;
  updateTrajectory(empty_cell);

  StructureAnalyzer const analyzer(empty_cell, 3.0, trajectory().getBondCutoffsSQ());
  EXPECT_TRUE(analyzer.distances().empty());
  EXPECT_TRUE(analyzer.angles().empty());
  EXPECT_TRUE(analyzer.dihedrals().empty());
  EXPECT_EQ(analyzer.neighborGraph().nodeCount(), 0);
}

} // namespace
} // namespace correlation::analysis

