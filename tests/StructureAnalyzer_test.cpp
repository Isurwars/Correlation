// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>

#include "../include/Atom.hpp"
#include "../include/Cell.hpp"
#include "../include/PhysicalData.hpp"
#include "../include/StructureAnalyzer.hpp"

// A test fixture for StructureAnalyzer tests.
class StructureAnalyzerTest : public ::testing::Test {};

TEST_F(StructureAnalyzerTest, FindsCorrectNeighborsForSilicon) {
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

  // Act: Calculate neighbors with a cutoff just beyond the first neighbor
  // shell. The first nearest neighbor distance in Si is (sqrt(3)/4)*a
  // approx 2.35 Å.
  StructureAnalyzer neighbors(si_cell, 3.0);
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

TEST_F(StructureAnalyzerTest, CalculatesCorrectAnglesForWater) {
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

  // Act: Calculate neighbors and angles.
  StructureAnalyzer neighbors(water_cell,
                              3.0); // Cutoff to include only H-O bonds
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

TEST_F(StructureAnalyzerTest, CalculatesCorrectAngleWithPBC) {
  // Arrange: Setup a system where the angle calculation requires the Minimum
  // Image Convention. Central atom B (index 1) at (0.5, 0.5, 0.5). Neighbor A
  // (index 0) at (3.5, 0.5, 0.5) -> PBC vector B->A is (-1.0, 0.0, 0.0).
  // Neighbor C (index 2) at (0.5, 3.5, 0.5) -> PBC vector B->C is (0.0, -1.0,
  // 0.0). The resulting angle A-B-C must be 90 degrees (pi/2 radians).

  const double side_length = 4.0;
  const double cutoff =
      1.1; // Cutoff > bond length (1.0) and < box half-side (2.0)
  // Calculate pi/2 explicitly using acos(-1.0) = pi.
  const double expected_angle_rad = std::acos(-1.0) / 2.0;

  Cell pbc_cell({side_length, side_length, side_length, 90.0, 90.0, 90.0});

  pbc_cell.addAtom("C", {3.5, 0.5, 0.5}); // Atom A
  pbc_cell.addAtom("C", {0.5, 0.5, 0.5}); // Atom B (Central)
  pbc_cell.addAtom("O", {0.5, 3.5, 0.5}); // Atom C

  // Act: Calculate neighbors and angles.
  StructureAnalyzer analyzer(pbc_cell, cutoff);
  const auto &angles = analyzer.angles();

  // Assert
  // Retrieve the element ID for indexing (should be 0).
  const int c_id = pbc_cell.findElement("C")->id.value;
  ASSERT_EQ(c_id, 0);

  // Retrieve the element ID for indexing (should be 0).
  const int o_id = pbc_cell.findElement("O")->id.value;
  ASSERT_EQ(o_id, 1);

  // The angle is stored in the [C][C][O] slot.
  const auto &cco_angles = angles[c_id][c_id][o_id];

  // There must be exactly one angle calculated (A-B-C).
  ASSERT_EQ(cco_angles.size(), 1);

  // The calculated angle should be pi/2 radians.
  EXPECT_NEAR(cco_angles[0], expected_angle_rad, 1e-6);
}
