// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"

#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::core;

class CellFunctionalTests : public ::testing::Test {};

TEST_F(CellFunctionalTests, BuildBCCLatticeAndVerifyDensity) {
  // Iron (Fe) has a BCC structure with a lattice parameter of ~2.866 Angstroms
  const double a = 2.866;
  Cell cell({{a, a, a, 90.0, 90.0, 90.0}});

  // BCC basis: (0,0,0) and (0.5, 0.5, 0.5) in fractional coordinates
  cell.addAtom("Fe", {0.0, 0.0, 0.0});
  cell.addAtom("Fe", {a * 0.5, a * 0.5, a * 0.5});

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_NEAR(cell.volume(), std::pow(a, 3), 1e-6);

  // Verify that wrapPositions doesn't move them if they are already inside
  cell.wrapPositions();
  EXPECT_NEAR(cell.atoms()[0].position().x(), 0.0, 1e-9);
  EXPECT_NEAR(cell.atoms()[1].position().x(), a * 0.5, 1e-9);
}

TEST_F(CellFunctionalTests, BuildFCCLatticeAndVerifyPBCDistances) {
  // Aluminum (Al) has an FCC structure with a lattice parameter of ~4.046 Angstroms
  const double a = 4.046;
  Cell cell({{a, a, a, 90.0, 90.0, 90.0}});

  // FCC basis
  cell.addAtom("Al", {0.0, 0.0, 0.0});
  cell.addAtom("Al", {0.0, a * 0.5, a * 0.5});
  cell.addAtom("Al", {a * 0.5, 0.0, a * 0.5});
  cell.addAtom("Al", {a * 0.5, a * 0.5, 0.0});

  EXPECT_EQ(cell.atomCount(), 4);

  // In FCC, the nearest neighbor distance is a/sqrt(2)
  const double expected_nn = a / std::sqrt(2.0);

  // Calculate distance between atom 0 and 1
  auto dist_vec = cell.atoms()[1].position() - cell.atoms()[0].position();
  auto mi_vec = cell.minimumImage(dist_vec);
  
  EXPECT_NEAR(correlation::math::norm(mi_vec), expected_nn, 1e-6);
}

TEST_F(CellFunctionalTests, VerifyWaterMoleculePBCStability) {
  // Simple water molecule in a large box
  Cell cell({{20.0, 20.0, 20.0, 90.0, 90.0, 90.0}});
  
  // O at origin, H atoms at typical distance/angle
  // OH distance ~ 0.96 A, HOH angle ~ 104.5
  const double oh_dist = 0.9584;
  const double hoh_angle_rad = 104.45 * (correlation::math::pi / 180.0);
  
  cell.addAtom("O", {10.0, 10.0, 10.0});
  cell.addAtom("H", {10.0 + oh_dist, 10.0, 10.0});
  cell.addAtom("H", {10.0 + oh_dist * std::cos(hoh_angle_rad), 
                     10.0 + oh_dist * std::sin(hoh_angle_rad), 
                     10.0});

  // Now move the molecule across the boundary
  for (auto &atom : const_cast<std::vector<Atom>&>(cell.atoms())) {
    auto pos = atom.position();
    pos.x() += 15.0; // Moves from 10 to 25, which should wrap to 5
    atom.setPosition(pos);
  }
  
  cell.wrapPositions();
  
  // Verify that the internal geometry (bond length/angle) is preserved
  const auto &atoms = cell.atoms();
  double d1 = distance(atoms[0], atoms[1]);
  double d2 = distance(atoms[0], atoms[2]);
  double ang = angle(atoms[0], atoms[1], atoms[2]);

  // Use minimum image for distance if they were wrapped differently
  auto v1 = cell.minimumImage(atoms[1].position() - atoms[0].position());
  auto v2 = cell.minimumImage(atoms[2].position() - atoms[0].position());
  
  EXPECT_NEAR(correlation::math::norm(v1), oh_dist, 1e-6);
  EXPECT_NEAR(correlation::math::norm(v2), oh_dist, 1e-6);
  
  // Angle function doesn't use PBC, so we must be careful. 
  // If we wrap them, they might be on opposite sides of the box.
  // We should calculate angle using minimum image vectors.
  double cos_theta = correlation::math::dot(v1, v2) / (correlation::math::norm(v1) * correlation::math::norm(v2));
  double calc_angle = std::acos(std::clamp(cos_theta, -1.0, 1.0));
  
  EXPECT_NEAR(calc_angle, hoh_angle_rad, 1e-6);
}

} // namespace correlation::testing
