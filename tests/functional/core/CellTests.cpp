// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"

#include <gtest/gtest.h>
#include <numbers>

namespace correlation::testing {

using namespace correlation::core;

namespace {
class CellFunctionalTests : public ::testing::Test {};

TEST_F(CellFunctionalTests, BuildBCCLatticeAndVerifyDensity) {
  // Iron (Fe) has a BCC structure with a lattice parameter of ~2.866 Angstroms
  const real_t lattice_parameter = 2.866;
  Cell cell({{lattice_parameter, lattice_parameter, lattice_parameter, 90.0, 90.0, 90.0}});

  // BCC basis: (0,0,0) and (0.5, 0.5, 0.5) in fractional coordinates
  cell.addAtom("Fe", {0.0, 0.0, 0.0});
  cell.addAtom("Fe", {lattice_parameter * 0.5, lattice_parameter * 0.5, lattice_parameter * 0.5});

  EXPECT_EQ(cell.atomCount(), 2);
  EXPECT_NEAR(cell.volume(), std::pow(lattice_parameter, 3), correlation::is_single_precision ? 1e-5 : 1e-6);

  // Verify that wrapPositions doesn't move them if they are already inside
  cell.wrapPositions();
  EXPECT_NEAR(cell.atoms()[0].position().x(), 0.0, correlation::is_single_precision ? 1e-6 : 1e-9);
  EXPECT_NEAR(cell.atoms()[1].position().x(), lattice_parameter * 0.5, correlation::is_single_precision ? 1e-6 : 1e-9);
}

TEST_F(CellFunctionalTests, BuildFCCLatticeAndVerifyPBCDistances) {
  // Aluminum (Al) has an FCC structure with a lattice parameter of ~4.046 Angstroms
  const real_t lattice_parameter = 4.046;
  Cell cell({{lattice_parameter, lattice_parameter, lattice_parameter, 90.0, 90.0, 90.0}});

  // FCC basis
  cell.addAtom("Al", {0.0, 0.0, 0.0});
  cell.addAtom("Al", {0.0, lattice_parameter * 0.5, lattice_parameter * 0.5});
  cell.addAtom("Al", {lattice_parameter * 0.5, 0.0, lattice_parameter * 0.5});
  cell.addAtom("Al", {lattice_parameter * 0.5, lattice_parameter * 0.5, 0.0});

  EXPECT_EQ(cell.atomCount(), 4);

  // In FCC, the nearest neighbor distance is a/sqrt(2)
  const auto expected_nn = static_cast<real_t>(lattice_parameter / std::numbers::sqrt2);

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
  const real_t oh_dist = 0.9584;
  const real_t hoh_angle_rad = 104.45 * (correlation::math::pi / 180.0);

  cell.addAtom("O", {10.0, 10.0, 10.0});
  cell.addAtom("H", {10.0 + oh_dist, 10.0, 10.0});
  cell.addAtom("H", {10.0 + oh_dist * std::cos(hoh_angle_rad), 10.0 + oh_dist * std::sin(hoh_angle_rad), 10.0});

  // Now move the molecule across the boundary
  for (auto &atom : const_cast<std::vector<Atom> &>(cell.atoms())) {
    auto pos = atom.position();
    pos.x() += 15.0; // Moves from 10 to 25, which should wrap to 5
    atom.setPosition(pos);
  }

  cell.wrapPositions();

  // Verify that the internal geometry (bond length/angle) is preserved
  const auto &atoms = cell.atoms();
  real_t const distance_1_2 = distance(atoms[0], atoms[1]);
  real_t const distance_1_3 = distance(atoms[0], atoms[2]);
  real_t const angle_1_2_3 = angle(atoms[0], atoms[1], atoms[2]);

  // Use minimum image for distance if they were wrapped differently
  auto image_1 = cell.minimumImage(atoms[1].position() - atoms[0].position());
  auto image_2 = cell.minimumImage(atoms[2].position() - atoms[0].position());

  EXPECT_NEAR(correlation::math::norm(image_1), oh_dist, correlation::is_single_precision ? 1e-5 : 1e-6);
  EXPECT_NEAR(correlation::math::norm(image_2), oh_dist, correlation::is_single_precision ? 1e-5 : 1e-6);

  // Angle function doesn't use PBC, so we must be careful.
  // If we wrap them, they might be on opposite sides of the box.
  // We should calculate angle using minimum image vectors.
  real_t const cos_theta =
      correlation::math::dot(image_1, image_2) / (correlation::math::norm(image_1) * correlation::math::norm(image_2));
  real_t const calc_angle = std::acos(std::clamp(cos_theta, static_cast<real_t>(-1.0), static_cast<real_t>(1.0)));

  EXPECT_NEAR(calc_angle, hoh_angle_rad, correlation::is_single_precision ? 1e-4 : 1e-6);
}

TEST_F(CellFunctionalTests, VerifyLatticeParameterGuards) {
  // 1. Negative/zero lengths
  EXPECT_THROW(Cell({{-1.0, 10.0, 10.0, 90.0, 90.0, 90.0}}), std::invalid_argument);
  EXPECT_THROW(Cell({{10.0, 0.0, 10.0, 90.0, 90.0, 90.0}}), std::invalid_argument);

  // 2. Invalid angles
  EXPECT_THROW(Cell({{10.0, 10.0, 10.0, 0.0, 90.0, 90.0}}), std::invalid_argument);
  EXPECT_THROW(Cell({{10.0, 10.0, 10.0, 90.0, 180.0, 90.0}}), std::invalid_argument);
  EXPECT_THROW(Cell({{10.0, 10.0, 10.0, 90.0, 90.0, -10.0}}), std::invalid_argument);

  // 3. NaN values
  const real_t nan_val = std::numeric_limits<real_t>::quiet_NaN();
  EXPECT_THROW(Cell({{nan_val, 10.0, 10.0, 90.0, 90.0, 90.0}}), std::invalid_argument);
  EXPECT_THROW(Cell({{10.0, 10.0, 10.0, nan_val, 90.0, 90.0}}), std::invalid_argument);
}
} // namespace
} // namespace correlation::testing
