/**
 * @file Scene3DTests.cpp
 * @brief Unit tests for Scene3D geometry extraction, unit cell wireframe, and CPK coloring.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "core/Cell.hpp"
#include "render3d/ElementColors.hpp"
#include "render3d/Scene3D.hpp"
#include <gtest/gtest.h>

namespace correlation::render3d {

TEST(Scene3DTests, ElementColoringPalette) {
  // Hydrogen (Z = 1) - ASE/Jmol white
  const auto c_h = getElementColor(1);
  EXPECT_EQ(c_h.r, 255);
  EXPECT_EQ(c_h.g, 255);
  EXPECT_EQ(c_h.b, 255);

  // Carbon (Z = 6) - ASE/Jmol gray (#909090)
  const auto c_c = getElementColor("C");
  EXPECT_EQ(c_c.r, 144);
  EXPECT_EQ(c_c.g, 144);
  EXPECT_EQ(c_c.b, 144);

  // Oxygen (Z = 8) - ASE/Jmol red
  const auto c_o = getElementColor("O");
  EXPECT_EQ(c_o.r, 255);
  EXPECT_EQ(c_o.g, 13);
  EXPECT_EQ(c_o.b, 13);

  // Gold (Z = 79)
  const auto c_au = getElementColor("Au");
  EXPECT_EQ(c_au.r, 255);
  EXPECT_EQ(c_au.g, 209);
  EXPECT_EQ(c_au.b, 35);

  // Superheavy: Oganesson (Z = 118)
  const auto c_og = getElementColor("Og");
  EXPECT_EQ(c_og.r, 252);
  EXPECT_EQ(c_og.g, 0);
  EXPECT_EQ(c_og.b, 0);
  EXPECT_EQ(c_og, getElementColor(118));

  // Symbol to atomic number mapping
  EXPECT_EQ(getAtomicNumber("H"), 1);
  EXPECT_EQ(getAtomicNumber("C"), 6);
  EXPECT_EQ(getAtomicNumber("Og"), 118);
  EXPECT_EQ(getAtomicNumber("UnknownX"), 0);

  // Out-of-bounds and unknown fallback (Deep Pink #FF1493)
  const auto c_unknown = getElementColor("UnknownX");
  EXPECT_EQ(c_unknown.r, 255);
  EXPECT_EQ(c_unknown.g, 20);
  EXPECT_EQ(c_unknown.b, 147);

  const auto c_zero = getElementColor(0);
  EXPECT_EQ(c_zero, c_unknown);

  const auto c_oob = getElementColor(119);
  EXPECT_EQ(c_oob, c_unknown);

  const float r_c = getElementDisplayRadius("C");
  EXPECT_GT(r_c, 0.1F);
  EXPECT_LT(r_c, 2.0F);

  // Numeric LAMMPS type strings ("1".."4") receive distinct categorical colors
  const auto c_t1 = getElementColor("1");
  const auto c_t2 = getElementColor("2");
  const auto c_t3 = getElementColor("3");
  const auto c_t4 = getElementColor("4");
  EXPECT_NE(c_t1, c_unknown);
  EXPECT_NE(c_t1, c_t2);
  EXPECT_NE(c_t2, c_t3);
  EXPECT_NE(c_t3, c_t4);
}

TEST(Scene3DTests, EmptySceneBounds) {
  Scene3D scene;
  scene.updateBounds();
  EXPECT_LE(scene.boundsMin().x(), scene.boundsMax().x());
  EXPECT_LE(scene.boundsMin().y(), scene.boundsMax().y());
  EXPECT_LE(scene.boundsMin().z(), scene.boundsMax().z());
}

TEST(Scene3DTests, BuildSceneFromCell) {
  core::Cell cell;
  const math::Matrix3<real_t> lattice(math::Vector3<real_t>(10.0, 0.0, 0.0),
                                      math::Vector3<real_t>(0.0, 10.0, 0.0),
                                      math::Vector3<real_t>(0.0, 0.0, 10.0));
  cell.updateLattice(lattice);

  cell.addAtom("C", math::Vector3<real_t>(5.0, 5.0, 5.0));
  cell.addAtom("H", math::Vector3<real_t>(5.0, 5.0, 6.09)); // 1.09 A bond distance

  const auto scene = Scene3D::fromCell(cell, 1.2, true);
  ASSERT_NE(scene, nullptr);

  // Validate atoms
  EXPECT_EQ(scene->atoms().size(), 2);
  EXPECT_EQ(scene->atoms()[0].color, getElementColor("C"));
  EXPECT_EQ(scene->atoms()[1].color, getElementColor("H"));

  // Validate unit cell wireframe (12 edges)
  EXPECT_EQ(scene->unitCellLines().size(), 12);

  // Validate C-H bond was generated
  EXPECT_EQ(scene->bonds().size(), 1);

  // Validate bounds enclose both unit cell and atoms
  EXPECT_LE(scene->boundsMin().x(), 0.0F);
  EXPECT_GE(scene->boundsMax().x(), 10.0F);
}

TEST(Scene3DTests, BuildSceneInCellFoldAtoms) {
  core::Cell cell;
  const math::Matrix3<real_t> lattice(math::Vector3<real_t>(10.0, 0.0, 0.0),
                                      math::Vector3<real_t>(0.0, 10.0, 0.0),
                                      math::Vector3<real_t>(0.0, 0.0, 10.0));
  cell.updateLattice(lattice);

  // Add atom outside unit cell
  cell.addAtom("C", math::Vector3<real_t>(12.5, -2.5, 23.5));

  // In Standard mode (in_cell = false), coordinates remain unwrapped
  const auto standard_scene = Scene3D::fromCell(cell, 1.2, false, false);
  ASSERT_NE(standard_scene, nullptr);
  EXPECT_FLOAT_EQ(standard_scene->atoms()[0].position.x(), 12.5F);
  EXPECT_FLOAT_EQ(standard_scene->atoms()[0].position.y(), -2.5F);
  EXPECT_FLOAT_EQ(standard_scene->atoms()[0].position.z(), 23.5F);

  // In In-Cell mode (in_cell = true), coordinates are folded into [0, 10)
  const auto in_cell_scene = Scene3D::fromCell(cell, 1.2, false, true);
  ASSERT_NE(in_cell_scene, nullptr);
  EXPECT_NEAR(in_cell_scene->atoms()[0].position.x(), 2.5F, 1e-4F);
  EXPECT_NEAR(in_cell_scene->atoms()[0].position.y(), 7.5F, 1e-4F);
  EXPECT_NEAR(in_cell_scene->atoms()[0].position.z(), 3.5F, 1e-4F);
}

TEST(Scene3DTests, BuildSceneInCellPeriodicBondsAndClipping) {
  core::Cell cell;
  const math::Matrix3<real_t> lattice(math::Vector3<real_t>(10.0, 0.0, 0.0),
                                      math::Vector3<real_t>(0.0, 10.0, 0.0),
                                      math::Vector3<real_t>(0.0, 0.0, 10.0));
  cell.updateLattice(lattice);

  // Two C atoms near opposite x faces: x=0.5 and x=9.5.
  // Direct Euclidean distance is 9.0 A (no bond without PBC).
  // Under PBC, distance is 1.0 A (well within covalent bond cutoff).
  cell.addAtom("C", math::Vector3<real_t>(0.5, 5.0, 5.0));
  cell.addAtom("C", math::Vector3<real_t>(9.5, 5.0, 5.0));

  // Standard mode: no bond
  const auto standard_scene = Scene3D::fromCell(cell, 1.2, true, false);
  ASSERT_NE(standard_scene, nullptr);
  EXPECT_EQ(standard_scene->bonds().size(), 0);

  // In-Cell mode: bond crosses x=0 and x=10 boundary faces, resulting in 2 boundary stubs
  const auto in_cell_scene = Scene3D::fromCell(cell, 1.2, true, true);
  ASSERT_NE(in_cell_scene, nullptr);
  EXPECT_EQ(in_cell_scene->bonds().size(), 2);

  // Verify stub 1 goes from x=0.5 to x=0.0
  const auto &bond_stub1 = in_cell_scene->bonds()[0];
  EXPECT_NEAR(bond_stub1.start.x(), 0.5F, 1e-4F);
  EXPECT_NEAR(bond_stub1.end.x(), 0.0F, 1e-4F);

  // Verify stub 2 goes from x=9.5 to x=10.0
  const auto &bond_stub2 = in_cell_scene->bonds()[1];
  EXPECT_NEAR(bond_stub2.start.x(), 9.5F, 1e-4F);
  EXPECT_NEAR(bond_stub2.end.x(), 10.0F, 1e-4F);
}

TEST(Scene3DTests, SetAtomScaleUpdatesRadiiAndClamps) {
  core::Cell cell;
  cell.addAtom("C", math::Vector3<real_t>(1.0, 2.0, 3.0));

  const auto scene = Scene3D::fromCell(cell, 1.15, false, false);
  ASSERT_NE(scene, nullptr);
  ASSERT_EQ(scene->atoms().size(), 1UL);
  EXPECT_FLOAT_EQ(scene->atomScale(), 0.35F);

  const float cov_r = scene->atoms()[0].covalent_radius;
  EXPECT_NEAR(scene->atoms()[0].radius, cov_r * 0.35F, 1e-5F);

  scene->setAtomScale(0.85F);
  EXPECT_FLOAT_EQ(scene->atomScale(), 0.85F);
  EXPECT_NEAR(scene->atoms()[0].radius, cov_r * 0.85F, 1e-5F);

  // Clamp upper bound to 2.0F and lower bound to 0.05F
  scene->setAtomScale(5.0F);
  EXPECT_FLOAT_EQ(scene->atomScale(), 2.0F);
  scene->setAtomScale(0.01F);
  EXPECT_FLOAT_EQ(scene->atomScale(), 0.05F);
}

} // namespace correlation::render3d
