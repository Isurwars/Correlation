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

} // namespace correlation::render3d
