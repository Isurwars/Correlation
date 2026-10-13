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
  const auto c_h = getElementColor(1);
  EXPECT_EQ(c_h.r, 245);
  EXPECT_EQ(c_h.g, 245);
  EXPECT_EQ(c_h.b, 245);

  const auto c_c = getElementColor("C");
  EXPECT_EQ(c_c.r, 64);
  EXPECT_EQ(c_c.g, 64);
  EXPECT_EQ(c_c.b, 64);

  const auto c_o = getElementColor("O");
  EXPECT_EQ(c_o.r, 255);
  EXPECT_EQ(c_o.g, 13);
  EXPECT_EQ(c_o.b, 13);

  const auto c_unknown = getElementColor("UnknownX");
  EXPECT_EQ(c_unknown.r, 255); // Fallback color (Deep Pink)

  const float r_c = getElementDisplayRadius("C");
  EXPECT_GT(r_c, 0.1f);
  EXPECT_LT(r_c, 2.0f);
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
  EXPECT_LE(scene->boundsMin().x(), 0.0f);
  EXPECT_GE(scene->boundsMax().x(), 10.0f);
}

} // namespace correlation::render3d
