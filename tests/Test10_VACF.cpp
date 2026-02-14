// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <gtest/gtest.h>
#include <vector>

#include "../include/Cell.hpp"
#include "../include/DistributionFunctions.hpp"
#include "../include/Trajectory.hpp"

// Test fixture for VACF and VDOS tests.
class Test10_VACF : public ::testing::Test {
protected:
  // No special setup needed for VACF usually, or different from Test06
};

TEST_F(Test10_VACF, CalculateVACF_and_VDOS) {
  Cell c({10, 10, 10, 90, 90, 90});
  c.addAtom("Ar", {0, 0, 0});
  Trajectory t;
  t.addFrame(c);
  t.addFrame(c); // Static
  t.calculateVelocities();

  DistributionFunctions df(c, 0.0, {{0.0}});

  df.calculateVACF(t, 1);
  EXPECT_NO_THROW(df.getHistogram("VACF"));
  const auto &vacf = df.getHistogram("VACF").partials.at("Total");
  // Should be 1.0 normalized (if constant 0 velocity? Wait 0 velocity -> ??)
  // If static, position constant -> velocity 0. Correlation of 0 with 0 is 0.
  // Let's give it velocity.

  Trajectory tMoving;
  Cell c1 = c;

  Cell c2({10, 10, 10, 90, 90, 90});
  // Atom 1 moves +1.0 in x
  c2.addAtom("Ar", {1.0, 0.0, 0.0});
  // Atom 2 moves -1.0 in x (balancing COM)
  c2.addAtom("Ar", {-1.0, 0.0, 0.0});

  Cell c3({10, 10, 10, 90, 90, 90});
  // Atom 1 moves to +2.0
  c3.addAtom("Ar", {2.0, 0.0, 0.0});
  // Atom 2 moves to -2.0
  c3.addAtom("Ar", {-2.0, 0.0, 0.0});

  // Need to update c1 (frame 0) to have 2 atoms at 0
  c1 = Cell({10, 10, 10, 90, 90, 90});
  c1.addAtom("Ar", {0.0, 0.0, 0.0});
  c1.addAtom("Ar", {0.0, 0.0, 0.0});

  tMoving.addFrame(c1);
  tMoving.addFrame(c2);
  tMoving.addFrame(c3);
  tMoving.setTimeStep(1.0);
  tMoving.calculateVelocities();

  df.calculateVACF(tMoving, 1);
  const auto &vacf2 = df.getHistogram("Normalized VACF").partials.at("Total");
  EXPECT_NEAR(vacf2[0], 1.0, 1e-6); // t=0
  EXPECT_NEAR(vacf2[1], 1.0, 1e-6); // t=1, const velocity

  // VDOS
  df.calculateVDOS();
  EXPECT_NO_THROW(df.getHistogram("VDOS"));
  const auto &vdos_hist = df.getHistogram("VDOS");
  EXPECT_TRUE(vdos_hist.partials.count("Frequency (cm-1)"));
}
