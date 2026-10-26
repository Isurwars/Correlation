// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "core/VelocityIntegrator.hpp"

#include <gtest/gtest.h>

namespace correlation::testing {

using namespace correlation::core;

TEST(VelocityIntegratorTests, IgnoresSingleFrameOrZeroTimeStep) {
  Cell cell({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  std::vector<Cell> frames{cell};

  // 1 frame
  VelocityIntegrator::calculateVelocities(frames, 1.0);
  EXPECT_NEAR(math::norm(frames[0].atoms()[0].velocity()), 0.0, 1e-6);

  // 2 frames with dt = 0.0
  Cell cell2 = cell;
  cell2.atoms()[0].setPosition({1.0, 0.0, 0.0});
  frames.push_back(cell2);

  VelocityIntegrator::calculateVelocities(frames, 0.0);
  EXPECT_NEAR(math::norm(frames[0].atoms()[0].velocity()), 0.0, 1e-6);
}

TEST(VelocityIntegratorTests, ComputesFiniteDifferencesCorrectly) {
  const real_t time_step = 2.0;

  Cell frame0({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  frame0.addAtom("Si", {0.0, 0.0, 0.0});

  Cell frame1({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  frame1.addAtom("Si", {2.0, 0.0, 0.0});

  Cell frame2({20.0, 20.0, 20.0, 90.0, 90.0, 90.0});
  frame2.addAtom("Si", {6.0, 0.0, 0.0});

  std::vector<Cell> frames{frame0, frame1, frame2};
  VelocityIntegrator::calculateVelocities(frames, time_step);

  // Frame 0: forward difference (frame1 - frame0) / time_step = (2.0 - 0.0) / 2.0 = 1.0
  EXPECT_NEAR(frames[0].atoms()[0].velocity().x(), 1.0, 1e-5);

  // Frame 1: central difference (frame2 - frame0) / (2 * time_step) = (6.0 - 0.0) / 4.0 = 1.5
  EXPECT_NEAR(frames[1].atoms()[0].velocity().x(), 1.5, 1e-5);

  // Frame 2: backward difference (frame2 - frame1) / time_step = (6.0 - 2.0) / 2.0 = 2.0
  EXPECT_NEAR(frames[2].atoms()[0].velocity().x(), 2.0, 1e-5);
}

TEST(VelocityIntegratorTests, IntegratorOperatesOnTrajectoryInstance) {
  Cell frame0({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  frame0.addAtom("Si", {0.0, 0.0, 0.0});

  Cell frame1({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
  frame1.addAtom("Si", {1.0, 0.0, 0.0});

  Trajectory traj({frame0, frame1}, 1.0);
  VelocityIntegrator::calculateVelocities(traj);

  EXPECT_NEAR(traj.getFrame(0).atoms()[0].velocity().x(), 1.0, 1e-5);
  EXPECT_NEAR(traj.getFrame(1).atoms()[0].velocity().x(), 1.0, 1e-5);
}

} // namespace correlation::testing
