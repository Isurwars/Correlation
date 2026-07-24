// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/LinearAlgebra.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {

using namespace correlation::core;

namespace {
class TrajectoryFunctionalTests : public ::testing::Test {
protected:
  static Cell createSingleAtomCell(const std::string &element, const math::Vector3R &pos) {
    Cell cell({{static_cast<real_t>(10.0), static_cast<real_t>(10.0), static_cast<real_t>(10.0),
                static_cast<real_t>(90.0), static_cast<real_t>(90.0), static_cast<real_t>(90.0)}});
    cell.addAtom(element, pos);
    return cell;
  }
};

TEST_F(TrajectoryFunctionalTests, VerifyPBCDiffusionVelocityCalculation) {
  // Simulate an atom diffusing across the periodic boundary in the X-direction
  // Box size is 10.0.
  // Frame 0: pos = 9.5
  // Frame 1: pos = 1.5 (crossing the right boundary)
  // Frame 2: pos = 3.0
  const real_t time_step = static_cast<real_t>(0.5); // fs
  std::vector<Cell> frames;
  frames.push_back(createSingleAtomCell("H", {9.5, 5.0, 5.0}));
  frames.push_back(createSingleAtomCell("H", {1.5, 5.0, 5.0}));
  frames.push_back(createSingleAtomCell("H", {3.0, 5.0, 5.0}));

  Trajectory traj(frames, time_step);
  EXPECT_EQ(traj.getFrameCount(), 3);
  EXPECT_DOUBLE_EQ(traj.getTimeStep(), time_step);

  traj.calculateVelocities();

  // For Frame 0:
  // Forward difference: displacement from 9.5 to 1.5.
  // Delta_r_raw = { -8.0, 0.0, 0.0 }
  // With PBC: minimumImage(-8.0) = 2.0 (the atom moved forward by 2.0 A across the boundary)
  // Velocity = displacement / dt = 2.0 / 0.5 = 4.0 A/fs
  EXPECT_NEAR(traj.getFrame(0).atoms()[0].velocity().x(), 4.0, 1e-6);
  EXPECT_NEAR(traj.getFrame(0).atoms()[0].velocity().y(), 0.0, 1e-6);
  EXPECT_NEAR(traj.getFrame(0).atoms()[0].velocity().z(), 0.0, 1e-6);

  // For Frame 1 (Central difference):
  // Displacement between Frame 2 (3.0) and Frame 0 (9.5).
  // Delta_r_raw = { -6.5, 0.0, 0.0 }
  // With PBC: minimumImage(-6.5) = 3.5
  // Velocity = displacement / (2 * dt) = 3.5 / 1.0 = 3.5 A/fs
  EXPECT_NEAR(traj.getFrame(1).atoms()[0].velocity().x(), 3.5, 1e-6);

  // For Frame 2:
  // Backward difference: displacement between Frame 2 (3.0) and Frame 1 (1.5).
  // Delta_r_raw = { 1.5, 0.0, 0.0 }
  // With PBC: minimumImage(1.5) = 1.5
  // Velocity = displacement / dt = 1.5 / 0.5 = 3.0 A/fs
  EXPECT_NEAR(traj.getFrame(2).atoms()[0].velocity().x(), 3.0, 1e-6);
}

TEST_F(TrajectoryFunctionalTests, VerifyValidationOfMismatchedFrames) {
  Trajectory traj;
  traj.addFrame(createSingleAtomCell("C", {0.0, 0.0, 0.0}));

  // 1. Mismatched atom count
  Cell cell_two_atoms({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  cell_two_atoms.addAtom("C", {1.0, 1.0, 1.0});
  cell_two_atoms.addAtom("C", {2.0, 2.0, 2.0});
  EXPECT_THROW(traj.addFrame(cell_two_atoms), std::runtime_error);

  // 2. Mismatched element type
  Cell cell_wrong_type = createSingleAtomCell("O", {0.0, 0.0, 0.0});
  EXPECT_THROW(traj.addFrame(cell_wrong_type), std::runtime_error);

  // 3. Mismatched element order (for multi-atom system)
  Trajectory multi_traj;
  Cell ref_cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  ref_cell.addAtom("C", {0.0, 0.0, 0.0});
  ref_cell.addAtom("H", {1.0, 1.0, 1.0});
  multi_traj.addFrame(ref_cell);

  Cell bad_order_cell({{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}});
  bad_order_cell.addAtom("H", {1.0, 1.0, 1.0});
  bad_order_cell.addAtom("C", {0.0, 0.0, 0.0});
  EXPECT_THROW(multi_traj.addFrame(bad_order_cell), std::runtime_error);
}

TEST_F(TrajectoryFunctionalTests, VerifyDeduplicationStatTracking) {
  // Test frame deduplication stats on a trajectory with repeated coordinates
  std::vector<Cell> frames = {
      createSingleAtomCell("Si", {1.0, 1.0, 1.0}), createSingleAtomCell("Si", {1.0, 1.0, 1.0}), // duplicate 1
      createSingleAtomCell("Si", {2.0, 2.0, 2.0}), createSingleAtomCell("Si", {2.0, 2.0, 2.0}), // duplicate 2
      createSingleAtomCell("Si", {2.0, 2.0, 2.0}),                                              // duplicate 3
      createSingleAtomCell("Si", {3.0, 3.0, 3.0})};

  Trajectory traj(frames, 1.0);
  // Deduplication happens in the constructor for vector-loaded trajectories
  EXPECT_EQ(traj.getFrameCount(), 3);
  EXPECT_EQ(traj.getRemovedFrameCount(), 3);

  EXPECT_NEAR(traj.getFrame(0).atoms()[0].position().x(), 1.0, 1e-9);
  EXPECT_NEAR(traj.getFrame(1).atoms()[0].position().x(), 2.0, 1e-9);
  EXPECT_NEAR(traj.getFrame(2).atoms()[0].position().x(), 3.0, 1e-9);
}

} // namespace
} // namespace correlation::testing
