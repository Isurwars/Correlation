/**
 * @file Structure3DControllerTests.cpp
 * @brief Unit tests for Structure3DController.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "AppWindow.h"
#include "app/core/AppOptions.hpp"
#include "app/viewmodel/Structure3DController.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <chrono>
#include <gtest/gtest.h>
#include <optional>
#include <stdexcept>
#include <thread>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

namespace correlation::app {
namespace {

core::Cell createTestCell(float offset = 0.0F) {
  core::Cell cell(correlation::math::Vector3<real_t>{15.0F, 0.0F, 0.0F},
                  correlation::math::Vector3<real_t>{0.0F, 15.0F, 0.0F},
                  correlation::math::Vector3<real_t>{0.0F, 0.0F, 15.0F});
  cell.addAtom("Si", correlation::math::Vector3<real_t>{offset, 0.0F, 0.0F});
  cell.addAtom("O", correlation::math::Vector3<real_t>{1.6F + offset, 0.0F, 0.0F});
  cell.addAtom("O", correlation::math::Vector3<real_t>{0.0F, 1.6F + offset, 0.0F});
  return cell;
}

class Structure3DControllerTests : public ::testing::Test {
public:
  [[nodiscard]] AppWindow &window() {
    if (!window_.has_value()) {
      throw std::runtime_error("Window is not initialized");
    }
    return **window_;
  }

protected:
  void SetUp() override {
#ifndef _WIN32
    setenv("SLINT_BACKEND", "software", 1);
#else
    _putenv_s("SLINT_BACKEND", "software");
#endif
    window_.emplace(AppWindow::create());
  }

private:
  std::optional<slint::ComponentHandle<AppWindow>> window_;
};

TEST_F(Structure3DControllerTests, InitializationAndDefaults) {
  const ProgramOptions options;
  const Structure3DController controller(window(), options);

  EXPECT_FALSE(controller.is3DModeActive());
  EXPECT_EQ(controller.scene(), nullptr);
  EXPECT_GT(controller.camera().fov(), 0.0F);
}

TEST_F(Structure3DControllerTests, ViewportResize) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  controller.handleViewportResized(1024.0F, 768.0F);
  EXPECT_NEAR(controller.camera().aspectRatio(), 1024.0F / 768.0F, 0.01F);

  // Resize below minimum threshold should be ignored
  const float aspect_before = controller.camera().aspectRatio();
  controller.handleViewportResized(5.0F, 5.0F);
  EXPECT_FLOAT_EQ(controller.camera().aspectRatio(), aspect_before);
}

TEST_F(Structure3DControllerTests, CameraNavigation) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  const float initial_azimuth = controller.camera().azimuth();
  controller.handleRotate(15.0F, 10.0F);
  EXPECT_NE(controller.camera().azimuth(), initial_azimuth);

  const auto initial_target = controller.camera().target();
  controller.handlePan(5.0F, -5.0F);
  EXPECT_NE(controller.camera().target().x(), initial_target.x());

  const float initial_distance = controller.camera().distance();
  controller.handleZoom(1.5F);
  EXPECT_LT(controller.camera().distance(), initial_distance);
}

TEST_F(Structure3DControllerTests, SetCellAndResetView) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  const auto cell = createTestCell();
  controller.setCell(cell);

  ASSERT_NE(controller.scene(), nullptr);
  EXPECT_EQ(controller.scene()->atoms().size(), 3UL);
  EXPECT_EQ(window().get_structure_3d_current_frame(), 0);
  EXPECT_EQ(window().get_structure_3d_total_frames(), 1);

  controller.handleResetView();
  EXPECT_GT(controller.camera().distance(), 0.0F);
}

TEST_F(Structure3DControllerTests, FeatureToggles) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  controller.handleToggleBonds(false);
  EXPECT_FALSE(window().get_structure_3d_bonds_enabled());

  controller.handleToggleBonds(true);
  EXPECT_TRUE(window().get_structure_3d_bonds_enabled());

  controller.handleToggleWireframe(false);
  EXPECT_FALSE(window().get_structure_3d_wireframe_enabled());

  controller.handleToggleWireframe(true);
  EXPECT_TRUE(window().get_structure_3d_wireframe_enabled());
}

TEST_F(Structure3DControllerTests, TrajectoryNavigation) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  const auto cell1 = createTestCell(0.0F);
  const auto cell2 = createTestCell(0.5F);
  const auto traj = std::make_shared<core::Trajectory>(std::vector<core::Cell>{cell1, cell2}, 1.0F);

  controller.setTrajectory(traj);
  EXPECT_EQ(window().get_structure_3d_total_frames(), 2);
  EXPECT_EQ(window().get_structure_3d_current_frame(), 0);
  ASSERT_NE(controller.scene(), nullptr);
  ASSERT_FALSE(controller.scene()->bonds().empty());
  EXPECT_NEAR(controller.scene()->bonds()[0].start.x(), 0.0F, 1e-4F);

  // Stepping forward updates both frame index and bond positions to frame 1
  controller.handleStepFrame(true);
  EXPECT_EQ(window().get_structure_3d_current_frame(), 1);
  ASSERT_FALSE(controller.scene()->bonds().empty());
  EXPECT_NEAR(controller.scene()->bonds()[0].start.x(), 0.5F, 1e-4F);

  // Stepping backward restores frame 0 and updates bonds back to frame 0
  controller.handleStepFrame(false);
  EXPECT_EQ(window().get_structure_3d_current_frame(), 0);
  ASSERT_FALSE(controller.scene()->bonds().empty());
  EXPECT_NEAR(controller.scene()->bonds()[0].start.x(), 0.0F, 1e-4F);

  // Seeking directly to frame 1 updates bond positions
  controller.handleSeekFrame(1);
  EXPECT_EQ(window().get_structure_3d_current_frame(), 1);
  ASSERT_FALSE(controller.scene()->bonds().empty());
  EXPECT_NEAR(controller.scene()->bonds()[0].start.x(), 0.5F, 1e-4F);

  // Scaling cutoff updates bonds in scene
  controller.updateBondCutoffScale(1.5);
  ASSERT_FALSE(controller.scene()->bonds().empty());
}

TEST_F(Structure3DControllerTests, ActivationAndRendering) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  const auto cell = createTestCell();
  controller.setCell(cell);

  controller.set3DModeActive(true);
  EXPECT_TRUE(controller.is3DModeActive());

  controller.requestRender(true);

  // Allow background worker thread to process at least one frame
  std::this_thread::sleep_for(std::chrono::milliseconds(50));

  controller.set3DModeActive(false);
  EXPECT_FALSE(controller.is3DModeActive());
}

} // namespace
} // namespace correlation::app
