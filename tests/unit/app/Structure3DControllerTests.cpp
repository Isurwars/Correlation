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
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <optional>
#include <stdexcept>
#include <thread>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>

using correlation::real_t;
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

  // In-cell mode defaults to true
  EXPECT_TRUE(controller.isInCell());
  EXPECT_TRUE(window().get_structure_3d_in_cell());

  controller.handleToggleInCell(false);
  EXPECT_FALSE(controller.isInCell());
  EXPECT_FALSE(window().get_structure_3d_in_cell());

  controller.handleToggleInCell(true);
  EXPECT_TRUE(controller.isInCell());
  EXPECT_TRUE(window().get_structure_3d_in_cell());

  // Isometric mode defaults to false (Perspective)
  EXPECT_FALSE(controller.isIsometric());
  EXPECT_FALSE(window().get_structure_3d_isometric());

  controller.handleToggleProjection(true);
  EXPECT_TRUE(controller.isIsometric());
  EXPECT_TRUE(controller.camera().isOrthographic());
  EXPECT_TRUE(window().get_structure_3d_isometric());

  controller.handleToggleProjection(false);
  EXPECT_FALSE(controller.isIsometric());
  EXPECT_FALSE(controller.camera().isOrthographic());
  EXPECT_FALSE(window().get_structure_3d_isometric());
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

TEST_F(Structure3DControllerTests, ConcurrentResizeAndTrajectoryRenderDoesNotCrash) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  const auto cell1 = createTestCell(0.0F);
  const auto cell2 = createTestCell(0.5F);
  const auto traj = std::make_shared<core::Trajectory>(std::vector<core::Cell>{cell1, cell2}, 1.0F);

  controller.setTrajectory(traj);
  controller.set3DModeActive(true);

  for (int i = 0; i < 20; ++i) {
    const float width = 640.0F + static_cast<float>((i % 5) * 32);
    const float height = 480.0F + static_cast<float>((i % 5) * 24);
    controller.handleViewportResized(width, height);
    controller.handleSeekFrame(i % 2);
    controller.requestRender(true);
    std::this_thread::sleep_for(std::chrono::milliseconds(2));
  }

  controller.set3DModeActive(false);
  EXPECT_FALSE(controller.is3DModeActive());
}

TEST_F(Structure3DControllerTests, SlabClippingAndResetCullsAndRestoresAtoms) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);
  controller.setCell(createTestCell());

  ASSERT_NE(controller.scene(), nullptr);
  ASSERT_GE(controller.scene()->atoms().size(), 2U);

  controller.setSlabClipConfig(SlabClipConfig{
      .enabled = true,
      .x_min = 0.0F,
      .x_max = 0.05F,
      .y_min = 0.0F,
      .y_max = 1.0F,
      .z_min = 0.0F,
      .z_max = 1.0F,
  });

  bool any_culled = false;
  for (const auto &atom : controller.scene()->atoms()) {
    if (atom.dim_factor <= 0.0F) {
      any_culled = true;
    }
  }
  EXPECT_TRUE(any_culled);

  controller.resetSlabClip();
  EXPECT_FALSE(controller.slabClipConfig().enabled);
  for (const auto &atom : controller.scene()->atoms()) {
    EXPECT_GT(atom.dim_factor, 0.0F);
  }
}

TEST_F(Structure3DControllerTests, RadialBrushAndCameraHudUpdateWindowState) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);
  controller.setCell(createTestCell());

  controller.handleRotate(15.0F, 10.0F);
  const auto hud = window().get_structure_3d_camera_hud();
  EXPECT_GT(hud.scale_bar_px, 0.0F);
  EXPECT_FALSE(std::string_view(hud.scale_bar_label.data()).empty());

  controller.setRadialBrushDistance(1.5);
  EXPECT_FALSE(std::string_view(window().get_structure_3d_radial_brush_text().data()).empty());

  controller.setRadialBrushDistance(std::nullopt);
  EXPECT_TRUE(std::string_view(window().get_structure_3d_radial_brush_text().data()).empty());
}

TEST_F(Structure3DControllerTests, DepthCueingToggleAndElementVisibility) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);
  controller.setCell(createTestCell());

  EXPECT_TRUE(controller.isDepthCueingEnabled());
  EXPECT_TRUE(window().get_structure_3d_depth_cueing());
  controller.handleToggleDepthCueing(false);
  EXPECT_FALSE(controller.isDepthCueingEnabled());
  EXPECT_FALSE(window().get_structure_3d_depth_cueing());

  // Element legend populated with O (2) and Si (1)
  const auto legend = window().get_structure_3d_element_legend();
  ASSERT_NE(legend, nullptr);
  EXPECT_EQ(legend->row_count(), 2UL);
  EXPECT_FALSE(window().get_structure_3d_has_hidden_elements());

  // Hide O atoms; Si (index 0) stays visible, O (indices 1, 2) become hidden
  controller.handleToggleElementVisibility("O");
  EXPECT_TRUE(window().get_structure_3d_has_hidden_elements());
  EXPECT_GT(controller.scene()->atoms()[0].dim_factor, 0.0F);
  EXPECT_FLOAT_EQ(controller.scene()->atoms()[1].dim_factor, 0.0F);
  EXPECT_FLOAT_EQ(controller.scene()->atoms()[2].dim_factor, 0.0F);

  // Restore all elements
  controller.handleShowAllElements();
  EXPECT_FALSE(window().get_structure_3d_has_hidden_elements());
  EXPECT_GT(controller.scene()->atoms()[1].dim_factor, 0.0F);
}

TEST_F(Structure3DControllerTests, HiResSnapshotExportWritesValidPng) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);
  controller.handleViewportResized(320.0F, 240.0F);
  controller.setCell(createTestCell());

  const auto tmp_png = std::filesystem::temp_directory_path() / "correlation_3d_snap_test.png";
  std::error_code err_code;
  std::filesystem::remove(tmp_png, err_code);

  std::uint32_t out_w = 0;
  std::uint32_t out_h = 0;
  EXPECT_TRUE(controller.exportSnapshotPng(tmp_png.string(), 0, 1, &out_w, &out_h));
  EXPECT_EQ(out_w, 320U);
  EXPECT_EQ(out_h, 240U);

  std::ifstream in_file(tmp_png, std::ios::binary);
  ASSERT_TRUE(in_file.is_open());
  std::array<char, 4> magic{};
  in_file.read(magic.data(), static_cast<std::streamsize>(magic.size()));
  EXPECT_EQ(static_cast<unsigned char>(magic[0]), 0x89U);
  EXPECT_EQ(magic[1], 'P');
  EXPECT_EQ(magic[2], 'N');
  EXPECT_EQ(magic[3], 'G');
  in_file.close();
  EXPECT_TRUE(std::filesystem::remove(tmp_png, err_code));
}

} // namespace
} // namespace correlation::app
