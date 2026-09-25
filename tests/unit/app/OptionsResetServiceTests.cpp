/**
 * @file OptionsResetServiceTests.cpp
 * @brief Unit tests for OptionsResetService.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "AppWindow.h"
#include "app/services/OptionsResetService.hpp"
#include "app/services/TrajectoryLoader.hpp"

#include <gtest/gtest.h>
#include <optional>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

namespace correlation::app {
namespace {

class OptionsResetServiceTests : public ::testing::Test {
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

TEST_F(OptionsResetServiceTests, ResetRDFAppliesMaterialSpecificDefaults) {
  auto &win = window();

  // Amorphous / Default (material_type == 0)
  auto opts = win.get_analysis_options();
  opts.material_type = 0;
  opts.r_max = "100.0";
  opts.r_bin_width = "5.0";
  win.set_analysis_options(opts);

  OptionsResetService::resetRDF(win);
  opts = win.get_analysis_options();
  EXPECT_EQ(opts.r_max, "20.00");
  EXPECT_EQ(opts.r_bin_width, "0.02");

  // Liquid (material_type == 1)
  opts.material_type = 1;
  win.set_analysis_options(opts);
  OptionsResetService::resetRDF(win);
  opts = win.get_analysis_options();
  EXPECT_EQ(opts.r_bin_width, "0.05");

  // Crystal (material_type == 2)
  opts.material_type = 2;
  win.set_analysis_options(opts);
  OptionsResetService::resetRDF(win);
  opts = win.get_analysis_options();
  EXPECT_EQ(opts.r_bin_width, "0.002");
}

TEST_F(OptionsResetServiceTests, ResetAngleAppliesDefaults) {
  auto &win = window();

  auto opts = win.get_analysis_options();
  opts.material_type = 0;
  win.set_analysis_options(opts);

  OptionsResetService::resetAngle(win);
  opts = win.get_analysis_options();
  EXPECT_EQ(opts.angle_bin_width, "0.25");
  EXPECT_EQ(opts.dihedral_bin_width, "0.25");

  // Crystal
  opts.material_type = 2;
  win.set_analysis_options(opts);
  OptionsResetService::resetAngle(win);
  opts = win.get_analysis_options();
  EXPECT_EQ(opts.angle_bin_width, "0.10");
  EXPECT_EQ(opts.dihedral_bin_width, "0.10");
}

TEST_F(OptionsResetServiceTests, ResetSQAppliesDefaults) {
  auto &win = window();

  auto opts = win.get_analysis_options();
  opts.material_type = 0;
  win.set_analysis_options(opts);

  OptionsResetService::resetSQ(win);
  opts = win.get_analysis_options();
  EXPECT_EQ(opts.q_max, "20.00");
  EXPECT_EQ(opts.r_int_max, "10.00");
  EXPECT_EQ(opts.q_bin_width, "0.02");

  opts.material_type = 2;
  win.set_analysis_options(opts);
  OptionsResetService::resetSQ(win);
  opts = win.get_analysis_options();
  EXPECT_EQ(opts.q_bin_width, "0.002");
}

TEST_F(OptionsResetServiceTests, ResetXRDAppliesDefaults) {
  auto &win = window();

  OptionsResetService::resetXRD(win);
  const auto opts = win.get_analysis_options();
  EXPECT_EQ(opts.xrd_radiation_preset, 0);
  EXPECT_EQ(opts.xrd_lambda, "1.5406");
  EXPECT_EQ(opts.xrd_theta_min, "10.0");
  EXPECT_EQ(opts.xrd_theta_max, "140.0");
  EXPECT_EQ(opts.xrd_bin_width, "0.05");
}

TEST_F(OptionsResetServiceTests, ResetRingsAndSmoothingAndAdvanced) {
  auto &win = window();

  OptionsResetService::resetRings(win);
  auto opts = win.get_analysis_options();
  EXPECT_EQ(opts.max_ring_size, "8");

  OptionsResetService::resetSmoothing(win);
  opts = win.get_analysis_options();
  EXPECT_TRUE(opts.smoothing_enabled);
  EXPECT_EQ(opts.smoothing_kernel, 0);
  EXPECT_EQ(opts.smoothing_sigma, "0.10");

  OptionsResetService::resetAdvanced(win);
  opts = win.get_analysis_options();
  EXPECT_EQ(opts.lef_cutoff, "5.00");
  EXPECT_EQ(opts.lef_sigma, "0.20");
  EXPECT_EQ(opts.hyper_samples, "10000");
}

TEST_F(OptionsResetServiceTests, ResetTrajectoryWithZeroAndNonZeroFrames) {
  auto &win = window();
  const correlation::app::TrajectoryLoader loader;

  // Frame count 0
  OptionsResetService::resetTrajectory(win, loader);
  auto opts = win.get_analysis_options();
  EXPECT_EQ(opts.min_frame, "1");
  EXPECT_EQ(opts.max_frame, "End");
  EXPECT_EQ(opts.frame_stride, "1");
}

TEST_F(OptionsResetServiceTests, ResetExportSettingsAndMaterialType) {
  auto &win = window();

  OptionsResetService::resetExportSettings(win);
  const auto cfg = win.get_export_config();
  EXPECT_EQ(cfg.size_preset, 0);
  EXPECT_EQ(cfg.palette, 0);
  EXPECT_EQ(cfg.font_scale, "1.0");
  EXPECT_TRUE(cfg.show_legend);
  EXPECT_TRUE(cfg.show_grid);
  EXPECT_FALSE(cfg.show_markers);

  auto opts = win.get_analysis_options();
  opts.material_type = 2;
  win.set_analysis_options(opts);
  OptionsResetService::resetMaterialType(win);
  opts = win.get_analysis_options();
  EXPECT_EQ(opts.material_type, 0);
}

} // namespace
} // namespace correlation::app
