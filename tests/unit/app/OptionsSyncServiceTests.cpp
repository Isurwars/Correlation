/**
 * @file OptionsSyncServiceTests.cpp
 * @brief Unit tests for OptionsSyncService.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "AppWindow.h"
#include "app/services/OptionsSyncService.hpp"

#include <gtest/gtest.h>
#include <optional>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

namespace correlation::app {
namespace {

class OptionsSyncServiceTests : public ::testing::Test {
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

TEST_F(OptionsSyncServiceTests, WriteToUIPopulatesAllFieldsCorrectly) {
  auto &win = window();

  ProgramOptions opt;
  opt.input_file = "/tmp/test_structure.xyz";
  opt.smoothing = true;
  opt.r_max = 12.5;
  opt.r_bin_width = 0.04;
  opt.q_max = 30.0;
  opt.q_bin_width = 0.05;
  opt.r_int_max = 22.0;
  opt.angle_bin_width = 0.45;
  opt.dihedral_bin_width = 0.55;
  opt.max_ring_size = 14;
  opt.smoothing_sigma = 0.15;
  opt.smoothing_kernel = correlation::math::KernelType::Gaussian;
  opt.material_type = 1;
  opt.lef_cutoff = 3.5;
  opt.lef_sigma = 0.12;
  opt.hyper_samples = 2500;
  opt.xrd_params.lambda = 1.7890;
  opt.xrd_params.theta_min = 10.0;
  opt.xrd_params.theta_max = 80.0;
  opt.xrd_params.bin_width = 0.08;
  opt.min_frame = 0;
  opt.max_frame = -1;
  opt.time_step = 1.25;
  opt.frame_stride = 2;

  OptionsSyncService::writeToUI(win, opt);

  EXPECT_EQ(win.get_in_file_text(), "/tmp/test_structure.xyz");
  const auto ui_opts = win.get_analysis_options();
  EXPECT_TRUE(ui_opts.smoothing_enabled);
  EXPECT_EQ(ui_opts.r_max, "12.50");
  EXPECT_EQ(ui_opts.r_bin_width, "0.04");
  EXPECT_EQ(ui_opts.q_max, "30.00");
  EXPECT_EQ(ui_opts.q_bin_width, "0.05");
  EXPECT_EQ(ui_opts.r_int_max, "22.00");
  EXPECT_EQ(ui_opts.angle_bin_width, "0.45");
  EXPECT_EQ(ui_opts.dihedral_bin_width, "0.55");
  EXPECT_EQ(ui_opts.max_ring_size, "14");
  EXPECT_EQ(ui_opts.smoothing_sigma, "0.15");
  EXPECT_EQ(ui_opts.smoothing_kernel, 0);
  EXPECT_EQ(ui_opts.material_type, 1);
  EXPECT_EQ(ui_opts.lef_cutoff, "3.50");
  EXPECT_EQ(ui_opts.lef_sigma, "0.12");
  EXPECT_EQ(ui_opts.hyper_samples, "2500");
  EXPECT_EQ(ui_opts.xrd_lambda, "1.7890");
  EXPECT_EQ(ui_opts.xrd_theta_min, "10.0");
  EXPECT_EQ(ui_opts.xrd_theta_max, "80.0");
  EXPECT_EQ(ui_opts.xrd_bin_width, "0.08");
  EXPECT_EQ(ui_opts.min_frame, "1");
  EXPECT_EQ(ui_opts.max_frame, "End");
  EXPECT_EQ(ui_opts.time_step, "1.25");
  EXPECT_EQ(ui_opts.frame_stride, "2");
}

TEST_F(OptionsSyncServiceTests, ReadFromUISuccessAndOutputBaseCalculation) {
  auto &win = window();

  ProgramOptions opt;
  opt.input_file = "/path/to/trajectory.xyz";
  opt.min_frame = 2; // 0-based index 2 -> UI frame 3
  opt.max_frame = 50;
  OptionsSyncService::writeToUI(win, opt);

  const auto res = OptionsSyncService::readFromUI(win, 100);
  ASSERT_TRUE(res.has_value());

  const auto &parsed = *res;
  EXPECT_EQ(parsed.input_file, "/path/to/trajectory.xyz");
  EXPECT_NE(parsed.output_file_base.find("trajectory"), std::string::npos);
  EXPECT_EQ(parsed.min_frame, 2);
  EXPECT_EQ(parsed.max_frame, 50);
}

TEST_F(OptionsSyncServiceTests, ReadFromUIFramePresetsStartAndEnd) {
  auto &win = window();

  auto opts = win.get_analysis_options();
  opts.min_frame = "Start";
  opts.max_frame = "End";
  opts.r_max = "10.0";
  opts.r_bin_width = "0.05";
  opts.q_max = "25.0";
  opts.q_bin_width = "0.05";
  win.set_analysis_options(opts);

  const auto res = OptionsSyncService::readFromUI(win, 50);
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->min_frame, 0);
  EXPECT_EQ(res->max_frame, -1);
}

TEST_F(OptionsSyncServiceTests, ReadFromUIValidationFailsWhenBinWidthExceedsMax) {
  auto &win = window();

  auto opts = win.get_analysis_options();
  opts.r_max = "5.0";
  opts.r_bin_width = "10.0"; // invalid: bin width > max
  opts.q_max = "20.0";
  opts.q_bin_width = "0.1";
  win.set_analysis_options(opts);

  auto res = OptionsSyncService::readFromUI(win, 10);
  EXPECT_FALSE(res.has_value());
  EXPECT_EQ(res.error(), "r_bin_width must be ≤ r_max");

  opts.r_bin_width = "0.05";
  opts.q_bin_width = "25.0"; // invalid: q_bin_width > q_max
  win.set_analysis_options(opts);
  res = OptionsSyncService::readFromUI(win, 10);
  EXPECT_FALSE(res.has_value());
  EXPECT_EQ(res.error(), "q_bin_width must be ≤ q_max");
}

TEST_F(OptionsSyncServiceTests, ReadFromUIValidationFailsWhenFramesInverted) {
  auto &win = window();

  auto opts = win.get_analysis_options();
  opts.r_max = "10.0";
  opts.r_bin_width = "0.05";
  opts.q_max = "20.0";
  opts.q_bin_width = "0.05";
  opts.min_frame = "10";
  opts.max_frame = "5"; // UI max frame 5 < min frame 10
  win.set_analysis_options(opts);

  const auto res = OptionsSyncService::readFromUI(win, 50);
  EXPECT_FALSE(res.has_value());
  EXPECT_EQ(res.error(), "End frame must be ≥ Start frame");
}

TEST_F(OptionsSyncServiceTests, UpdateActiveGroupFlagsSetsBooleans) {
  auto &win = window();

  ProgramOptions opt;
  opt.active_calculators["RDF"] = true;
  opt.active_calculators["SQ"] = true;
  opt.active_calculators["BAD"] = false;
  opt.active_calculators["Rings"] = false;

  OptionsSyncService::updateActiveGroupFlags(win, opt);
  EXPECT_TRUE(win.get_has_radial_active());
  EXPECT_TRUE(win.get_has_scattering_active());
}

} // namespace
} // namespace correlation::app
