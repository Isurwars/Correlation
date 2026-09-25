/**
 * @file OptionsResetService.cpp
 * @brief Implementation of OptionsResetService.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/OptionsResetService.hpp"

#include <format>
#include <string>

namespace correlation::app {

void OptionsResetService::resetRDF(AppWindow &window) {
  auto opts = window.get_analysis_options();
  opts.r_max = slint::SharedString(std::format("{:.2f}", AppDefaults::R_MAX));
  if (opts.material_type == 2) {
    opts.r_bin_width = slint::SharedString(std::format("{:.3f}", AppDefaults::R_BIN_WIDTH_CRYSTAL));
  } else if (opts.material_type == 1) {
    opts.r_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::R_BIN_WIDTH_LIQUID));
  } else {
    opts.r_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::R_BIN_WIDTH));
  }
  window.set_analysis_options(opts);
}

void OptionsResetService::resetAngle(AppWindow &window) {
  auto opts = window.get_analysis_options();
  if (opts.material_type == 2) {
    opts.angle_bin_width =
        slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH_CRYSTAL));
    opts.dihedral_bin_width =
        slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH_CRYSTAL));
  } else if (opts.material_type == 1) {
    opts.angle_bin_width =
        slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH_LIQUID));
    opts.dihedral_bin_width =
        slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH_LIQUID));
  } else {
    opts.angle_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH));
    opts.dihedral_bin_width =
        slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH));
  }
  window.set_analysis_options(opts);
}

void OptionsResetService::resetSQ(AppWindow &window) {
  auto opts = window.get_analysis_options();
  opts.q_max = slint::SharedString(std::format("{:.2f}", AppDefaults::Q_MAX));
  opts.r_int_max = slint::SharedString(std::format("{:.2f}", AppDefaults::R_INT_MAX));
  if (opts.material_type == 2) {
    opts.q_bin_width = slint::SharedString(std::format("{:.3f}", AppDefaults::Q_BIN_WIDTH_CRYSTAL));
  } else if (opts.material_type == 1) {
    opts.q_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::Q_BIN_WIDTH_LIQUID));
  } else {
    opts.q_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::Q_BIN_WIDTH));
  }
  window.set_analysis_options(opts);
}

void OptionsResetService::resetXRD(AppWindow &window) {
  auto opts = window.get_analysis_options();
  opts.xrd_radiation_preset = 0;
  opts.xrd_lambda = slint::SharedString(std::format("{:.4f}", AppDefaults::XRD_LAMBDA));
  opts.xrd_theta_min = slint::SharedString(std::format("{:.1f}", AppDefaults::XRD_THETA_MIN));
  opts.xrd_theta_max = slint::SharedString(std::format("{:.1f}", AppDefaults::XRD_THETA_MAX));
  opts.xrd_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::XRD_BIN_WIDTH));
  window.set_analysis_options(opts);
}

void OptionsResetService::resetRings(AppWindow &window) {
  auto opts = window.get_analysis_options();
  opts.max_ring_size = slint::SharedString(std::to_string(ProgramOptions{}.max_ring_size));
  window.set_analysis_options(opts);
}

void OptionsResetService::resetSmoothing(AppWindow &window) {
  auto opts = window.get_analysis_options();
  opts.smoothing_enabled = ProgramOptions{}.smoothing;
  opts.smoothing_kernel = static_cast<int>(AppDefaults::SMOOTHING_KERNEL);
  if (opts.material_type == 2) {
    opts.smoothing_sigma =
        slint::SharedString(std::format("{:.2f}", AppDefaults::SMOOTHING_SIGMA_CRYSTAL));
  } else if (opts.material_type == 1) {
    opts.smoothing_sigma =
        slint::SharedString(std::format("{:.2f}", AppDefaults::SMOOTHING_SIGMA_LIQUID));
  } else {
    opts.smoothing_sigma = slint::SharedString(std::format("{:.2f}", AppDefaults::SMOOTHING_SIGMA));
  }
  window.set_analysis_options(opts);
}

void OptionsResetService::resetAdvanced(AppWindow &window) {
  auto opts = window.get_analysis_options();
  opts.lef_cutoff = slint::SharedString(std::format("{:.2f}", AppDefaults::LEF_CUTOFF));
  opts.lef_sigma = slint::SharedString(std::format("{:.2f}", AppDefaults::LEF_SIGMA));
  opts.hyper_samples = slint::SharedString(std::to_string(ProgramOptions{}.hyper_samples));
  window.set_analysis_options(opts);
}

void OptionsResetService::resetTrajectory(AppWindow &window, const AppBackend &backend) {
  auto opts = window.get_analysis_options();
  if (backend.getFrameCount() > 0) {
    opts.time_step = slint::SharedString(std::format("{:.2f}", backend.getRecommendedTimeStep()));
    opts.min_frame = "1";
    opts.max_frame = slint::SharedString(std::to_string(backend.getFrameCount()));
  } else {
    opts.time_step = slint::SharedString(std::format("{:.2f}", AppDefaults::TIME_STEP));
    opts.min_frame = "1";
    opts.max_frame = "End";
  }
  opts.frame_stride = "1";
  window.set_analysis_options(opts);
}

void OptionsResetService::resetExportSettings(AppWindow &window) {
  ExportConfig cfg;
  cfg.size_preset = 0;
  cfg.palette = 0;
  cfg.font_scale = "1.0";
  cfg.line_width = "3.0";
  cfg.marker_size = "3.5";
  cfg.show_legend = true;
  cfg.show_grid = true;
  cfg.show_markers = false;
  cfg.fill_area = false;
  window.set_export_config(cfg);
}

void OptionsResetService::resetMaterialType(AppWindow &window) {
  auto opts = window.get_analysis_options();
  opts.material_type = 0;
  window.set_analysis_options(opts);
}

} // namespace correlation::app
