/**
 * @file InputValidator.cpp
 * @brief Implementation of InputValidator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#if defined(_WIN32)
#define NOMINMAX
#include <Windows.h>
#endif

#include "AppWindow.h"
#include "app/AppController.hpp"
#include "app/InputValidator.hpp"

#include <algorithm>
#include <cctype>
#include <format>
#include <string>

namespace correlation::app {

namespace {

std::string to_lower_str(const std::string &str) {
  std::string data = str;
  std::ranges::transform(data, data.begin(), [](const unsigned char chr) { return std::tolower(chr); });
  return data;
}

bool is_positive_float(const std::string &str, float &val) {
  try {
    size_t idx = 0;
    float parsed_val = std::stof(str, &idx);
    if (idx < str.size() || parsed_val <= 0.0F) {
      return false;
    }
    val = parsed_val;
    return true;
  } catch (...) {
    return false;
  }
}

bool is_non_negative_float(const std::string &str, float &val) {
  try {
    size_t idx = 0;
    float parsed_val = std::stof(str, &idx);
    if (idx < str.size() || parsed_val < 0.0F) {
      return false;
    }
    val = parsed_val;
    return true;
  } catch (...) {
    return false;
  }
}

bool is_positive_int(const std::string &str, int &val) {
  try {
    size_t idx = 0;
    int parsed_val = std::stoi(str, &idx);
    if (idx < str.size() || parsed_val <= 0) {
      return false;
    }
    val = parsed_val;
    return true;
  } catch (...) {
    return false;
  }
}

bool validate_min_frame_input(const std::string &frame_s, int total_frames, int &frame_val,
                              slint::SharedString &error_out) {
  std::string frame_lower = to_lower_str(frame_s);

  if (frame_lower == "start" || frame_s.empty()) {
    frame_val = 0;
    return true;
  }

  if (frame_lower == "end") {
    frame_val = total_frames > 0 ? total_frames - 1 : 0;
    return true;
  }

  int parsed_val = 0;
  if (!is_positive_int(frame_s, parsed_val)) {
    error_out = "Must be positive integer, 'Start', or 'End'";
    return false;
  }

  if (total_frames > 0 && parsed_val > total_frames) {
    error_out = slint::SharedString(std::format("Must be ≤ total frames ({})", total_frames));
    return false;
  }

  frame_val = parsed_val - 1;
  return true;
}

bool validate_max_frame_input(const std::string &frame_s, int total_frames, int &frame_val,
                              slint::SharedString &error_out) {
  std::string frame_lower = to_lower_str(frame_s);

  if (frame_lower == "end" || frame_s.empty()) {
    frame_val = total_frames > 0 ? total_frames - 1 : -1;
    return true;
  }

  if (frame_lower == "start") {
    frame_val = 0;
    return true;
  }

  int parsed_val = 0;
  if (!is_positive_int(frame_s, parsed_val)) {
    error_out = "Must be positive integer or 'End'";
    return false;
  }

  if (total_frames > 0 && parsed_val > total_frames) {
    error_out = slint::SharedString(std::format("Must be ≤ total frames ({})", total_frames));
    return false;
  }

  frame_val = parsed_val - 1;
  return true;
}

} // namespace

InputValidator::InputValidator(::AppWindow &window, AppBackend &backend, AppController &controller)
    : window_(&window), backend_(&backend), controller_(&controller) {}

bool InputValidator::validateRadialAndScattering(AppErrors &errs, float &r_max_val, float &q_max_val) {
  bool valid = true;

  std::string r_max_s = window_->get_analysis_options().r_max.data();
  if (!is_positive_float(r_max_s, r_max_val)) {
    errs.r_max_error = "Must be a positive number";
    valid = false;
  }

  float r_bin_val = 0.0F;
  std::string r_bin_s = window_->get_analysis_options().r_bin_width.data();
  if (!is_positive_float(r_bin_s, r_bin_val)) {
    errs.r_bin_error = "Must be a positive number";
    valid = false;
  } else if (r_max_val > 0.0F && r_bin_val > r_max_val) {
    errs.r_bin_error = "Must be ≤ r_max";
    valid = false;
  }

  std::string q_max_s = window_->get_analysis_options().q_max.data();
  if (!is_positive_float(q_max_s, q_max_val)) {
    errs.q_max_error = "Must be a positive number";
    valid = false;
  }

  float q_bin_val = 0.0F;
  std::string q_bin_s = window_->get_analysis_options().q_bin_width.data();
  if (!is_positive_float(q_bin_s, q_bin_val)) {
    errs.q_bin_error = "Must be a positive number";
    valid = false;
  } else if (q_max_val > 0.0F && q_bin_val > q_max_val) {
    errs.q_bin_error = "Must be ≤ q_max";
    valid = false;
  }

  float r_int_max_val = 0.0F;
  std::string r_int_max_s = window_->get_analysis_options().r_int_max.data();
  if (!is_positive_float(r_int_max_s, r_int_max_val)) {
    errs.r_int_max_error = "Must be a positive number";
    valid = false;
  }

  return valid;
}

bool InputValidator::validateAngularAndRings(AppErrors &errs) {
  bool valid = true;

  float angle_bin_val = 0.0F;
  std::string angle_bin_s = window_->get_analysis_options().angle_bin_width.data();
  if (!is_positive_float(angle_bin_s, angle_bin_val)) {
    errs.angle_bin_error = "Must be a positive number";
    valid = false;
  } else if (angle_bin_val > 180.0F) {
    errs.angle_bin_error = "Must be ≤ 180°";
    valid = false;
  }

  float dihedral_bin_val = 0.0F;
  std::string dihedral_bin_s = window_->get_analysis_options().dihedral_bin_width.data();
  if (!is_positive_float(dihedral_bin_s, dihedral_bin_val)) {
    errs.dihedral_bin_error = "Must be a positive number";
    valid = false;
  } else if (dihedral_bin_val > 360.0F) {
    errs.dihedral_bin_error = "Must be ≤ 360°";
    valid = false;
  }

  int max_ring_val = 0;
  std::string max_ring_s = window_->get_analysis_options().max_ring_size.data();
  if (!is_positive_int(max_ring_s, max_ring_val) || max_ring_val < 3) {
    errs.max_ring_error = "Must be an integer ≥ 3";
    valid = false;
  }

  return valid;
}

bool InputValidator::validateOtherAnalysisOptions(AppErrors &errs) {
  bool valid = true;

  float smoothing_sigma_val = 0.0F;
  std::string smoothing_sigma_s = window_->get_analysis_options().smoothing_sigma.data();
  if (!is_non_negative_float(smoothing_sigma_s, smoothing_sigma_val)) {
    errs.smoothing_sigma_error = "Must be a non-negative number";
    valid = false;
  }

  float time_step_val = 0.0F;
  std::string time_step_s = window_->get_analysis_options().time_step.data();
  if (!is_positive_float(time_step_s, time_step_val)) {
    errs.time_step_error = "Must be a positive number";
    valid = false;
  }

  float lef_cutoff_val = 0.0F;
  std::string lef_cutoff_s = window_->get_analysis_options().lef_cutoff.data();
  if (!is_positive_float(lef_cutoff_s, lef_cutoff_val)) {
    errs.lef_cutoff_error = "Must be a positive number";
    valid = false;
  }

  float lef_sigma_val = 0.0F;
  std::string lef_sigma_s = window_->get_analysis_options().lef_sigma.data();
  if (!is_positive_float(lef_sigma_s, lef_sigma_val)) {
    errs.lef_sigma_error = "Must be a positive number";
    valid = false;
  }

  int hyper_samples_val = 0;
  std::string hyper_samples_s = window_->get_analysis_options().hyper_samples.data();
  if (!is_positive_int(hyper_samples_s, hyper_samples_val)) {
    errs.hyper_samples_error = "Must be a positive integer";
    valid = false;
  }

  return valid;
}

bool InputValidator::validateFrames(AppErrors &errs) {
  bool valid = true;

  int min_frame_val = -1;
  int max_frame_val = -1;
  int total_frames = window_->get_num_frames();

  std::string min_frame_s = window_->get_analysis_options().min_frame.data();
  bool min_frame_valid = validate_min_frame_input(min_frame_s, total_frames, min_frame_val, errs.min_frame_error);
  if (!min_frame_valid) {
    valid = false;
  }

  std::string max_frame_s = window_->get_analysis_options().max_frame.data();
  bool max_frame_valid = validate_max_frame_input(max_frame_s, total_frames, max_frame_val, errs.max_frame_error);
  if (!max_frame_valid) {
    valid = false;
  }

  if (min_frame_valid && max_frame_valid && min_frame_val >= 0 && max_frame_val >= 0) {
    if (min_frame_val > max_frame_val) {
      errs.min_frame_error = "Start frame must be ≤ End frame";
      errs.max_frame_error = "End frame must be ≥ Start frame";
      valid = false;
    }
  }

  return valid;
}

bool InputValidator::validateExportConfig(AppErrors &errs) {
  bool valid = true;

  float font_scale_val = 0.0F;
  std::string font_scale_s = window_->get_export_config().font_scale.data();
  if (!is_positive_float(font_scale_s, font_scale_val)) {
    errs.export_font_scale_error = "Must be a positive number";
    valid = false;
  }

  float line_width_val = 0.0F;
  std::string line_width_s = window_->get_export_config().line_width.data();
  if (!is_positive_float(line_width_s, line_width_val)) {
    errs.export_line_width_error = "Must be a positive number";
    valid = false;
  }

  return valid;
}

bool InputValidator::validateInputs() {
  bool valid = true;
  auto errs = window_->get_app_errors();

  errs.r_max_error = "";
  errs.r_bin_error = "";
  errs.q_max_error = "";
  errs.q_bin_error = "";
  errs.r_int_max_error = "";
  errs.angle_bin_error = "";
  errs.dihedral_bin_error = "";
  errs.max_ring_error = "";
  errs.smoothing_sigma_error = "";
  errs.time_step_error = "";
  errs.min_frame_error = "";
  errs.max_frame_error = "";
  errs.export_font_scale_error = "";
  errs.export_line_width_error = "";
  errs.lef_cutoff_error = "";
  errs.lef_sigma_error = "";

  float r_max_val = 0.0F;
  float q_max_val = 0.0F;

  if (!validateRadialAndScattering(errs, r_max_val, q_max_val)) {
    valid = false;
  }
  if (!validateAngularAndRings(errs)) {
    valid = false;
  }
  if (!validateOtherAnalysisOptions(errs)) {
    valid = false;
  }
  if (!validateFrames(errs)) {
    valid = false;
  }
  if (!validateExportConfig(errs)) {
    valid = false;
  }

  bool has_errors = !valid;
  if (window_->get_has_validation_errors() != has_errors) {
    window_->set_has_validation_errors(has_errors);
  }

  auto current_errs = window_->get_app_errors();
  if (current_errs.r_max_error != errs.r_max_error || current_errs.r_bin_error != errs.r_bin_error ||
      current_errs.q_max_error != errs.q_max_error || current_errs.q_bin_error != errs.q_bin_error ||
      current_errs.r_int_max_error != errs.r_int_max_error || current_errs.angle_bin_error != errs.angle_bin_error ||
      current_errs.dihedral_bin_error != errs.dihedral_bin_error ||
      current_errs.max_ring_error != errs.max_ring_error ||
      current_errs.smoothing_sigma_error != errs.smoothing_sigma_error ||
      current_errs.time_step_error != errs.time_step_error || current_errs.min_frame_error != errs.min_frame_error ||
      current_errs.max_frame_error != errs.max_frame_error ||
      current_errs.export_font_scale_error != errs.export_font_scale_error ||
      current_errs.export_line_width_error != errs.export_line_width_error ||
      current_errs.lef_cutoff_error != errs.lef_cutoff_error || current_errs.lef_sigma_error != errs.lef_sigma_error) {
    window_->set_app_errors(errs);
  }

  return valid;
}

} // namespace correlation::app
