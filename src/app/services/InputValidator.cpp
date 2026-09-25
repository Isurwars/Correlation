/**
 * @file InputValidator.cpp
 * @brief Implementation of InputValidator coordinating UI validation rules.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

#include "AppWindow.h"
#include "app/AppController.hpp"
#include "app/InputValidator.hpp"
#include "app/ValidationRuleService.hpp"

namespace correlation::app {

InputValidator::InputValidator(::AppWindow &window, AppController &controller)
    : window_(&window), controller_(&controller) {}

InputValidator::InputValidator(::AppWindow &window, AppBackend & /*backend*/,
                               AppController &controller)
    : InputValidator(window, controller) {}

bool InputValidator::validateRadialAndScattering(AppErrors &errs, float &r_max_val,
                                                 float &q_max_val) {
  bool valid = true;
  const auto opts = window_->get_analysis_options();

  if (const auto res = ValidationRuleService::parsePositiveFloat(opts.r_max.data()); res) {
    r_max_val = *res;
  } else {
    errs.r_max_error = slint::SharedString(res.error());
    valid = false;
  }

  if (const auto res = ValidationRuleService::parsePositiveFloat(opts.r_bin_width.data()); res) {
    if (const auto bin_check =
            ValidationRuleService::validateBinWithinMax(*res, r_max_val, "r_max");
        !bin_check) {
      errs.r_bin_error = slint::SharedString(bin_check.error());
      valid = false;
    }
  } else {
    errs.r_bin_error = slint::SharedString(res.error());
    valid = false;
  }

  if (const auto res = ValidationRuleService::parsePositiveFloat(opts.q_max.data()); res) {
    q_max_val = *res;
  } else {
    errs.q_max_error = slint::SharedString(res.error());
    valid = false;
  }

  if (const auto res = ValidationRuleService::parsePositiveFloat(opts.q_bin_width.data()); res) {
    if (const auto bin_check =
            ValidationRuleService::validateBinWithinMax(*res, q_max_val, "q_max");
        !bin_check) {
      errs.q_bin_error = slint::SharedString(bin_check.error());
      valid = false;
    }
  } else {
    errs.q_bin_error = slint::SharedString(res.error());
    valid = false;
  }

  if (const auto res = ValidationRuleService::parsePositiveFloat(opts.r_int_max.data()); !res) {
    errs.r_int_max_error = slint::SharedString(res.error());
    valid = false;
  }

  return valid;
}

bool InputValidator::validateXrdOptions(AppErrors &errs) {
  bool valid = true;
  const auto opts = window_->get_analysis_options();

  if (const auto res = ValidationRuleService::parsePositiveFloat(opts.xrd_lambda.data()); !res) {
    errs.xrd_lambda_error = slint::SharedString(res.error());
    valid = false;
  }

  float theta_min = -1.0F;
  float theta_max = -1.0F;
  if (const auto res = ValidationRuleService::parseNonNegativeFloat(opts.xrd_theta_min.data());
      res) {
    theta_min = *res;
  } else {
    errs.xrd_theta_min_error = slint::SharedString(res.error());
    valid = false;
  }

  if (const auto res = ValidationRuleService::parsePositiveFloat(opts.xrd_theta_max.data()); res) {
    theta_max = *res;
  } else {
    errs.xrd_theta_max_error = slint::SharedString(res.error());
    valid = false;
  }

  if (theta_min >= 0.0F && theta_max > 0.0F) {
    if (const auto theta_check = ValidationRuleService::validateXrdTheta(theta_min, theta_max);
        !theta_check) {
      errs.xrd_theta_max_error = slint::SharedString(theta_check.error());
      valid = false;
    }
  }

  if (const auto res = ValidationRuleService::parsePositiveFloat(opts.xrd_bin_width.data()); !res) {
    errs.xrd_bin_width_error = slint::SharedString(res.error());
    valid = false;
  }

  return valid;
}

bool InputValidator::validateAngularAndRings(AppErrors &errs) {
  bool valid = true;
  const auto opts = window_->get_analysis_options();

  if (const auto res = ValidationRuleService::parsePositiveFloat(opts.angle_bin_width.data());
      res) {
    if (const auto angle_check = ValidationRuleService::validateAngleDegrees(*res); !angle_check) {
      errs.angle_bin_error = slint::SharedString(angle_check.error());
      valid = false;
    }
  } else {
    errs.angle_bin_error = slint::SharedString(res.error());
    valid = false;
  }

  if (const auto res = ValidationRuleService::parsePositiveFloat(opts.dihedral_bin_width.data());
      res) {
    if (*res > 360.0F) {
      errs.dihedral_bin_error = "Must be ≤ 360°";
      valid = false;
    }
  } else {
    errs.dihedral_bin_error = slint::SharedString(res.error());
    valid = false;
  }

  if (const auto res = ValidationRuleService::parsePositiveInt(opts.max_ring_size.data()); res) {
    if (*res < 3) {
      errs.max_ring_error = "Must be an integer ≥ 3";
      valid = false;
    }
  } else {
    errs.max_ring_error = "Must be an integer ≥ 3";
    valid = false;
  }

  return valid;
}

bool InputValidator::validateOtherAnalysisOptions(AppErrors &errs) {
  bool valid = true;
  const auto opts = window_->get_analysis_options();

  if (const auto res = ValidationRuleService::parseNonNegativeFloat(opts.smoothing_sigma.data());
      !res) {
    errs.smoothing_sigma_error = "Must be a non-negative number";
    valid = false;
  }

  if (const auto res = ValidationRuleService::parsePositiveFloat(opts.time_step.data()); !res) {
    errs.time_step_error = slint::SharedString(res.error());
    valid = false;
  }

  if (const auto res = ValidationRuleService::parsePositiveFloat(opts.lef_cutoff.data()); !res) {
    errs.lef_cutoff_error = slint::SharedString(res.error());
    valid = false;
  }

  if (const auto res = ValidationRuleService::parsePositiveFloat(opts.lef_sigma.data()); !res) {
    errs.lef_sigma_error = slint::SharedString(res.error());
    valid = false;
  }

  if (const auto res = ValidationRuleService::parsePositiveInt(opts.hyper_samples.data()); !res) {
    errs.hyper_samples_error = slint::SharedString(res.error());
    valid = false;
  }

  return valid;
}

bool InputValidator::validateFrames(AppErrors &errs) {
  bool valid = true;
  const int total_frames = window_->get_num_frames();
  const auto opts = window_->get_analysis_options();

  int min_frame_val = -1;
  int max_frame_val = -1;

  if (const auto res = ValidationRuleService::parseMinFrame(opts.min_frame.data(), total_frames);
      res) {
    min_frame_val = *res;
  } else {
    errs.min_frame_error = slint::SharedString(res.error());
    valid = false;
  }

  if (const auto res = ValidationRuleService::parseMaxFrame(opts.max_frame.data(), total_frames);
      res) {
    max_frame_val = *res;
  } else {
    errs.max_frame_error = slint::SharedString(res.error());
    valid = false;
  }

  if (min_frame_val >= 0 && max_frame_val >= 0 && min_frame_val > max_frame_val) {
    errs.min_frame_error = "Start frame must be ≤ End frame";
    errs.max_frame_error = "End frame must be ≥ Start frame";
    valid = false;
  }

  if (const auto res = ValidationRuleService::parseFrameStride(opts.frame_stride.data()); !res) {
    errs.frame_stride_error = slint::SharedString(res.error());
    valid = false;
  }

  return valid;
}

bool InputValidator::validateExportConfig(AppErrors &errs) {
  bool valid = true;
  const auto export_cfg = window_->get_export_config();

  if (const auto res = ValidationRuleService::parsePositiveFloat(export_cfg.font_scale.data());
      !res) {
    errs.export_font_scale_error = slint::SharedString(res.error());
    valid = false;
  }

  if (const auto res = ValidationRuleService::parsePositiveFloat(export_cfg.line_width.data());
      !res) {
    errs.export_line_width_error = slint::SharedString(res.error());
    valid = false;
  }

  if (const auto res = ValidationRuleService::parsePositiveFloat(export_cfg.marker_size.data());
      !res) {
    errs.export_marker_size_error = slint::SharedString(res.error());
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
  errs.frame_stride_error = "";
  errs.export_font_scale_error = "";
  errs.export_line_width_error = "";
  errs.export_marker_size_error = "";
  errs.lef_cutoff_error = "";
  errs.lef_sigma_error = "";
  errs.xrd_lambda_error = "";
  errs.xrd_theta_min_error = "";
  errs.xrd_theta_max_error = "";
  errs.xrd_bin_width_error = "";

  float r_max_val = 0.0F;
  float q_max_val = 0.0F;

  if (!validateRadialAndScattering(errs, r_max_val, q_max_val)) {
    valid = false;
  }
  if (!validateXrdOptions(errs)) {
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

  const bool has_errors = !valid;
  if (window_->get_has_validation_errors() != has_errors) {
    window_->set_has_validation_errors(has_errors);
  }

  const auto current_errs = window_->get_app_errors();
  if (current_errs.r_max_error != errs.r_max_error ||
      current_errs.r_bin_error != errs.r_bin_error ||
      current_errs.q_max_error != errs.q_max_error ||
      current_errs.q_bin_error != errs.q_bin_error ||
      current_errs.r_int_max_error != errs.r_int_max_error ||
      current_errs.angle_bin_error != errs.angle_bin_error ||
      current_errs.dihedral_bin_error != errs.dihedral_bin_error ||
      current_errs.max_ring_error != errs.max_ring_error ||
      current_errs.smoothing_sigma_error != errs.smoothing_sigma_error ||
      current_errs.time_step_error != errs.time_step_error ||
      current_errs.min_frame_error != errs.min_frame_error ||
      current_errs.max_frame_error != errs.max_frame_error ||
      current_errs.frame_stride_error != errs.frame_stride_error ||
      current_errs.export_font_scale_error != errs.export_font_scale_error ||
      current_errs.export_line_width_error != errs.export_line_width_error ||
      current_errs.export_marker_size_error != errs.export_marker_size_error ||
      current_errs.lef_cutoff_error != errs.lef_cutoff_error ||
      current_errs.lef_sigma_error != errs.lef_sigma_error ||
      current_errs.xrd_lambda_error != errs.xrd_lambda_error ||
      current_errs.xrd_theta_min_error != errs.xrd_theta_min_error ||
      current_errs.xrd_theta_max_error != errs.xrd_theta_max_error ||
      current_errs.xrd_bin_width_error != errs.xrd_bin_width_error) {
    window_->set_app_errors(errs);
  }

  return valid;
}

} // namespace correlation::app
