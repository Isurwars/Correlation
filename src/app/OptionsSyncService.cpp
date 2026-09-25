/**
 * @file OptionsSyncService.cpp
 * @brief Implementation of OptionsSyncService.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/OptionsSyncService.hpp"
#include "app/AppBackend.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <format>
#include <memory>
#include <string>
#include <type_traits>

namespace correlation::app {

namespace {

template <typename T> T safeParse(const slint::SharedString &str, T default_value) {
  try {
    if constexpr (std::is_same_v<T, float>) {
      return std::stof(str.data());
    } else if constexpr (std::is_same_v<T, real_t>) {
      return std::stod(str.data());
    } else {
      return default_value;
    }
  } catch (const std::exception &) {
    return default_value;
  }
}

[[nodiscard]] std::string toLowerStr(std::string_view str) {
  std::string data(str);
  std::ranges::transform(data, data.begin(),
                         [](unsigned char chr) { return static_cast<char>(std::tolower(chr)); });
  return data;
}

void collectActiveCalculators(const std::shared_ptr<slint::Model<CalculatorGroup>> &groups,
                              ProgramOptions &opt) {
  if (groups == nullptr) {
    return;
  }
  for (size_t gi = 0; gi < groups->row_count(); ++gi) {
    const auto maybe_group = groups->row_data(gi);
    if (!maybe_group.has_value() || maybe_group->calculators == nullptr) {
      continue;
    }
    const auto &group = maybe_group.value();
    for (size_t ci = 0; ci < group.calculators->row_count(); ++ci) {
      const auto maybe_calc = group.calculators->row_data(ci);
      if (!maybe_calc.has_value()) {
        continue;
      }
      const auto &calc = maybe_calc.value();
      opt.active_calculators[std::string(calc.id.data())] = calc.enabled;
    }
  }
}

void parseFrameSelection(const AnalysisOptions &ui_opts, size_t frame_count, ProgramOptions &opt) {
  try {
    const std::string min_s = ui_opts.min_frame.data();
    const std::string min_s_lower = toLowerStr(min_s);
    if (min_s_lower == "start" || min_s.empty()) {
      opt.min_frame = 0;
    } else if (min_s_lower == "end") {
      opt.min_frame = std::max(0, static_cast<int>(frame_count) - 1);
    } else {
      opt.min_frame = std::max(0, std::stoi(min_s) - 1);
    }
  } catch (const std::exception &) {
    opt.min_frame = 0;
  }

  try {
    const std::string max_s = ui_opts.max_frame.data();
    const std::string max_s_lower = toLowerStr(max_s);
    if (max_s_lower == "end" || max_s.empty()) {
      opt.max_frame = -1;
    } else if (max_s_lower == "start") {
      opt.max_frame = 1;
    } else {
      opt.max_frame = std::stoi(max_s);
    }
  } catch (const std::exception &) {
    opt.max_frame = -1;
  }

  opt.time_step = safeParse(ui_opts.time_step, opt.time_step);

  try {
    const std::string stride_s = ui_opts.frame_stride.data();
    opt.frame_stride = std::max(1, std::stoi(stride_s));
  } catch (const std::exception &) {
    opt.frame_stride = 1;
  }
}

} // namespace

void OptionsSyncService::updateActiveGroupFlags(AppWindow &window, const ProgramOptions &opts) {
  const auto &calculators =
      ::correlation::calculators::CalculatorFactory::instance().getCalculators();

  bool has_radial = false;
  bool has_scattering = false;
  bool has_angular = false;
  bool has_rings = false;

  for (const auto &calc : calculators) {
    const std::string_view grp = calc->getGroup();
    bool enabled = true;
    auto calc_iter = opts.active_calculators.find(std::string(calc->getName()));
    if (calc_iter != opts.active_calculators.end()) {
      enabled = calc_iter->second;
    }
    if (enabled) {
      if (grp == "Radial") {
        has_radial = true;
      } else if (grp == "Scattering") {
        has_scattering = true;
      } else if (grp == "Angular") {
        has_angular = true;
      } else if (grp == "Rings") {
        has_rings = true;
      }
    }
  }

  window.set_has_radial_active(has_radial);
  window.set_has_scattering_active(has_scattering);
  window.set_has_angular_active(has_angular);
  window.set_has_rings_active(has_rings);
}

void OptionsSyncService::writeToUI(AppWindow &window, const ProgramOptions &options,
                                   const AppBackend & /*backend*/) {
  window.set_in_file_text(slint::SharedString(options.input_file));

  auto opts = window.get_analysis_options();
  opts.smoothing_enabled = options.smoothing;
  opts.r_max = slint::SharedString(std::format("{:.2f}", options.r_max));
  opts.r_bin_width = slint::SharedString(std::format("{:.2f}", options.r_bin_width));
  opts.q_max = slint::SharedString(std::format("{:.2f}", options.q_max));
  opts.q_bin_width = slint::SharedString(std::format("{:.2f}", options.q_bin_width));
  opts.r_int_max = slint::SharedString(std::format("{:.2f}", options.r_int_max));
  opts.angle_bin_width = slint::SharedString(std::format("{:.2f}", options.angle_bin_width));
  opts.dihedral_bin_width = slint::SharedString(std::format("{:.2f}", options.dihedral_bin_width));
  opts.max_ring_size = slint::SharedString(std::to_string(options.max_ring_size));
  opts.smoothing_sigma = slint::SharedString(std::format("{:.2f}", options.smoothing_sigma));
  opts.smoothing_kernel = static_cast<int>(options.smoothing_kernel);
  opts.material_type = options.material_type;
  opts.lef_cutoff = slint::SharedString(std::format("{:.2f}", options.lef_cutoff));
  opts.lef_sigma = slint::SharedString(std::format("{:.2f}", options.lef_sigma));
  opts.hyper_samples = slint::SharedString(std::to_string(options.hyper_samples));
  opts.xrd_radiation_preset = 0;
  opts.xrd_lambda = slint::SharedString(std::format("{:.4f}", options.xrd_params.lambda));
  opts.xrd_theta_min = slint::SharedString(std::format("{:.1f}", options.xrd_params.theta_min));
  opts.xrd_theta_max = slint::SharedString(std::format("{:.1f}", options.xrd_params.theta_max));
  opts.xrd_bin_width = slint::SharedString(std::format("{:.2f}", options.xrd_params.bin_width));

  opts.min_frame = slint::SharedString(std::to_string(options.min_frame + 1));
  if (options.max_frame == -1) {
    opts.max_frame = "End";
  } else {
    opts.max_frame = slint::SharedString(std::to_string(options.max_frame));
  }
  opts.time_step = slint::SharedString(std::format("{:.2f}", options.time_step));
  opts.frame_stride = slint::SharedString(std::to_string(options.frame_stride));

  window.set_analysis_options(opts);
  updateActiveGroupFlags(window, options);
}

std::expected<ProgramOptions, std::string>
OptionsSyncService::readFromUI(const AppWindow &window, size_t frame_count,
                               const correlation::analysis::BondCutoffMatrix &bond_cutoffs) {
  ProgramOptions opt;
  const std::string input_path_str = window.get_in_file_text().data();
  const std::filesystem::path full_path(input_path_str);
  std::filesystem::path output_path = full_path.parent_path() / full_path.stem();
  opt.input_file = input_path_str;
  opt.output_file_base = output_path.make_preferred().string();
  opt.smoothing = true;

  const auto ui_opts = window.get_analysis_options();
  opt.r_max = safeParse(ui_opts.r_max, opt.r_max);
  opt.r_bin_width = safeParse(ui_opts.r_bin_width, opt.r_bin_width);
  opt.q_max = safeParse(ui_opts.q_max, opt.q_max);
  opt.q_bin_width = safeParse(ui_opts.q_bin_width, opt.q_bin_width);
  opt.r_int_max = safeParse(ui_opts.r_int_max, opt.r_int_max);
  opt.angle_bin_width = safeParse(ui_opts.angle_bin_width, opt.angle_bin_width);
  opt.dihedral_bin_width = safeParse(ui_opts.dihedral_bin_width, opt.dihedral_bin_width);
  opt.max_ring_size =
      static_cast<size_t>(safeParse(ui_opts.max_ring_size, static_cast<real_t>(opt.max_ring_size)));
  opt.hyper_samples =
      static_cast<size_t>(safeParse(ui_opts.hyper_samples, static_cast<real_t>(opt.hyper_samples)));

  collectActiveCalculators(window.get_calculator_groups(), opt);

  opt.smoothing_sigma = safeParse(ui_opts.smoothing_sigma, opt.smoothing_sigma);
  opt.smoothing_kernel = static_cast<correlation::math::KernelType>(ui_opts.smoothing_kernel);
  opt.material_type = ui_opts.material_type;
  opt.lef_cutoff = safeParse(ui_opts.lef_cutoff, opt.lef_cutoff);
  opt.lef_sigma = safeParse(ui_opts.lef_sigma, opt.lef_sigma);
  opt.xrd_params.lambda = safeParse(ui_opts.xrd_lambda, opt.xrd_params.lambda);
  opt.xrd_params.theta_min = safeParse(ui_opts.xrd_theta_min, opt.xrd_params.theta_min);
  opt.xrd_params.theta_max = safeParse(ui_opts.xrd_theta_max, opt.xrd_params.theta_max);
  opt.xrd_params.bin_width = safeParse(ui_opts.xrd_bin_width, opt.xrd_params.bin_width);

  parseFrameSelection(ui_opts, frame_count, opt);

  if (opt.max_frame >= 0 && opt.min_frame >= opt.max_frame) {
    return std::unexpected("End frame must be ≥ Start frame");
  }
  if (opt.r_max > 0.0 && opt.r_bin_width > opt.r_max) {
    return std::unexpected("r_bin_width must be ≤ r_max");
  }
  if (opt.q_max > 0.0 && opt.q_bin_width > opt.q_max) {
    return std::unexpected("q_bin_width must be ≤ q_max");
  }

  opt.bond_cutoffs = bond_cutoffs;
  return opt;
}

} // namespace correlation::app
