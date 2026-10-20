/**
 * @file BondCutoffController.cpp
 * @brief Implementation of BondCutoffController.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/BondCutoffController.hpp"
#include "analysis/CutoffSanityChecker.hpp"
#include "analysis/FirstMinimumDetector.hpp"
#include "app/formatters/BondCutoffMapper.hpp"
#include "app/services/AnalysisDispatcher.hpp"
#include "app/viewmodel/PlotController.hpp"
#include "app/viewmodel/PlotViewportNavigator.hpp"
#include "physics/PhysicalData.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <memory>
#include <string>
#include <vector>

namespace correlation::app {

BondCutoffController::BondCutoffController(AppWindow &window, TrajectoryLoader &loader,
                                           ProgramOptions &options, AnalysisDispatcher *dispatcher)
    : window_(window), loader_(loader), options_(options), dispatcher_(dispatcher) {}

BondCutoffController::BondCutoffController(AppWindow &window, TrajectoryLoader &loader,
                                           [[maybe_unused]] BondCutoffService &cutoff_service,
                                           ProgramOptions &options)
    : BondCutoffController(window, loader, options, nullptr) {}

void BondCutoffController::setBondCutoffs() {
  if (loader_.cell() == nullptr) {
    return;
  }

  const auto &elements = loader_.cell()->elements();
  const auto entries = BondCutoffMapper::createDefaultCutoffEntries(elements);

  const auto slint_cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();
  for (const auto &entry : entries) {
    slint_cutoffs->push_back({
        .element1 = slint::SharedString(entry.element1),
        .element2 = slint::SharedString(entry.element2),
        .min_distance = slint::SharedString(entry.min_distance),
        .max_distance = slint::SharedString(entry.max_distance),
    });
  }

  window_.set_bond_cutoffs(slint_cutoffs);
  window_.set_bond_cutoffs_reset_trigger(window_.get_bond_cutoffs_reset_trigger() + 1);
  options_.bond_cutoffs = getBondCutoffs();
  BondCutoffService::setBondCutoffs(loader_.trajectoryMut(), options_.bond_cutoffs);

  initializeCalibrationPairs();
}

correlation::analysis::BondCutoffMatrix BondCutoffController::getBondCutoffs() const {
  const auto slint_cutoffs = window_.get_bond_cutoffs();
  if (loader_.cell() == nullptr || slint_cutoffs == nullptr) {
    return {};
  }

  const auto &elements = loader_.cell()->elements();
  std::vector<CutoffEntry> entries;
  entries.reserve(slint_cutoffs->row_count());

  for (size_t k = 0; k < slint_cutoffs->row_count(); ++k) {
    auto maybe_item = slint_cutoffs->row_data(k);
    if (!maybe_item.has_value()) {
      continue;
    }
    const auto &item = maybe_item.value();
    entries.push_back(CutoffEntry{
        .element1 = std::string(item.element1.data()),
        .element2 = std::string(item.element2.data()),
        .min_distance = std::string(item.min_distance.data()),
        .max_distance = std::string(item.max_distance.data()),
    });
  }

  return BondCutoffMapper::parseCutoffMatrix(entries, elements);
}

void BondCutoffController::applyScaledCutoffs(float scale_factor) {
  if (loader_.cell() == nullptr || scale_factor <= 0.0F) {
    return;
  }
  options_.bond_cutoffs = BondCutoffService::applyScaledBondCutoffs(
      loader_.trajectoryMut(), static_cast<real_t>(scale_factor));
  syncCutoffsToUI();
}

void BondCutoffController::setUniformCutoff(float max_cutoff) {
  if (loader_.cell() == nullptr || max_cutoff <= 0.0F) {
    return;
  }
  options_.bond_cutoffs = BondCutoffService::setUniformBondCutoff(
      loader_.cell(), loader_.trajectoryMut(), 0.0, static_cast<real_t>(max_cutoff));
  syncCutoffsToUI();
}

namespace {

[[nodiscard]] real_t safeGetCovalentRadius(const std::string &symbol) {
  try {
    return physics::getCovalentRadius(symbol);
  } catch (const std::out_of_range &) {
    return static_cast<real_t>(1.5);
  }
}

struct CovalentCutoffInput {
  std::string symbol_a;
  std::string symbol_b;
  real_t radius_a = 0.0;
  real_t radius_b = 0.0;
  size_t row_idx = 0;
  float factor = 1.0F;
  FactorBound bound = FactorBound::Min;
};

[[nodiscard]] BondCutoff
createCovalentCutoffRow(const CovalentCutoffInput &input,
                        const std::shared_ptr<slint::Model<BondCutoff>> &slint_cutoffs) {
  const real_t sum_radii = input.radius_a + input.radius_b;
  slint::SharedString min_str;
  slint::SharedString max_str;

  if (input.bound == FactorBound::Min) {
    const real_t min_dist = sum_radii * static_cast<real_t>(input.factor);
    min_str = slint::SharedString(std::format("{:.2f}", min_dist));
    if (slint_cutoffs != nullptr && input.row_idx < slint_cutoffs->row_count()) {
      const auto opt = slint_cutoffs->row_data(input.row_idx);
      max_str = opt.has_value() ? opt->max_distance
                                : slint::SharedString(std::format(
                                      "{:.2f}", sum_radii * AppDefaults::BOND_MAX_FACTOR));
    } else {
      max_str =
          slint::SharedString(std::format("{:.2f}", sum_radii * AppDefaults::BOND_MAX_FACTOR));
    }
  } else {
    const real_t max_dist = sum_radii * static_cast<real_t>(input.factor);
    max_str = slint::SharedString(std::format("{:.2f}", max_dist));
    if (slint_cutoffs != nullptr && input.row_idx < slint_cutoffs->row_count()) {
      const auto opt = slint_cutoffs->row_data(input.row_idx);
      min_str = opt.has_value() ? opt->min_distance
                                : slint::SharedString(std::format(
                                      "{:.2f}", sum_radii * AppDefaults::BOND_MIN_FACTOR));
    } else {
      min_str =
          slint::SharedString(std::format("{:.2f}", sum_radii * AppDefaults::BOND_MIN_FACTOR));
    }
  }

  return BondCutoff{
      .element1 = slint::SharedString(input.symbol_a),
      .element2 = slint::SharedString(input.symbol_b),
      .min_distance = min_str,
      .max_distance = max_str,
  };
}

[[nodiscard]] const correlation::analysis::Histogram *
findRdfHistogram(const AnalysisDispatcher *dispatcher) {
  if (dispatcher == nullptr) {
    return nullptr;
  }
  if (const auto *hist = dispatcher->getHistogram("Radial Distribution Function");
      hist != nullptr) {
    return hist;
  }
  if (const auto *hist = dispatcher->getHistogram("RDF"); hist != nullptr) {
    return hist;
  }
  for (const auto &[name, hist_entry] : dispatcher->getHistograms()) {
    if (name.contains("RDF") || name.contains("Distribution")) {
      return &hist_entry;
    }
  }
  return nullptr;
}

[[nodiscard]] std::string resolveTargetPartial(const correlation::analysis::Histogram &hist,
                                               std::string_view label) {
  std::string target(label);
  if (hist.partials.contains(target)) {
    return target;
  }
  const auto hyphen = target.find('-');
  if (hyphen != std::string::npos) {
    std::string swapped = target.substr(hyphen + 1) + "-" + target.substr(0, hyphen);
    if (hist.partials.contains(swapped)) {
      return swapped;
    }
  }
  if (hist.partials.contains("Total")) {
    return "Total";
  }
  return target;
}

void syncBondCutoffsRow(const AppWindow &window, int pair_index, const slint::SharedString &r_min,
                        const slint::SharedString &r_max) {
  const auto slint_cutoffs = window.get_bond_cutoffs();
  if (slint_cutoffs == nullptr || static_cast<size_t>(pair_index) >= slint_cutoffs->row_count()) {
    return;
  }
  const auto maybe_row = slint_cutoffs->row_data(pair_index);
  if (maybe_row.has_value()) {
    auto row = maybe_row.value();
    row.min_distance = r_min;
    row.max_distance = r_max;
    slint_cutoffs->set_row_data(pair_index, row);
  }
}

[[nodiscard]] std::optional<real_t> parseDistanceValue(std::string_view text) {
  if (text.empty() || text == "—") {
    return std::nullopt;
  }
  try {
    return static_cast<real_t>(std::stod(std::string(text)));
  } catch (const std::exception &) {
    return std::nullopt;
  }
}

} // namespace

void BondCutoffController::applyCovalentFactor(float factor, FactorBound bound) {
  if (loader_.cell() == nullptr || factor <= 0.0F) {
    return;
  }
  const auto &elements = loader_.cell()->elements();
  const auto num_elements = elements.size();
  const auto slint_cutoffs = window_.get_bond_cutoffs();

  const auto new_cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();
  size_t row_index = 0;
  for (size_t i = 0; i < num_elements; ++i) {
    const real_t radius_a = safeGetCovalentRadius(elements[i].symbol);
    for (size_t j = i; j < num_elements; ++j) {
      const real_t radius_b = safeGetCovalentRadius(elements[j].symbol);
      new_cutoffs->push_back(createCovalentCutoffRow(
          CovalentCutoffInput{
              .symbol_a = elements[i].symbol,
              .symbol_b = elements[j].symbol,
              .radius_a = radius_a,
              .radius_b = radius_b,
              .row_idx = row_index,
              .factor = factor,
              .bound = bound,
          },
          slint_cutoffs));
      ++row_index;
    }
  }

  window_.set_bond_cutoffs(new_cutoffs);
  window_.set_bond_cutoffs_reset_trigger(window_.get_bond_cutoffs_reset_trigger() + 1);
  options_.bond_cutoffs = getBondCutoffs();
  BondCutoffService::setBondCutoffs(loader_.trajectoryMut(), options_.bond_cutoffs);
  syncCutoffsToUI();
}

void BondCutoffController::applyGlobalCutoff(float global_cutoff) {
  setUniformCutoff(global_cutoff);
}

void BondCutoffController::initializeCalibrationPairs() {
  if (loader_.cell() == nullptr) {
    return;
  }

  const auto &elements = loader_.cell()->elements();
  const auto entries = BondCutoffMapper::createDefaultCutoffEntries(elements);

  const auto calib_pairs = std::make_shared<slint::VectorModel<CalibrationPair>>();
  for (size_t idx = 0; idx < entries.size(); ++idx) {
    const auto &entry = entries[idx];
    calib_pairs->push_back(CalibrationPair{
        .pair_index = static_cast<int>(idx),
        .label = slint::SharedString(std::format("{}-{}", entry.element1, entry.element2)),
        .r_min = slint::SharedString(entry.min_distance),
        .r_max = slint::SharedString(entry.max_distance),
        .r_peak = slint::SharedString("—"),
        .r_min_px = 0.0F,
        .r_max_px = 0.0F,
        .r_peak_px = 0.0F,
        .is_valid = true,
    });
  }

  window_.set_calibration_pairs(calib_pairs);
  updateCalibrationPixels(last_canvas_w_, last_canvas_h_);
  evaluateSanityState();
}

void BondCutoffController::setCutoffRange(int pair_index, float r_min, float r_max) {
  if (loader_.cell() == nullptr || pair_index < 0) {
    return;
  }

  const auto calib_pairs = window_.get_calibration_pairs();
  if (calib_pairs != nullptr && static_cast<size_t>(pair_index) < calib_pairs->row_count()) {
    const auto maybe_pair = calib_pairs->row_data(pair_index);
    if (maybe_pair.has_value()) {
      auto pair = maybe_pair.value();
      pair.r_min = slint::SharedString(std::format("{:.2f}", r_min));
      pair.r_max = slint::SharedString(std::format("{:.2f}", r_max));
      calib_pairs->set_row_data(pair_index, pair);
    }
  }

  syncBondCutoffsRow(window_, pair_index, slint::SharedString(std::format("{:.2f}", r_min)),
                     slint::SharedString(std::format("{:.2f}", r_max)));

  options_.bond_cutoffs = getBondCutoffs();
  BondCutoffService::setBondCutoffs(loader_.trajectoryMut(), options_.bond_cutoffs);
  updateCalibrationPixels(last_canvas_w_, last_canvas_h_);
  evaluateSanityState();
}

void BondCutoffController::resetCalibrationPair(int pair_index) {
  if (loader_.cell() == nullptr) {
    return;
  }
  const auto &elements = loader_.cell()->elements();
  const auto entries = BondCutoffMapper::createDefaultCutoffEntries(elements);
  if (pair_index < 0 || static_cast<size_t>(pair_index) >= entries.size()) {
    return;
  }

  const auto &entry = entries[pair_index];
  const auto calib_pairs = window_.get_calibration_pairs();
  if (calib_pairs != nullptr && static_cast<size_t>(pair_index) < calib_pairs->row_count()) {
    const auto maybe_pair = calib_pairs->row_data(pair_index);
    if (maybe_pair.has_value()) {
      auto pair = maybe_pair.value();
      pair.r_min = slint::SharedString(entry.min_distance);
      pair.r_max = slint::SharedString(entry.max_distance);
      pair.r_peak = slint::SharedString("—");
      calib_pairs->set_row_data(pair_index, pair);
    }
  }

  syncBondCutoffsRow(window_, pair_index, slint::SharedString(entry.min_distance),
                     slint::SharedString(entry.max_distance));

  options_.bond_cutoffs = getBondCutoffs();
  BondCutoffService::setBondCutoffs(loader_.trajectoryMut(), options_.bond_cutoffs);
  updateCalibrationPixels(last_canvas_w_, last_canvas_h_);
  evaluateSanityState();
}

void BondCutoffController::autoDetectCutoffs(int pair_index) {
  if (loader_.cell() == nullptr || dispatcher_ == nullptr) {
    return;
  }

  const auto calib_pairs = window_.get_calibration_pairs();
  if (calib_pairs == nullptr || pair_index < 0 ||
      static_cast<size_t>(pair_index) >= calib_pairs->row_count()) {
    return;
  }

  const auto maybe_pair = calib_pairs->row_data(pair_index);
  if (!maybe_pair.has_value()) {
    return;
  }
  auto pair = maybe_pair.value();

  const auto *hist = findRdfHistogram(dispatcher_);
  if (hist == nullptr || hist->bins.empty()) {
    return;
  }

  const std::string target_partial = resolveTargetPartial(*hist, pair.label.data());
  correlation::analysis::DetectionOptions opts;
  opts.min_search_radius = static_cast<real_t>(1.0);
  opts.max_search_radius = static_cast<real_t>(6.0);

  auto result =
      correlation::analysis::FirstMinimumDetector::detectCutoffs(*hist, target_partial, opts);
  if (result.has_value()) {
    const auto &cutoff = result.value();
    pair.r_min = slint::SharedString(std::format("{:.2f}", cutoff.r_min));
    pair.r_max = slint::SharedString(std::format("{:.2f}", cutoff.r_max));
    pair.r_peak = slint::SharedString(std::format("{:.2f}", cutoff.r_peak));
    calib_pairs->set_row_data(pair_index, pair);

    syncBondCutoffsRow(window_, pair_index, pair.r_min, pair.r_max);

    options_.bond_cutoffs = getBondCutoffs();
    BondCutoffService::setBondCutoffs(loader_.trajectoryMut(), options_.bond_cutoffs);
    updateCalibrationPixels(last_canvas_w_, last_canvas_h_);
    evaluateSanityState();
  }
}

void BondCutoffController::dragCutoffHandle(int pair_index, bool is_max, float mouse_x) {
  const auto calib_pairs = window_.get_calibration_pairs();
  if (calib_pairs == nullptr || pair_index < 0 ||
      static_cast<size_t>(pair_index) >= calib_pairs->row_count()) {
    return;
  }

  const auto maybe_pair = calib_pairs->row_data(pair_index);
  if (!maybe_pair.has_value()) {
    return;
  }
  const auto &pair = maybe_pair.value();

  float canvas_w = last_canvas_w_;
  float canvas_h = last_canvas_h_;
  if (plot_controller_ != nullptr) {
    if (plot_controller_->getLastPlotWidth() > 10.0F) {
      canvas_w = plot_controller_->getLastPlotWidth();
    }
    if (plot_controller_->getLastPlotHeight() > 10.0F) {
      canvas_h = plot_controller_->getLastPlotHeight();
    }
  }

  correlation::analysis::Histogram fallback_hist;
  fallback_hist.bins = {0.0, 15.0};
  const correlation::analysis::Histogram *hist = &fallback_hist;
  if (dispatcher_ != nullptr) {
    const auto *found = dispatcher_->getHistogram("Radial Distribution Function");
    if (found != nullptr && !found->bins.empty()) {
      hist = found;
    }
  }

  correlation::plotters::PlotConfig config;
  if (plot_controller_ != nullptr) {
    config = plot_controller_->getActivePlotConfig();
  }

  const auto [data_r, unused_y] = PlotViewportNavigator::screenToData(
      mouse_x, canvas_h / 2.0F, canvas_w, canvas_h, *hist, config);

  const real_t current_min =
      parseDistanceValue(pair.r_min.data()).value_or(static_cast<real_t>(0.5));
  const real_t current_max =
      parseDistanceValue(pair.r_max.data()).value_or(static_cast<real_t>(3.0));

  if (is_max) {
    const real_t new_max =
        std::clamp(data_r, current_min + static_cast<real_t>(0.05), static_cast<real_t>(25.0));
    setCutoffRange(pair_index, static_cast<float>(current_min), static_cast<float>(new_max));
  } else {
    const real_t new_min =
        std::clamp(data_r, static_cast<real_t>(0.1), current_max - static_cast<real_t>(0.05));
    setCutoffRange(pair_index, static_cast<float>(new_min), static_cast<float>(current_max));
  }
}

void BondCutoffController::updateCalibrationPixels(float canvas_w, float canvas_h) {
  if (canvas_w > 10.0F) {
    last_canvas_w_ = canvas_w;
  }
  if (canvas_h > 10.0F) {
    last_canvas_h_ = canvas_h;
  }

  const auto calib_pairs = window_.get_calibration_pairs();
  if (calib_pairs == nullptr || calib_pairs->row_count() == 0) {
    return;
  }

  correlation::analysis::Histogram fallback_hist;
  fallback_hist.bins = {0.0, 15.0};
  const correlation::analysis::Histogram *hist = &fallback_hist;
  if (dispatcher_ != nullptr) {
    const auto *found = dispatcher_->getHistogram("Radial Distribution Function");
    if (found != nullptr && !found->bins.empty()) {
      hist = found;
    }
  }

  correlation::plotters::PlotConfig config;
  if (plot_controller_ != nullptr) {
    config = plot_controller_->getActivePlotConfig();
  }

  for (size_t k = 0; k < calib_pairs->row_count(); ++k) {
    const auto maybe_pair = calib_pairs->row_data(k);
    if (!maybe_pair.has_value()) {
      continue;
    }
    auto pair = maybe_pair.value();

    if (const auto min_val = parseDistanceValue(pair.r_min.data()); min_val.has_value()) {
      const auto [screen_x, unused_y] = PlotViewportNavigator::dataToScreen(
          *min_val, 0.0, last_canvas_w_, last_canvas_h_, *hist, config);
      pair.r_min_px = static_cast<float>(screen_x);
    } else {
      pair.r_min_px = 0.0F;
    }

    if (const auto max_val = parseDistanceValue(pair.r_max.data()); max_val.has_value()) {
      const auto [screen_x, unused_y] = PlotViewportNavigator::dataToScreen(
          *max_val, 0.0, last_canvas_w_, last_canvas_h_, *hist, config);
      pair.r_max_px = static_cast<float>(screen_x);
    } else {
      pair.r_max_px = 0.0F;
    }

    if (const auto peak_val = parseDistanceValue(pair.r_peak.data()); peak_val.has_value()) {
      const auto [screen_x, unused_y] = PlotViewportNavigator::dataToScreen(
          *peak_val, 0.0, last_canvas_w_, last_canvas_h_, *hist, config);
      pair.r_peak_px = static_cast<float>(screen_x);
    } else {
      pair.r_peak_px = -1.0F;
    }

    calib_pairs->set_row_data(k, pair);
  }
}

void BondCutoffController::syncCutoffsToUI() {
  if (loader_.cell() == nullptr) {
    evaluateSanityState();
    return;
  }

  const auto &elements = loader_.cell()->elements();
  const auto num_elements = elements.size();
  const auto slint_cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();
  const auto calib_pairs = std::make_shared<slint::VectorModel<CalibrationPair>>();

  int pair_idx = 0;
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      const real_t min_d = std::sqrt(options_.bond_cutoffs[i][j].min_sq);
      const real_t max_d = std::sqrt(options_.bond_cutoffs[i][j].max_sq);
      const auto min_str = std::format("{:.2f}", min_d);
      const auto max_str = std::format("{:.2f}", max_d);

      slint_cutoffs->push_back(BondCutoff{
          .element1 = slint::SharedString(elements[i].symbol),
          .element2 = slint::SharedString(elements[j].symbol),
          .min_distance = slint::SharedString(min_str),
          .max_distance = slint::SharedString(max_str),
      });

      calib_pairs->push_back(CalibrationPair{
          .pair_index = pair_idx,
          .label =
              slint::SharedString(std::format("{}-{}", elements[i].symbol, elements[j].symbol)),
          .r_min = slint::SharedString(min_str),
          .r_max = slint::SharedString(max_str),
          .r_peak = slint::SharedString("—"),
          .r_min_px = 0.0F,
          .r_max_px = 0.0F,
          .r_peak_px = 0.0F,
          .is_valid = true,
      });
      ++pair_idx;
    }
  }

  window_.set_bond_cutoffs(slint_cutoffs);
  window_.set_bond_cutoffs_reset_trigger(window_.get_bond_cutoffs_reset_trigger() + 1);
  window_.set_calibration_pairs(calib_pairs);
  updateCalibrationPixels(last_canvas_w_, last_canvas_h_);
  evaluateSanityState();
}

void BondCutoffController::evaluateSanityState() {
  if (loader_.cell() == nullptr) {
    window_.set_sanity_badge_state({
        .severity = 0,
        .title = slint::SharedString("No Structure"),
        .message = slint::SharedString("Load a structure or trajectory to validate bond cutoffs."),
        .has_fix = false,
        .fix_action = 0,
        .fix_label = slint::SharedString(""),
    });
    return;
  }

  const auto &elements = loader_.cell()->elements();
  const auto report = correlation::analysis::CutoffSanityChecker::checkCutoffs(
      options_.bond_cutoffs, elements, loader_.cell());

  window_.set_sanity_badge_state({
      .severity = static_cast<int>(report.severity),
      .title = slint::SharedString(report.title),
      .message = slint::SharedString(report.message),
      .has_fix = report.has_fix,
      .fix_action = static_cast<int>(report.fix_action),
      .fix_label = slint::SharedString(report.fix_label),
  });
}

void BondCutoffController::handleFixSanityIssue(int fix_action) {
  if (loader_.cell() == nullptr) {
    return;
  }

  const auto action = static_cast<correlation::analysis::SanityFixAction>(fix_action);
  if (action == correlation::analysis::SanityFixAction::AutoDetectFirstMin) {
    int active_idx = window_.get_active_calibration_pair_index();
    const auto calib_pairs = window_.get_calibration_pairs();
    if (calib_pairs != nullptr && calib_pairs->row_count() > 0) {
      if (active_idx < 0 || static_cast<size_t>(active_idx) >= calib_pairs->row_count()) {
        active_idx = 0;
      }
      autoDetectCutoffs(active_idx);
    }
    return;
  }

  correlation::analysis::CutoffSanityChecker::applyFix(
      options_.bond_cutoffs, loader_.cell()->elements(), loader_.cell(), action);
  BondCutoffService::setBondCutoffs(loader_.trajectoryMut(), options_.bond_cutoffs);

  syncCutoffsToUI();
}

} // namespace correlation::app
