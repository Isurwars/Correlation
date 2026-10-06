/**
 * @file AnalysisRunner.cpp
 * @brief Implementation of AnalysisRunner.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/services/AnalysisRunner.hpp"
#include "AppWindow.h"
#include "app/core/AppController.hpp"
#include "app/services/AnalysisDispatcher.hpp"
#include "app/services/InputValidator.hpp"
#include "app/viewmodel/PlotController.hpp"
#include <algorithm>

namespace correlation::app {

AnalysisRunner::AnalysisRunner(::AppWindow &window, TrajectoryLoader &loader,
                               AnalysisDispatcher &dispatcher, ProgramOptions &options,
                               AppController &controller)
    : window_(window), loader_(loader), dispatcher_(dispatcher), options_(options),
      controller_(controller) {}

AnalysisRunner::~AnalysisRunner() {
  if (analysis_thread_.joinable()) {
    analysis_thread_.join();
  }
}

void AnalysisRunner::updateProgress(float progress, const std::string &msg) {
  progress = std::clamp(progress, 0.0F, 1.0F);
  slint::invoke_from_event_loop([progress, msg, this]() {
    window_.set_progress(progress);
    if (!msg.empty()) {
      window_.set_analysis_status_text(slint::SharedString(msg));
    }
  });
}

namespace {

[[nodiscard]] std::string determineStatusMessage(const std::string &err, bool is_cancelled) {
  if (!err.empty()) {
    return err;
  }
  if (is_cancelled) {
    return "Analysis Cancelled.";
  }
  return {AppDefaults::MSG_ANALYSIS_ENDED};
}

} // namespace

void AnalysisRunner::joinPreviousThreadAsync() {
  if (!analysis_thread_.joinable()) {
    return;
  }

  std::jthread old_thread = std::move(analysis_thread_);
  std::jthread([thread_to_join = std::move(old_thread)]() mutable {
    if (thread_to_join.joinable()) {
      thread_to_join.join();
    }
  }).detach();
}

std::string AnalysisRunner::executeAnalysis() {
  if (loader_.trajectory() == nullptr || loader_.getFrameCount() == 0) {
    return {AppDefaults::MSG_ANALYSIS_ABORTED};
  }

  const auto result = dispatcher_.runAnalysis(*loader_.trajectoryMut(), options_);
  if (!result) {
    return result.error();
  }
  return {};
}

void AnalysisRunner::handleAnalysisCompletion(const std::string &err) {
  const bool cancelled = dispatcher_.isCancelled();
  const bool successful = err.empty() && !cancelled;

  window_.set_analysis_running(false);
  window_.set_analysis_status_text(slint::SharedString(determineStatusMessage(err, cancelled)));
  window_.set_analysis_done(successful);
  window_.set_progress(successful ? 1.0F : 0.0F);

  if (!successful) {
    return;
  }

  controller_.getPlotController()->populatePlotList();
  controller_.populateExportAnalyses();
  if (!dispatcher_.getAvailableHistogramNames().empty()) {
    controller_.getPlotController()->requestPlotUpdate(0, true);
  }
}

void AnalysisRunner::handleRunAnalysis() {
  if (!controller_.getInputValidator()->validateInputs()) {
    return;
  }

  window_.set_analysis_done(false); // Reset done state
  window_.set_analysis_running(true);
  window_.set_progress(0.0F);
  window_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_RUNNING_ANALYSIS));

  // Create a ProgramOptions object from the UI properties
  options_ = controller_.handleOptionsfromUI();

  // Set the progress callback
  dispatcher_.setProgressCallback(
      [this](float progress, const std::string &msg) { updateProgress(progress, msg); });

  // Run analysis in a separate thread asynchronously without blocking GUI event loop
  joinPreviousThreadAsync();

  analysis_thread_ = std::jthread([this]() {
    const std::string err = executeAnalysis();
    slint::invoke_from_event_loop([this, err]() { handleAnalysisCompletion(err); });
  });
}

} // namespace correlation::app
