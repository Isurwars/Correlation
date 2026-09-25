/**
 * @file AnalysisRunner.cpp
 * @brief Implementation of AnalysisRunner.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/AnalysisRunner.hpp"
#include "AppWindow.h"
#include "app/AppBackend.hpp"
#include "app/AppController.hpp"
#include "app/InputValidator.hpp"
#include "app/PlotController.hpp"
#include <algorithm>

namespace correlation::app {

AnalysisRunner::AnalysisRunner(::AppWindow &window, TrajectoryLoader &loader,
                               AnalysisDispatcher &dispatcher, ProgramOptions &options,
                               AppController &controller)
    : window_(window), loader_(loader), dispatcher_(dispatcher), options_(options),
      controller_(controller) {}

AnalysisRunner::AnalysisRunner(::AppWindow &window, AppBackend &backend, AppController &controller)
    : AnalysisRunner(window, backend.loader(), backend.dispatcher(), backend.options(),
                     controller) {}

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
  if (analysis_thread_.joinable()) {
    std::thread old_thread = std::move(analysis_thread_);
    std::thread([thread_to_join = std::move(old_thread)]() mutable {
      if (thread_to_join.joinable()) {
        thread_to_join.join();
      }
    }).detach();
  }

  analysis_thread_ = std::thread([this]() {
    std::string err;
    if (loader_.trajectory() == nullptr || loader_.getFrameCount() == 0) {
      err = AppDefaults::MSG_ANALYSIS_ABORTED;
    } else {
      auto result = dispatcher_.runAnalysis(*loader_.trajectoryMut(), options_);
      if (!result) {
        err = result.error();
      }
    }

    slint::invoke_from_event_loop([this, err]() {
      window_.set_analysis_running(false);
      if (err.empty()) {
        if (dispatcher_.isCancelled()) {
          window_.set_analysis_status_text(slint::SharedString("Analysis Cancelled."));
        } else {
          window_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_ANALYSIS_ENDED));
        }
      } else {
        window_.set_analysis_status_text(slint::SharedString(err));
      }
      window_.set_analysis_done(true);
      window_.set_progress(1.0F);

      // Populate the plot dropdown and auto-preview the first histogram
      controller_.getPlotController()->populatePlotList();
      if (!dispatcher_.getAvailableHistogramNames().empty()) {
        controller_.getPlotController()->requestPlotUpdate(0, true);
      }
    });
  });

  // Detach or move is not enough, we need to keep the thread object alive.
  // We keep it as a member variable.
}

} // namespace correlation::app
