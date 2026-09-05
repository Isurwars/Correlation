/**
 * @file AnalysisRunner.cpp
 * @brief Implementation of AnalysisRunner.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/AnalysisRunner.hpp"
#include "AppWindow.h"
#include "app/AppController.hpp"
#include "app/InputValidator.hpp"
#include "app/PlotController.hpp"
#include <algorithm>

namespace correlation::app {

AnalysisRunner::AnalysisRunner(::AppWindow &window, AppBackend &backend, AppController &controller)
    : window_(window), backend_(backend), controller_(controller) {}

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
  backend_.setOptions(controller_.handleOptionsfromUI());

  // Set the progress callback
  backend_.setProgressCallback([this](float progress, const std::string &msg) { updateProgress(progress, msg); });

  // Run analysis in a separate thread asynchronously without blocking GUI event loop
  if (analysis_thread_.joinable()) {
    std::thread old_thread = std::move(analysis_thread_);
    std::thread([t = std::move(old_thread)]() mutable {
      if (t.joinable()) {
        t.join();
      }
    }).detach();
  }

  analysis_thread_ = std::thread([this]() {
    auto result = backend_.run_analysis();
    std::string err = result ? "" : result.error();

    slint::invoke_from_event_loop([this, err]() {
      window_.set_analysis_running(false);
      if (err.empty()) {
        if (backend_.is_cancelled()) {
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
      if (!backend_.getAvailableHistogramNames().empty()) {
        controller_.getPlotController()->requestPlotUpdate(0, true);
      }
    });
  });

  // Detach or move is not enough, we need to keep the thread object alive.
  // We keep it as a member variable.
}

} // namespace correlation::app
