/**
 * @file AnalysisRunner.hpp
 * @brief Analysis execution logic.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "app/core/AppOptions.hpp"
#include "app/services/TrajectoryLoader.hpp"
#include <string>
#include <thread>

class AppWindow;

namespace correlation::app {

class AnalysisDispatcher;
class AppController;

/**
 * @class AnalysisRunner
 * @brief Manages the background thread and logic for running analysis workflows.
 */
class AnalysisRunner {
public:
  /**
   * @brief Constructs the AnalysisRunner with injected services.
   * @param[in,out] window Reference to the UI window.
   * @param[in,out] loader Reference to the trajectory loader.
   * @param[in,out] dispatcher Reference to the analysis dispatcher.
   * @param[in,out] options Reference to program options.
   * @param[in,out] controller Reference to the main AppController.
   */
  AnalysisRunner(::AppWindow &window, TrajectoryLoader &loader, AnalysisDispatcher &dispatcher,
                 ProgramOptions &options, AppController &controller);

  /**
   * @brief Destructor. Ensures analysis thread is joined.
   */
  ~AnalysisRunner();

  /**
   * @brief Handles triggering the trajectory analysis execution.
   */
  void handleRunAnalysis();

  /**
   * @brief Updates the UI progress indicator and status text.
   * @param[in] progress Completion fraction [0.0 - 1.0].
   * @param[in] msg Status update message string.
   */
  void updateProgress(float progress, const std::string &msg);

private:
  ::AppWindow &window_;
  TrajectoryLoader &loader_;
  AnalysisDispatcher &dispatcher_;
  ProgramOptions &options_;
  AppController &controller_;

  std::jthread analysis_thread_; ///< Handle for the background analysis computation.

  /**
   * @brief Asynchronously joins any running background thread to prevent blocking the UI.
   */
  void joinPreviousThreadAsync();

  /**
   * @brief Executes the analysis workflow synchronously within the worker thread.
   * @return Error message if analysis failed or was aborted, empty string otherwise.
   */
  [[nodiscard]] std::string executeAnalysis();

  /**
   * @brief Updates the UI state and models once the analysis has completed.
   * @param[in] err Error message from execution, if any.
   */
  void handleAnalysisCompletion(const std::string &err);
};

} // namespace correlation::app
