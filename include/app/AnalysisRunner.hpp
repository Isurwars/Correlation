/**
 * @file AnalysisRunner.hpp
 * @brief Analysis execution logic.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "app/AppBackend.hpp"
#include <string>
#include <thread>

class AppWindow;

namespace correlation::app {

class AppController;

/**
 * @class AnalysisRunner
 * @brief Manages the background thread and logic for running analysis workflows.
 */
class AnalysisRunner {
public:
  /**
   * @brief Constructs the AnalysisRunner.
   * @param[in,out] window Reference to the UI window.
   * @param[in,out] backend Reference to the application backend.
   * @param[in,out] controller Reference to the main AppController.
   */
  AnalysisRunner(::AppWindow &window, AppBackend &backend, AppController &controller);

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
  AppBackend &backend_;
  AppController &controller_;

  std::thread analysis_thread_; ///< Handle for the background analysis computation.
};

} // namespace correlation::app
