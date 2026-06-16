/**
 * @file AnalysisRunner.hpp
 * @brief Analysis execution logic.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "AppWindow.h"
#include "app/AppBackend.hpp"
#include <thread>
#include <string>

namespace correlation::app {

class AppController;

/**
 * @class AnalysisRunner
 * @brief Manages the background thread and logic for running analysis.
 */
class AnalysisRunner {
public:
  AnalysisRunner(AppWindow &window, AppBackend &backend, AppController &controller);
  ~AnalysisRunner();

  void handleRunAnalysis();
  void updateProgress(float progress, const std::string &msg);

private:
  AppWindow &window_;
  AppBackend &backend_;
  AppController &controller_;

  std::thread analysis_thread_; ///< Handle for the background analysis computation.
};

} // namespace correlation::app
