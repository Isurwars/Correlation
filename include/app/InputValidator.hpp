/**
 * @file InputValidator.hpp
 * @brief Input validation and CLI command generation logic.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "AppWindow.h"
#include "app/AppBackend.hpp"

namespace correlation::app {

class AppController; // Forward declaration

/**
 * @class InputValidator
 * @brief Handles UI input validation and equivalent CLI command generation.
 */
class InputValidator {
public:
  /**
   * @brief Constructs the InputValidator.
   * @param window Reference to the UI window.
   * @param backend Reference to the application backend.
   * @param controller Reference to the main AppController.
   */
  InputValidator(AppWindow &window, AppBackend &backend, AppController &controller);

  /**
   * @brief Validates all numeric input fields and pushes error states to the UI.
   * @return true if all inputs are valid, false otherwise.
   */
  bool validateInputs();

  /**
   * @brief Generates and updates the equivalent CLI command in the UI.
   */
  void updateCliCommand();

  /**
   * @brief Copies the equivalent CLI command to the clipboard.
   */
  void handleCopyCliCommand();

private:
  AppWindow *window_;
  AppBackend *backend_;
  AppController *controller_;
};

} // namespace correlation::app
