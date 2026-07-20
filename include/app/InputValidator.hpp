/**
 * @file InputValidator.hpp
 * @brief Input validation and CLI command generation logic.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "app/AppBackend.hpp"

class AppWindow;

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
   * @param[in,out] window Reference to the UI window.
   * @param[in,out] backend Reference to the application backend.
   * @param[in,out] controller Reference to the main AppController.
   */
  InputValidator(::AppWindow &window, AppBackend &backend, AppController &controller);

  /**
   * @brief Validates all numeric input fields and pushes error states to the UI.
   * @return true if all inputs are valid, false otherwise.
   */
  [[nodiscard]] bool validateInputs();

  /**
   * @brief Generates and updates the equivalent CLI command in the UI.
   */
  void updateCliCommand();

  /**
   * @brief Copies the equivalent CLI command to the system clipboard.
   */
  void handleCopyCliCommand();

private:
  /**
   * @brief Validates radial distribution and scattering options.
   * @param[out] errs AppErrors structure for reporting invalid field states.
   * @param[out] r_max_val Evaluated maximum radial cutoff.
   * @param[out] q_max_val Evaluated maximum reciprocal space momentum.
   * @return true if valid, false otherwise.
   */
  [[nodiscard]] bool validateRadialAndScattering(AppErrors &errs, float &r_max_val, float &q_max_val);

  /**
   * @brief Validates angular and ring distribution options.
   * @param[out] errs AppErrors structure for reporting invalid field states.
   * @return true if valid, false otherwise.
   */
  [[nodiscard]] bool validateAngularAndRings(AppErrors &errs);

  /**
   * @brief Validates Local Entropy and Hyperuniformity parameters.
   * @param[out] errs AppErrors structure for reporting invalid field states.
   * @return true if valid, false otherwise.
   */
  [[nodiscard]] bool validateOtherAnalysisOptions(AppErrors &errs);

  /**
   * @brief Validates trajectory frame indexing bounds.
   * @param[out] errs AppErrors structure for reporting invalid field states.
   * @return true if valid, false otherwise.
   */
  [[nodiscard]] bool validateFrames(AppErrors &errs);

  /**
   * @brief Validates plot export layout settings.
   * @param[out] errs AppErrors structure for reporting invalid field states.
   * @return true if valid, false otherwise.
   */
  [[nodiscard]] bool validateExportConfig(AppErrors &errs);

  ::AppWindow *window_;
  AppBackend *backend_;
  AppController *controller_;
};

} // namespace correlation::app
