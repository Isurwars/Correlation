/**
 * @file InputValidator.hpp
 * @brief Input validation logic.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "AppWindow.h"

class AppWindow;

namespace correlation::app {

class AppController; // Forward declaration

/**
 * @class InputValidator
 * @brief Handles UI input validation.
 */
class InputValidator {
public:
  /**
   * @brief Constructs the InputValidator.
   * @param[in,out] window Reference to the UI window.
   * @param[in,out] controller Reference to the main AppController.
   */
  InputValidator(::AppWindow &window, AppController &controller);

  /**
   * @brief Validates all numeric input fields and pushes error states to the UI.
   * @return true if all inputs are valid, false otherwise.
   */
  [[nodiscard]] bool validateInputs();

private:
  /**
   * @brief Validates radial distribution and scattering options.
   * @param[in] opts Analysis options from the UI.
   * @param[out] errs AppErrors structure for reporting invalid field states.
   * @param[out] r_max_val Evaluated maximum radial cutoff.
   * @param[out] q_max_val Evaluated maximum reciprocal space momentum.
   * @return true if valid, false otherwise.
   */
  [[nodiscard]] static bool validateRadialAndScattering(const AnalysisOptions &opts,
                                                        AppErrors &errs, float &r_max_val,
                                                        float &q_max_val);

  /**
   * @brief Validates powder X-ray diffraction options.
   * @param[in] opts Analysis options from the UI.
   * @param[out] errs AppErrors structure for reporting invalid field states.
   * @return true if valid, false otherwise.
   */
  [[nodiscard]] static bool validateXrdOptions(const AnalysisOptions &opts, AppErrors &errs);

  /**
   * @brief Validates angular and ring distribution options.
   * @param[in] opts Analysis options from the UI.
   * @param[out] errs AppErrors structure for reporting invalid field states.
   * @return true if valid, false otherwise.
   */
  [[nodiscard]] static bool validateAngularAndRings(const AnalysisOptions &opts, AppErrors &errs);

  /**
   * @brief Validates Local Entropy and Hyperuniformity parameters.
   * @param[in] opts Analysis options from the UI.
   * @param[out] errs AppErrors structure for reporting invalid field states.
   * @return true if valid, false otherwise.
   */
  [[nodiscard]] static bool validateOtherAnalysisOptions(const AnalysisOptions &opts,
                                                         AppErrors &errs);

  /**
   * @brief Validates trajectory frame indexing bounds.
   * @param[in] opts Analysis options from the UI.
   * @param[out] errs AppErrors structure for reporting invalid field states.
   * @return true if valid, false otherwise.
   */
  [[nodiscard]] bool validateFrames(const AnalysisOptions &opts, AppErrors &errs);

  /**
   * @brief Validates plot export layout settings.
   * @param[out] errs AppErrors structure for reporting invalid field states.
   * @return true if valid, false otherwise.
   */
  [[nodiscard]] bool validateExportConfig(AppErrors &errs);

  ::AppWindow *window_;
  AppController *controller_;
};

} // namespace correlation::app
