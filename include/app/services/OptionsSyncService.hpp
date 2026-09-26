/**
 * @file OptionsSyncService.hpp
 * @brief Stateless service synchronizing ProgramOptions with Slint UI properties.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "AppWindow.h"
#include "app/core/AppOptions.hpp"

#include <cstddef>
#include <expected>
#include <string>

namespace correlation::app {

/**
 * @class OptionsSyncService
 * @brief Stateless bidirectional synchronization between ProgramOptions and UI analysis models.
 */
class OptionsSyncService {
public:
  /**
   * @brief Reads and validates analysis options from the UI window into a ProgramOptions struct.
   * @param[in] window Source UI window.
   * @param[in] frame_count Total available frames in loaded trajectory.
   * @param[in] bond_cutoffs Active bond cutoff matrix to attach.
   * @return Parsed and validated ProgramOptions on success, or an error string on failure.
   */
  [[nodiscard]] static std::expected<ProgramOptions, std::string>
  readFromUI(const AppWindow &window, size_t frame_count,
             const correlation::analysis::BondCutoffMatrix &bond_cutoffs = {});

  /**
   * @brief Writes ProgramOptions fields to the target UI window.
   * @param[in,out] window Target UI window.
   * @param[in] options Source program options to serialize.
   */
  static void writeToUI(AppWindow &window, const ProgramOptions &options);

  /**
   * @brief Updates boolean UI visibility/active flags for calculator groups.
   * @param[in,out] window Target UI window.
   * @param[in] opts Source program options containing active calculators.
   */
  static void updateActiveGroupFlags(AppWindow &window, const ProgramOptions &opts);
};

} // namespace correlation::app
