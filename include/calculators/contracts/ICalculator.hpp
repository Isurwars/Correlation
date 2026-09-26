/**
 * @file ICalculator.hpp
 * @brief Base contract for all structural analysis calculators.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <string_view>

namespace correlation::calculators {

/**
 * @brief Primary interface defining metadata and capability queries for all calculators.
 */
class ICalculator {
public:
  ICalculator() = default;
  virtual ~ICalculator() = default;
  ICalculator(const ICalculator &) = delete;
  ICalculator &operator=(const ICalculator &) = delete;
  ICalculator(ICalculator &&) = delete;
  ICalculator &operator=(ICalculator &&) = delete;

  /**
   * @brief Returns the unique identifier/name of the calculator.
   * @return The full name of the calculator (e.g. "Radial Distribution Function").
   */
  [[nodiscard]] virtual std::string_view getName() const = 0;

  /**
   * @brief Returns a short, UI-friendly name of the calculator (e.g. "g_r", "S_q").
   * @return The short name/abbreviation of the calculator.
   */
  [[nodiscard]] virtual std::string_view getShortName() const = 0;

  /**
   * @brief Returns the UI group this calculator belongs to (e.g., "Radial", "Angular", "Dynamic").
   * @return The group name for UI categorization.
   */
  [[nodiscard]] virtual std::string_view getGroup() const = 0;

  /**
   * @brief Returns a brief description of the calculator's purpose.
   * @return A human-readable description string.
   */
  [[nodiscard]] virtual std::string_view getDescription() const = 0;

  /**
   * @brief Check if this calculator runs per-frame (e.g. RDF, PAD).
   * @return True if it calculates properties for individual snapshots.
   */
  [[nodiscard]] virtual bool isFrameCalculator() const = 0;

  /**
   * @brief Check if this calculator runs on the whole trajectory (e.g. VACF, MSD).
   * @return True if it calculates time-dependent or multi-frame properties.
   */
  [[nodiscard]] virtual bool isTrajectoryCalculator() const = 0;

  /**
   * @brief Check if the calculator has all required models, parameters, and resources configured.
   * @return True if configured and runnable; false if prerequisites are missing.
   */
  [[nodiscard]] virtual bool isConfigured() const { return true; }
};

} // namespace correlation::calculators
