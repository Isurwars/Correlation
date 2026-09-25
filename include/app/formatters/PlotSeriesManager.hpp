/**
 * @file PlotSeriesManager.hpp
 * @brief Manages curve visibility, custom colors, pinned runs, and toggle models for plots.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "AppWindow.h"
#include "analysis/DistributionFunctions.hpp"
#include "plotters/PlotTypes.hpp"
#include <slint.h>

#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace correlation::app {

/**
 * @struct PinnedRun
 * @brief Represents a snapshot of analytical histograms preserved for comparison overlay.
 */
struct PinnedRun {
  std::string label; ///< Identifier label for the run (e.g., "Run 1").
  std::map<std::string, correlation::analysis::Histogram>
      histograms; ///< Map of plot name to histogram.
};

/**
 * @class PlotSeriesManager
 * @brief Manages curve visibility, custom palette colors, pinned comparison runs, and Slint toggle
 * models.
 */
class PlotSeriesManager {
public:
  /**
   * @brief Default constructor.
   */
  PlotSeriesManager() = default;

  /**
   * @brief Toggles the visibility of a curve by its toggle index.
   * @param[in] curve_id Index corresponding to the current toggle item.
   * @param[in] visible Target visibility state.
   */
  void setCurveVisible(int curve_id, bool visible);

  /**
   * @brief Sets the visibility for all currently tracked toggle curves.
   * @param[in] visible Target visibility state for all curves.
   */
  void setAllCurvesVisible(bool visible);

  /**
   * @brief Checks if a specific curve key is visible.
   * @param[in] key Curve key name.
   * @param[in] default_val Default state returned if key is not explicitly registered.
   * @return True if curve is visible, false otherwise.
   */
  [[nodiscard]] bool isCurveVisible(const std::string &key, bool default_val = true) const;

  /**
   * @brief Retrieves the complete curve visibility map.
   * @return Const reference to the visibility map.
   */
  [[nodiscard]] const std::map<std::string, bool> &getCurveVisibilityMap() const noexcept {
    return curve_visibility_map_;
  }

  /**
   * @brief Sets a custom hex color for a curve by its toggle index.
   * @param[in] curve_id Index corresponding to the current toggle item.
   * @param[in] color_hex Hex color string (e.g. "#RRGGBB").
   */
  void setCustomColor(int curve_id, const std::string &color_hex);

  /**
   * @brief Retrieves the map of custom curve colors.
   * @return Const reference to custom colors map.
   */
  [[nodiscard]] const std::map<std::string, std::string> &getCustomColors() const noexcept {
    return custom_curve_colors_;
  }

  /**
   * @brief Pins the provided histogram dataset as a comparison overlay run.
   * @param[in] hists Map of plot names to histograms.
   */
  void pinCurrentRun(const std::map<std::string, correlation::analysis::Histogram> &hists);

  /**
   * @brief Clears all pinned comparison runs.
   */
  void clearPinnedRuns() noexcept;

  /**
   * @brief Retrieves all pinned comparison runs.
   * @return Const reference to vector of pinned runs.
   */
  [[nodiscard]] const std::vector<PinnedRun> &getPinnedRuns() const noexcept {
    return pinned_runs_;
  }

  /**
   * @brief Retrieves the count of pinned comparison runs.
   * @return Total pinned runs.
   */
  [[nodiscard]] size_t getPinnedRunsCount() const noexcept { return pinned_runs_.size(); }

  /**
   * @brief Sets whether the overlaid difference curve Y_diff is shown.
   * @param[in] show True to show difference curve, false to hide.
   */
  void setShowDifference(bool show) noexcept { show_difference_curve_ = show; }

  /**
   * @brief Queries whether the overlaid difference curve is visible.
   * @return True if difference curve is enabled.
   */
  [[nodiscard]] bool shouldShowDifference() const noexcept { return show_difference_curve_; }

  /**
   * @brief Retrieves the currently ordered list of curve toggle keys.
   * @return Const reference to current toggle keys vector.
   */
  [[nodiscard]] const std::vector<std::string> &getCurrentToggleKeys() const noexcept {
    return current_toggle_keys_;
  }

  /**
   * @brief Generates the Slint VectorModel of CurveToggleItem for UI presentation.
   * @param[in] hist Pointer to the current active analytical histogram.
   * @param[in] config Active plot rendering configuration for default palette colors.
   * @return Shared pointer to Slint VectorModel of CurveToggleItem.
   */
  [[nodiscard]] std::shared_ptr<slint::VectorModel<CurveToggleItem>>
  generateToggleItems(const correlation::analysis::Histogram *hist,
                      const correlation::plotters::PlotConfig &config);

  /**
   * @brief Resets all series state (visibility overrides, custom colors, pinned runs, and flags).
   */
  void reset() noexcept;

private:
  std::map<std::string, bool> curve_visibility_map_;
  std::map<std::string, std::string> custom_curve_colors_;
  std::vector<std::string> current_toggle_keys_;
  std::vector<PinnedRun> pinned_runs_;
  bool show_difference_curve_{false};
};

} // namespace correlation::app
