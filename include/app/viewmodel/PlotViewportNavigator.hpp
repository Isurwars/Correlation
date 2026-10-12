/**
 * @file PlotViewportNavigator.hpp
 * @brief Viewport navigation, coordinate projection, and zoom state management.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "plotters/SvgPlotter.hpp"

#include <optional>
#include <utility>

namespace correlation::app {

/**
 * @struct ViewportLimits
 * @brief Optional manual boundary limits representing active zoom window.
 */
struct ViewportLimits {
  std::optional<real_t> x_min; ///< Lower X bound in data units
  std::optional<real_t> x_max; ///< Upper X bound in data units
  std::optional<real_t> y_min; ///< Lower Y bound in data units
  std::optional<real_t> y_max; ///< Upper Y bound in data units
};

/**
 * @class PlotViewportNavigator
 * @brief Manages 2D plot coordinate space, active zoom window, and interactive pan/zoom
 * transitions.
 */
class PlotViewportNavigator {
public:
  /**
   * @brief Computes standard or manual nice axis scale for the X dimension.
   * @param[in] hist Active histogram data.
   * @param[in] config Base plotting configuration.
   * @return Formatted NiceScale with ticks and boundaries.
   */
  [[nodiscard]] static correlation::plotters::detail::NiceScale
  computeXScale(const correlation::analysis::Histogram &hist,
                const correlation::plotters::PlotConfig &config);

  /**
   * @brief Computes standard or manual nice axis scale for the Y dimension across all visible
   * partials.
   * @param[in] hist Active histogram data.
   * @param[in] config Base plotting configuration.
   * @return Formatted NiceScale with ticks and boundaries.
   */
  [[nodiscard]] static correlation::plotters::detail::NiceScale
  computeYScale(const correlation::analysis::Histogram &hist,
                const correlation::plotters::PlotConfig &config);

  /**
   * @brief Projects screen-space coordinates to data-space coordinates.
   * @param[in] screen_x Screen horizontal coordinate in pixels.
   * @param[in] screen_y Screen vertical coordinate in pixels.
   * @param[in] canvas_w Total canvas width in pixels.
   * @param[in] canvas_h Total canvas height in pixels.
   * @param[in] hist Active histogram data.
   * @param[in] config Base plotting configuration.
   * @return Pair containing (data_x, data_y).
   */
  [[nodiscard]] static std::pair<real_t, real_t>
  screenToData(float screen_x, float screen_y, float canvas_w, float canvas_h,
               const correlation::analysis::Histogram &hist,
               const correlation::plotters::PlotConfig &config);

  /**
   * @brief Checks if a custom zoom window is currently active.
   * @return True if any zoom coordinate boundary is set.
   */
  [[nodiscard]] bool hasActiveZoom() const noexcept;

  /**
   * @brief Applies rectangular selection box zoom.
   * @param[in] start_x Starting drag X coordinate.
   * @param[in] start_y Starting drag Y coordinate.
   * @param[in] end_x Ending drag X coordinate.
   * @param[in] end_y Ending drag Y coordinate.
   * @param[in] canvas_w Canvas width in pixels.
   * @param[in] canvas_h Canvas height in pixels.
   * @param[in] hist Active histogram data.
   * @param[in] config Base plotting configuration.
   */
  void applyZoomRect(float start_x, float start_y, float end_x, float end_y, float canvas_w,
                     float canvas_h, const correlation::analysis::Histogram &hist,
                     const correlation::plotters::PlotConfig &config);

  /**
   * @brief Applies viewport translation pan by delta pixel displacement.
   * @param[in] delta_x Screen horizontal pixel shift.
   * @param[in] delta_y Screen vertical pixel shift.
   * @param[in] mouse_x Reference cursor X coordinate.
   * @param[in] mouse_y Reference cursor Y coordinate.
   * @param[in] canvas_w Canvas width in pixels.
   * @param[in] canvas_h Canvas height in pixels.
   * @param[in] hist Active histogram data.
   * @param[in] config Base plotting configuration.
   */
  void applyPan(float delta_x, float delta_y, float mouse_x, float mouse_y, float canvas_w,
                float canvas_h, const correlation::analysis::Histogram &hist,
                const correlation::plotters::PlotConfig &config);

  /**
   * @brief Applies stepped or continuous zoom centered around mouse cursor.
   * @param[in] mouse_x Focus cursor X coordinate in pixels.
   * @param[in] mouse_y Focus cursor Y coordinate in pixels.
   * @param[in] delta Zoom displacement factor (>0 zooms in, <0 zooms out).
   * @param[in] canvas_w Canvas width in pixels.
   * @param[in] canvas_h Canvas height in pixels.
   * @param[in] hist Active histogram data.
   * @param[in] config Base plotting configuration.
   */
  void applyWheelZoom(float mouse_x, float mouse_y, float delta, float canvas_w, float canvas_h,
                      const correlation::analysis::Histogram &hist,
                      const correlation::plotters::PlotConfig &config);

  /**
   * @brief Resets zoom state to full automatic bounds.
   */
  void resetZoom() noexcept;

  /**
   * @brief Returns currently active zoom window limits.
   * @return ViewportLimits structure.
   */
  [[nodiscard]] ViewportLimits limits() const noexcept;

  /**
   * @brief Sets explicit manual viewport limits.
   * @param[in] limits Manual boundary limits to apply.
   */
  void setLimits(ViewportLimits limits) noexcept;

  /**
   * @brief Injects active zoom limits into target PlotConfig.
   * @param[in,out] config Configuration structure to enrich.
   */
  void applyToConfig(correlation::plotters::PlotConfig &config) const noexcept;

private:
  std::optional<real_t> zoom_x_min_;
  std::optional<real_t> zoom_x_max_;
  std::optional<real_t> zoom_y_min_;
  std::optional<real_t> zoom_y_max_;
};

} // namespace correlation::app
