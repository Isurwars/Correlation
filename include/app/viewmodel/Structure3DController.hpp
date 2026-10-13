/**
 * @file Structure3DController.hpp
 * @brief Viewmodel controller coordinating 3D molecular visualization and Slint MVVM.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "app/core/AppOptions.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "render3d/Camera3D.hpp"
#include "render3d/IRenderEngine.hpp"
#include "render3d/LodController.hpp"
#include "render3d/Scene3D.hpp"
#include "render3d/TrajectoryStreamer.hpp"
#include <slint.h>

#include <atomic>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <thread>

class AppWindow;

namespace correlation::app {

/**
 * @class Structure3DController
 * @brief Coordinates 3D orbit interactions, scene streaming, and asynchronous off-thread rendering.
 */
class Structure3DController {
public:
  /**
   * @brief Constructs controller injected with UI window and options.
   * @param window Slint AppWindow reference.
   * @param options Application program options.
   */
  Structure3DController(::AppWindow &window, const ProgramOptions &options);

  /**
   * @brief Destructor. Cleanly stops render worker thread.
   */
  ~Structure3DController();

  /**
   * @brief Schedules an asynchronous 3D frame redraw.
   * @param immediate If true, signals worker without debounce delay.
   */
  void requestRender(bool immediate = false);

  /**
   * @brief Updates viewport dimensions.
   * @param width Layout width in pixels.
   * @param height Layout height in pixels.
   */
  void handleViewportResized(float width, float height);

  /**
   * @brief Handles mouse orbit rotation deltas.
   * @param delta_x Horizontal drag delta in pixels.
   * @param delta_y Vertical drag delta in pixels.
   */
  void handleRotate(float delta_x, float delta_y);

  /**
   * @brief Handles mouse panning deltas.
   * @param delta_x Horizontal drag delta in pixels.
   * @param delta_y Vertical drag delta in pixels.
   */
  void handlePan(float delta_x, float delta_y);

  /**
   * @brief Handles scroll wheel zoom delta.
   * @param delta Wheel scroll delta factor.
   */
  void handleZoom(float delta);

  /**
   * @brief Resets camera view to auto-fit system bounding box.
   */
  void handleResetView();

  /**
   * @brief Toggles rendering of covalent bonds.
   * @param enabled True to render bonds.
   */
  void handleToggleBonds(bool enabled);

  /**
   * @brief Toggles rendering of periodic simulation unit cell wireframe.
   * @param enabled True to render wireframe box.
   */
  void handleToggleWireframe(bool enabled);

  /**
   * @brief Seeks to target frame in active trajectory.
   * @param frame_index Zero-based frame index.
   */
  void handleSeekFrame(int frame_index);

  /**
   * @brief Steps forward or backward through trajectory frames.
   * @param forward True for next frame, false for previous frame.
   */
  void handleStepFrame(bool forward);

  /**
   * @brief Updates active atomic structure from loaded Cell.
   * @param cell New atomic simulation cell.
   */
  void setCell(const core::Cell &cell);

  /**
   * @brief Updates active trajectory dataset.
   * @param trajectory Shared pointer to loaded trajectory.
   */
  void setTrajectory(std::shared_ptr<const core::Trajectory> trajectory);

  /**
   * @brief Sets whether the 3D Structure tab is currently active.
   * @param active True if 3D viewport is visible on screen.
   */
  void set3DModeActive(bool active);

  [[nodiscard]] bool is3DModeActive() const noexcept { return is_active_; }
  [[nodiscard]] const render3d::Camera3D &camera() const noexcept { return camera_; }
  [[nodiscard]] std::shared_ptr<const render3d::Scene3D> scene() const noexcept { return scene_; }

private:
  void renderLoop();
  void dispatchFrame(const uint8_t *pixel_data, uint32_t width, uint32_t height);

  ::AppWindow &window_;
  const ProgramOptions &options_;

  std::unique_ptr<render3d::IRenderEngine> engine_;
  std::shared_ptr<render3d::Scene3D> scene_;
  render3d::Camera3D camera_;
  render3d::LodController lod_controller_;
  render3d::TrajectoryStreamer streamer_;
  render3d::ViewportDimensions viewport_dims_{.width = 800, .height = 600};

  bool is_active_{false};
  bool render_bonds_{true};
  bool render_wireframe_{true};
  bool is_interacting_{false};

  std::mutex mutex_;
  std::condition_variable cv_;
  bool render_pending_{false};
  std::atomic<bool> stop_worker_{false};
  std::shared_ptr<std::atomic<bool>> alive_;
  std::jthread render_worker_;
};

} // namespace correlation::app
