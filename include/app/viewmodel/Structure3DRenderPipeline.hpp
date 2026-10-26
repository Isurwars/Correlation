/**
 * @file Structure3DRenderPipeline.hpp
 * @brief Background asynchronous rendering pipeline and LOD controller for 3D structure
 * visualization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "render3d/Camera3D.hpp"
#include "render3d/IRenderEngine.hpp"
#include "render3d/LodController.hpp"
#include "render3d/Scene3D.hpp"
#include <slint.h>

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>

namespace correlation::app {

/**
 * @struct RenderPipelineSnapshot
 * @brief Thread-safe capture of parameters required to render a 3D structure frame.
 */
struct RenderPipelineSnapshot {
  render3d::Camera3D camera;                ///< Camera viewing transform snapshot.
  std::shared_ptr<render3d::Scene3D> scene; ///< Molecular scene geometry to render.
  bool is_active{false};                    ///< Viewport active visibility flag.
  bool render_bonds{true};                  ///< Covalent bond rendering toggle.
  bool render_wireframe{true};              ///< Unit cell wireframe toggle.
  bool is_interacting{false};               ///< True during active user orbit/pan drag.
};

/**
 * @class Structure3DRenderPipeline
 * @brief Manages background 3D offscreen rendering worker thread, LOD evaluation, and Slint image
 * dispatch.
 */
class Structure3DRenderPipeline {
public:
  using FrameDispatchedCallback = std::function<void(const slint::Image &)>;

  /**
   * @brief Constructs pipeline and starts background rendering worker thread.
   * @param[in] callback Callback invoked with rendered Slint image on Slint event loop.
   * @param[in] alive Liveness flag to avoid dispatching after parent controller destruction.
   * @param[in] engine Optional custom IRenderEngine instance (defaults to auto backend).
   */
  explicit Structure3DRenderPipeline(FrameDispatchedCallback callback,
                                     std::shared_ptr<std::atomic<bool>> alive = nullptr,
                                     std::unique_ptr<render3d::IRenderEngine> engine = nullptr);

  /**
   * @brief Destructor. Signals stop condition and joins worker thread safely.
   */
  ~Structure3DRenderPipeline();

  Structure3DRenderPipeline(const Structure3DRenderPipeline &) = delete;
  Structure3DRenderPipeline &operator=(const Structure3DRenderPipeline &) = delete;
  Structure3DRenderPipeline(Structure3DRenderPipeline &&) = delete;
  Structure3DRenderPipeline &operator=(Structure3DRenderPipeline &&) = delete;

  /**
   * @brief Initializes the rendering engine and framebuffer with dimensions.
   * @param[in] dimensions Viewport dimensions in pixels.
   */
  void initialize(const render3d::ViewportDimensions &dimensions);

  /**
   * @brief Resizes the underlying framebuffer and viewport.
   * @param[in] dimensions New viewport dimensions.
   */
  void resize(const render3d::ViewportDimensions &dimensions);

  /**
   * @brief Updates active render parameters and triggers an asynchronous frame redraw.
   * @param[in] snapshot Snapshot bundle containing camera, scene, and feature toggles.
   * @param[in] immediate If true, signals worker without debounce delay.
   */
  void updateSnapshot(const RenderPipelineSnapshot &snapshot, bool immediate = false);

  /**
   * @brief Schedules an asynchronous 3D frame redraw with currently stored parameters.
   * @param[in] immediate If true, signals worker without debounce delay.
   */
  void requestRender(bool immediate = false);

  /**
   * @brief Current viewport pixel dimensions.
   * @return Viewport dimensions.
   */
  [[nodiscard]] render3d::ViewportDimensions dimensions() const noexcept;

  /**
   * @brief Checks if the pipeline viewport is currently active.
   * @return True if active.
   */
  [[nodiscard]] bool isActive() const noexcept;

  /**
   * @brief Sets viewport active status.
   * @param[in] active True if viewport is visible and receiving updates.
   */
  void setActive(bool active) noexcept;

private:
  void renderLoop();
  void dispatchFrame(const std::uint8_t *pixel_data, std::uint32_t width, std::uint32_t height);

  FrameDispatchedCallback callback_;
  std::shared_ptr<std::atomic<bool>> alive_;
  std::unique_ptr<render3d::IRenderEngine> engine_;
  render3d::LodController lod_controller_;
  render3d::ViewportDimensions viewport_dims_{.width = 800, .height = 600};

  mutable std::mutex mutex_;
  std::condition_variable cv_;
  RenderPipelineSnapshot current_snapshot_;
  bool render_pending_{false};
  std::atomic<bool> stop_worker_{false};
  std::jthread render_worker_;
};

} // namespace correlation::app
