/**
 * @file Structure3DController.hpp
 * @brief Viewmodel controller coordinating 3D molecular visualization and Slint MVVM.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "app/core/AppOptions.hpp"
#include "app/viewmodel/PhaseSearchTypes.hpp"
#include "app/viewmodel/SteinhardtLRUCache.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "plotters/ColormapService.hpp"
#include "render3d/Camera3D.hpp"
#include "render3d/IRenderEngine.hpp"
#include "render3d/LodController.hpp"
#include "render3d/Scene3D.hpp"
#include "render3d/TrajectoryStreamer.hpp"
#include <slint.h>

#include <atomic>
#include <condition_variable>
#include <cstdint>
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
   * @param[in] window Slint AppWindow reference.
   * @param[in] options Application program options.
   */
  Structure3DController(::AppWindow &window, const ProgramOptions &options);

  /**
   * @brief Destructor. Cleanly stops render worker thread.
   */
  ~Structure3DController();

  /**
   * @brief Schedules an asynchronous 3D frame redraw.
   * @param[in] immediate If true, signals worker without debounce delay.
   */
  void requestRender(bool immediate = false);

  /**
   * @brief Updates viewport dimensions.
   * @param[in] width Layout width in pixels.
   * @param[in] height Layout height in pixels.
   */
  void handleViewportResized(float width, float height);

  /**
   * @brief Handles mouse orbit rotation deltas.
   * @param[in] delta_x Horizontal drag delta in pixels.
   * @param[in] delta_y Vertical drag delta in pixels.
   */
  void handleRotate(float delta_x, float delta_y);

  /**
   * @brief Handles mouse panning deltas.
   * @param[in] delta_x Horizontal drag delta in pixels.
   * @param[in] delta_y Vertical drag delta in pixels.
   */
  void handlePan(float delta_x, float delta_y);

  /**
   * @brief Handles scroll wheel zoom delta.
   * @param[in] delta Wheel scroll delta factor.
   */
  void handleZoom(float delta);

  /**
   * @brief Resets camera view to auto-fit system bounding box.
   */
  void handleResetView();

  /**
   * @brief Toggles rendering of covalent bonds.
   * @param[in] enabled True to render bonds.
   */
  void handleToggleBonds(bool enabled);

  /**
   * @brief Toggles rendering of periodic simulation unit cell wireframe.
   * @param[in] enabled True to render wireframe box.
   */
  void handleToggleWireframe(bool enabled);

  /**
   * @brief Toggles in-cell PBC coordinate wrapping and boundary-clipped bonds.
   * @param[in] enabled True to wrap atoms and clip crossing bonds to unit cell faces.
   */
  void handleToggleInCell(bool enabled);

  /**
   * @brief Seeks to target frame in active trajectory.
   * @param[in] frame_index Zero-based frame index.
   */
  void handleSeekFrame(int frame_index);

  /**
   * @brief Steps forward or backward through trajectory frames.
   * @param[in] forward True for next frame, false for previous frame.
   */
  void handleStepFrame(bool forward);

  /**
   * @brief Queries whether in-cell PBC wrapping mode is active.
   * @return True if in-cell wrapping and boundary clipping are enabled.
   */
  [[nodiscard]] bool isInCell() const noexcept { return in_cell_mode_; }

  /**
   * @brief Updates active atomic structure from loaded Cell.
   * @param[in] cell New atomic simulation cell.
   */
  void setCell(const core::Cell &cell);

  /**
   * @brief Updates active trajectory dataset.
   * @param[in] trajectory Shared pointer to loaded trajectory.
   */
  void setTrajectory(std::shared_ptr<const core::Trajectory> trajectory);

  /**
   * @brief Sets whether the 3D Structure tab is currently active.
   * @param[in] active True if 3D viewport is visible on screen.
   */
  void set3DModeActive(bool active);

  /**
   * @brief Updates the covalent radius scaling factor used for bond generation.
   * @param[in] cutoff_scale Multiplier for covalent radius bond cutoff.
   */
  void updateBondCutoffScale(double cutoff_scale);

  /**
   * @brief Queries whether the 3D Structure tab is active.
   * @return True if 3D viewport is visible and receiving updates.
   */
  [[nodiscard]] bool is3DModeActive() const noexcept { return is_active_; }

  /**
   * @brief Gets reference to the active orbit camera.
   * @return Const reference to Camera3D instance.
   */
  [[nodiscard]] const render3d::Camera3D &camera() const noexcept { return camera_; }

  /**
   * @brief Gets shared pointer to the active render scene.
   * @return Shared pointer to const Scene3D instance.
   */
  [[nodiscard]] std::shared_ptr<const render3d::Scene3D> scene() const noexcept { return scene_; }

  /**
   * @brief Enables or disables phase search halo visualization.
   * @param[in] enabled True to activate phase halos and non-matching dimming.
   */
  void setPhaseSearchEnabled(bool enabled);

  /**
   * @brief Queries whether phase search halo visualization is active.
   * @return True if phase halos are enabled.
   */
  [[nodiscard]] bool isPhaseSearchEnabled() const noexcept;

  /**
   * @brief Sets the active Steinhardt order parameter descriptor.
   * @param[in] descriptor Active descriptor (e.g. Q6Bar, W6Hat).
   */
  void setPhaseDescriptor(PhaseDescriptor descriptor);

  /**
   * @brief Gets current active Steinhardt order parameter descriptor.
   * @return Active descriptor enum value.
   */
  [[nodiscard]] PhaseDescriptor phaseDescriptor() const noexcept;

  /**
   * @brief Applies a crystal structure preset and auto-configures thresholds.
   * @param[in] preset Preset enum value (Liquid, FCC, HCP, BCC, Icosahedral).
   */
  void setCrystalPreset(CrystalPreset preset);

  /**
   * @brief Gets current crystal preset mode.
   * @return Preset enum value.
   */
  [[nodiscard]] CrystalPreset crystalPreset() const noexcept;

  /**
   * @brief Sets filtering threshold boundaries.
   * @param[in] min_threshold Lower scalar threshold.
   * @param[in] max_threshold Upper scalar threshold.
   */
  void setPhaseThresholds(float min_threshold, float max_threshold);

  /**
   * @brief Gets current lower scalar threshold.
   * @return Min threshold.
   */
  [[nodiscard]] float phaseMinThreshold() const noexcept;

  /**
   * @brief Gets current upper scalar threshold.
   * @return Max threshold.
   */
  [[nodiscard]] float phaseMaxThreshold() const noexcept;

  /**
   * @brief Sets scientific colormap for halo gradient.
   * @param[in] colormap Colormap enum value (Viridis, Plasma, Coolwarm, Turbo).
   */
  void setPhaseColormap(plotters::ColormapType colormap);

  /**
   * @brief Gets active colormap type.
   * @return Colormap type.
   */
  [[nodiscard]] plotters::ColormapType phaseColormap() const noexcept;

  /**
   * @brief Sets halo alpha glow intensity [0.1, 1.0].
   * @param[in] intensity Halo alpha factor.
   */
  void setPhaseHaloIntensity(float intensity);

  /**
   * @brief Gets current halo intensity factor.
   * @return Halo intensity.
   */
  [[nodiscard]] float phaseHaloIntensity() const noexcept;

  /**
   * @brief Sets halo radius multiplier ratio [1.1, 2.0].
   * @param[in] ratio Ratio of halo radius over visual sphere radius.
   */
  void setPhaseHaloRadiusRatio(float ratio);

  /**
   * @brief Gets current halo radius ratio.
   * @return Radius ratio.
   */
  [[nodiscard]] float phaseHaloRadiusRatio() const noexcept;

  /**
   * @brief Toggles ghosting/dimming of atoms outside threshold boundaries.
   * @param[in] dim True to apply dimming to non-matching atoms.
   */
  /**
   * @brief Toggles ghosting/dimming of atoms outside threshold boundaries.
   * @param[in] dim True to apply dimming to non-matching atoms.
   */
  void setPhaseDimNonMatching(bool dim);

  /**
   * @brief Queries whether non-matching atom dimming is enabled.
   * @return True if dimming is enabled.
   */
  [[nodiscard]] bool isPhaseDimNonMatching() const noexcept;

  /**
   * @brief Toggles hiding/culling of atoms outside threshold boundaries.
   * @param[in] hide True to cull non-matching atoms.
   */
  void setPhaseHideNonMatching(bool hide);

  /**
   * @brief Queries whether non-matching atom hiding is enabled.
   * @return True if hiding is enabled.
   */
  [[nodiscard]] bool isPhaseHideNonMatching() const noexcept;

  /**
   * @brief Gets snapshot copy of entire phase search configuration bundle.
   * @return PhaseSearchSettings struct.
   */
  [[nodiscard]] PhaseSearchSettings phaseSearchSettings() const noexcept;

  /**
   * @brief Retrieves latest computed 1D distribution histogram.
   * @return PhaseHistogram telemetry struct.
   */
  [[nodiscard]] PhaseHistogram phaseHistogram() const;

  /**
   * @brief Enables or disables 2D joint distribution filtering mode.
   * @param[in] is_2d True to activate 2D bounding box phase search.
   */
  void setPhase2DMode(bool is_2d);

  /**
   * @brief Queries whether 2D phase search mode is active.
   * @return True if 2D mode is enabled.
   */
  [[nodiscard]] bool isPhase2DMode() const noexcept;

  /**
   * @brief Sets active Y-axis descriptor for 2D mode.
   * @param[in] desc_y Y-axis order descriptor.
   */
  void setPhaseDescriptorY(PhaseDescriptor desc_y);

  /**
   * @brief Gets active Y-axis order descriptor.
   * @return Active Y descriptor enum value.
   */
  [[nodiscard]] PhaseDescriptor phaseDescriptorY() const noexcept;

  /**
   * @brief Sets 2D bounding box thresholds [Xmin, Xmax] x [Ymin, Ymax].
   * @param[in] min_x Lower scalar threshold on X.
   * @param[in] max_x Upper scalar threshold on X.
   * @param[in] min_y Lower scalar threshold on Y.
   * @param[in] max_y Upper scalar threshold on Y.
   */
  void setPhaseThresholds2D(float min_x, float max_x, float min_y, float max_y);

  /**
   * @brief Gets current lower scalar threshold on Y.
   * @return Min Y threshold.
   */
  [[nodiscard]] float phaseMinThresholdY() const noexcept;

  /**
   * @brief Gets current upper scalar threshold on Y.
   * @return Max Y threshold.
   */
  [[nodiscard]] float phaseMaxThresholdY() const noexcept;

  /**
   * @brief Retrieves latest computed 64x64 2D joint distribution heatmap.
   * @return PhaseHeatmap2D telemetry struct.
   */
  [[nodiscard]] PhaseHeatmap2D phaseHeatmap() const;

  /**
   * @brief Exports the active 2D heatmap to a PNG image file.
   * @param[in] filepath Destination file path (defaults to heatmap_snapshot.png if empty).
   * @return True if image file was successfully written.
   */
  bool exportHeatmapPng(const std::string &filepath = "");

  /**
   * @brief Pushes current phase search settings, colormap texture, and histogram telemetry to Slint
   * UI.
   */
  void syncPhaseSearchToWindow();

  /**
   * @brief Accesses reference to the frame phase descriptors LRU cache.
   * @return PhaseDescriptorsLRUCache reference.
   */
  [[nodiscard]] PhaseDescriptorsLRUCache &lruCache() noexcept { return lru_cache_; }

  /**
   * @brief Accesses const reference to the frame phase descriptors LRU cache.
   * @return Const PhaseDescriptorsLRUCache reference.
   */
  [[nodiscard]] const PhaseDescriptorsLRUCache &lruCache() const noexcept { return lru_cache_; }

  /**
   * @brief Synchronously computes or retrieves phase descriptors for a frame.
   * @param[in] frame_index Zero-based frame index.
   * @return Shared pointer to FramePhaseDescriptors, or nullptr if unavailable.
   */
  [[nodiscard]] std::shared_ptr<const FramePhaseDescriptors>
  getOrComputeDescriptors(std::size_t frame_index);

  /**
   * @brief Synchronously computes or retrieves Steinhardt parameters for a frame.
   * @param[in] frame_index Zero-based frame index.
   * @return Shared pointer to SteinhardtParams, or nullptr if unavailable.
   */
  [[nodiscard]] std::shared_ptr<const calculators::SteinhardtCalculator::SteinhardtParams>
  getOrComputeSteinhardt(std::size_t frame_index);

private:
  void renderLoop();
  void phaseWorkerLoop();
  [[nodiscard]] int64_t waitForNextPhaseFrame();
  void processTargetPhaseFrame(int64_t target_frame);
  void prewarmAdjacentFrames(int64_t target_frame);
  void requestPhaseCalculation(std::size_t frame_index);
  void updateScenePhaseVisuals();
  void applyPhaseVisualsToScene(render3d::Scene3D &target_scene,
                                const FramePhaseDescriptors *descriptors) const;
  [[nodiscard]] PhaseHistogram
  computePhaseHistogram(const FramePhaseDescriptors *descriptors) const;
  [[nodiscard]] PhaseHeatmap2D computePhaseHeatmap(const FramePhaseDescriptors *descriptors) const;
  [[nodiscard]] static slint::Image renderHeatmapImage(const PhaseHeatmap2D &heatmap,
                                                       plotters::ColormapType colormap);
  [[nodiscard]] std::shared_ptr<const FramePhaseDescriptors>
  computeDescriptorsSync(std::size_t frame_index);
  void dispatchFrame(const std::uint8_t *pixel_data, std::uint32_t width, std::uint32_t height);

  template <typename F> void safeInvokeOnWindow(F &&action) {
    if (std::this_thread::get_id() == main_thread_id_) {
      action();
    } else {
      slint::invoke_from_event_loop([alive = alive_, callback = std::forward<F>(action)] mutable {
        if (alive && alive->load()) {
          callback();
        }
      });
    }
  }

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
  bool in_cell_mode_{false};
  bool is_interacting_{false};
  double bond_cutoff_scale_{1.15};
  core::Cell current_cell_;
  std::shared_ptr<const core::Trajectory> current_trajectory_;

  PhaseSearchSettings phase_settings_;
  mutable std::mutex phase_mutex_;
  PhaseDescriptorsLRUCache lru_cache_{32};
  PhaseHistogram latest_histogram_;
  PhaseHeatmap2D latest_heatmap_;
  std::condition_variable phase_cv_;
  std::atomic<bool> stop_phase_worker_{false};
  std::atomic<int64_t> requested_phase_frame_{-1};

  std::mutex mutex_;
  std::condition_variable cv_;
  bool render_pending_{false};
  std::atomic<bool> stop_worker_{false};
  std::shared_ptr<std::atomic<bool>> alive_;
  std::thread::id main_thread_id_{std::this_thread::get_id()};
  std::jthread phase_worker_;
  std::jthread render_worker_;
};

} // namespace correlation::app
