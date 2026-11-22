/**
 * @file Structure3DController.hpp
 * @brief Viewmodel controller coordinating 3D molecular visualization and Slint MVVM.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/PerFrameDistributionStore.hpp"
#include "app/core/AppOptions.hpp"
#include "app/viewmodel/PhaseSearchService.hpp"
#include "app/viewmodel/PhaseSearchTypes.hpp"
#include "app/viewmodel/SteinhardtLRUCache.hpp"
#include "app/viewmodel/Structure3DRenderPipeline.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "plotters/ColormapService.hpp"
#include "render3d/Camera3D.hpp"
#include "render3d/MiniPlotRenderer.hpp"
#include "render3d/Scene3D.hpp"
#include "render3d/TrajectoryStreamer.hpp"
#include <slint.h>

#include <atomic>
#include <chrono>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <unordered_set>

class AppWindow;

namespace correlation::app {

/**
 * @struct SlabClipConfig
 * @brief Normalized [0, 1] X/Y/Z bounding box clipping configuration for the 3D scene.
 */
struct SlabClipConfig {
  bool enabled{false};
  float x_min{0.0F};
  float x_max{1.0F};
  float y_min{0.0F};
  float y_max{1.0F};
  float z_min{0.0F};
  float z_max{1.0F};
};

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
   * @brief Sets camera orientation to a standard orthographic or isometric preset.
   * @param[in] preset_idx 0: Top (XY), 1: Front (XZ), 2: Side (YZ), 3: Isometric.
   */
  void handleSetCameraPreset(int preset_idx);

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
   * @brief Queries whether isometric (orthographic) projection mode is active.
   * @return True if isometric projection is enabled.
   */
  [[nodiscard]] bool isIsometric() const noexcept {
    const std::scoped_lock lock(mutex_);
    return camera_.isOrthographic();
  }

  /**
   * @brief Toggles between perspective and isometric (orthographic) 3D camera projection.
   * @param[in] isometric True for isometric/orthographic, false for perspective.
   */
  void handleToggleProjection(bool isometric);

  /**
   * @brief Toggles atmospheric depth cueing (fog) in Vulkan and CPU renderers.
   * @param[in] enabled True to enable depth cueing.
   */
  void handleToggleDepthCueing(bool enabled);

  /**
   * @brief Queries whether depth cueing is enabled.
   * @return True if depth cueing is active.
   */
  [[nodiscard]] bool isDepthCueingEnabled() const noexcept {
    const std::scoped_lock lock(mutex_);
    return depth_cueing_;
  }

  /**
   * @brief Toggles screen-space ambient occlusion (SSAO) contact darkening on Vulkan.
   * @param[in] enabled True to enable SSAO.
   */
  void handleToggleSsao(bool enabled);

  /**
   * @brief Queries whether SSAO is enabled.
   * @return True if SSAO is enabled.
   */
  [[nodiscard]] bool isSsaoEnabled() const noexcept {
    const std::scoped_lock lock(mutex_);
    return ssao_enabled_;
  }

  /**
   * @brief Sets visual-only periodic supercell replication counts (1..3 per axis).
   */
  void handleSetVisualSupercell(int rep_a, int rep_b, int rep_c);

  /**
   * @brief Resets visual supercell replication back to 1x1x1.
   */
  void handleResetVisualSupercell();

  /**
   * @brief Returns current visual supercell replication counts [nx, ny, nz].
   */
  [[nodiscard]] std::array<int, 3> visualSupercellDims() const noexcept {
    const std::scoped_lock lock(mutex_);
    return supercell_dims_;
  }

  /**
   * @brief Exports all frames of the loaded trajectory as a PNG-sequence ZIP or Animated GIF.
   * @param[in] filepath Destination .zip or .gif path.
   * @param[in] as_gif True for GIF89a animation, false for PNG-sequence ZIP archive.
   * @param[in] fps Playback frame rate for GIF timing (default 15).
   * @param[in] bg_mode Background mode (0 = Dark, 1 = White, 2 = Transparent).
   * @param[out] out_frames Optional count of exported frames.
   * @return True if export succeeded.
   */
  [[nodiscard]] bool exportTrajectoryAnimation(const std::string &filepath, bool as_gif,
                                               int fps = 15, int bg_mode = 0,
                                               std::size_t *out_frames = nullptr);

  /**
   * @brief Smoothly centers and zooms the camera onto the currently picked atom.
   */
  void handleFocusPickedAtom();

  /**
   * @brief Advances active smooth camera interpolation by one frame tick.
   */
  void tickCameraAnimation();

  /**
   * @brief Toggles hide/show visibility for all atoms of a chemical element symbol.
   * @param[in] symbol Chemical element symbol (e.g. "Cu", "Si").
   */
  void handleToggleElementVisibility(const std::string &symbol);

  /**
   * @brief Restores visibility of all chemical elements.
   */
  void handleShowAllElements();

  /**
   * @brief Exports a high-resolution offscreen 3D viewport snapshot to a PNG file.
   * @param[in] filepath Destination PNG file path.
   * @param[in] scale_idx Resolution multiplier index (0 = 1x, 1 = 2x, 2 = 4x).
   * @param[in] bg_mode Background mode (0 = Viewport dark, 1 = White, 2 = Transparent).
   * @param[out] out_width Optional exported pixel width.
   * @param[out] out_height Optional exported pixel height.
   * @return True if PNG file was written successfully.
   */
  [[nodiscard]] bool exportSnapshotPng(const std::string &filepath, int scale_idx = 1,
                                       int bg_mode = 0, std::uint32_t *out_width = nullptr,
                                       std::uint32_t *out_height = nullptr);

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
   * @brief Sets or updates the trajectory-wide per-frame distribution store.
   * @param[in] store Shared pointer to PerFrameDistributionStore.
   */
  void setDistributionStore(std::shared_ptr<const analysis::PerFrameDistributionStore> store);

  /**
   * @brief Returns active per-frame distribution store if available.
   * @return Shared pointer to const PerFrameDistributionStore, or nullptr.
   */
  [[nodiscard]] std::shared_ptr<const analysis::PerFrameDistributionStore>
  distributionStore() const noexcept {
    const std::scoped_lock lock(mutex_);
    return distribution_store_;
  }

  /**
   * @brief Selects active distribution metric displayed in 3D viewer mini-chart HUD.
   * @param[in] metric_index 0 = RDF g(r), 1 = PAD θ, 2 = S(Q).
   */
  void handleSelectDistributionMetric(int metric_index);

  /**
   * @brief Toggles visibility of the distribution mini-chart HUD.
   * @param[in] visible True to display mini-chart HUD overlay.
   */
  void handleToggleDistributionHUD(bool visible);

  /**
   * @brief Queries whether the distribution mini-graph HUD is currently visible.
   * @return True if HUD is visible, false otherwise.
   */
  [[nodiscard]] bool isDistributionHUDVisible() const noexcept {
    const std::scoped_lock lock(mutex_);
    return dist_hud_visible_;
  }

  /**
   * @brief Sets custom plot configuration / limits for a distribution metric on the HUD.
   * @param[in] metric Target distribution metric.
   * @param[in] config Custom mini-plot configuration.
   */
  void setMetricPlotLimits(analysis::DistributionMetric metric,
                           const render3d::MiniPlotConfig &config);

  /**
   * @brief Renders and synchronizes the active frame mini-chart with the Slint UI.
   */
  void syncDistributionHUDToWindow();

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
  [[nodiscard]] std::shared_ptr<const render3d::Scene3D> scene() const noexcept {
    const std::scoped_lock lock(mutex_);
    return scene_;
  }

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
   * @brief Sets halo radius multiplier ratio [1.0, 5.0].
   * @param[in] ratio Ratio of halo radius over visual sphere radius.
   */
  void setPhaseHaloRadiusRatio(float ratio);

  /**
   * @brief Gets current halo radius ratio.
   * @return Radius ratio.
   */
  [[nodiscard]] float phaseHaloRadiusRatio() const noexcept;

  /**
   * @brief Sets uniform atom ball radius scaling multiplier [0.05, 2.0].
   * @param[in] scale Atom radius scale multiplier.
   */
  void setAtomScale(float scale);

  /**
   * @brief Gets current uniform atom ball radius scaling multiplier.
   * @return Atom radius scale multiplier.
   */
  [[nodiscard]] float atomScale() const noexcept;

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
  [[nodiscard]] PhaseDescriptorsLRUCache &lruCache() noexcept { return phase_service_->lruCache(); }

  /**
   * @brief Accesses const reference to the frame phase descriptors LRU cache.
   * @return Const PhaseDescriptorsLRUCache reference.
   */
  [[nodiscard]] const PhaseDescriptorsLRUCache &lruCache() const noexcept {
    return phase_service_->lruCache();
  }

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

  /**
   * @brief Synchronizes projected axis triad and Angstrom scale bar HUD to the Slint window.
   */
  void syncCameraHudToWindow();

  /**
   * @brief Performs ray-sphere atom picking at normalized viewport coordinates [0, 1].
   * @param[in] norm_x Normalized horizontal click position (0 = left, 1 = right).
   * @param[in] norm_y Normalized vertical click position (0 = top, 1 = bottom).
   */
  void handlePickAtom(float norm_x, float norm_y);

  /**
   * @brief Clears the currently inspected atom selection and removes inspection highlights.
   */
  void clearPickedAtom();

  /**
   * @brief Picks an atom for the 3D distance/angle ruler (up to 3 atoms, Shift+Click).
   * @param[in] norm_x Normalized horizontal click position (0 = left, 1 = right).
   * @param[in] norm_y Normalized vertical click position (0 = top, 1 = bottom).
   */
  void handleMeasureAtom(float norm_x, float norm_y);

  /**
   * @brief Clears active 3D distance/angle measurement atoms and guide lines.
   */
  void clearMeasurement();

  /**
   * @brief Returns currently measured atom indices (0..3).
   * @return Copy of measured atom index list.
   */
  [[nodiscard]] std::vector<std::size_t> measuredAtomIndices() const {
    const std::scoped_lock lock(mutex_);
    return measured_atom_indices_;
  }

  /**
   * @brief Returns the currently picked atom index, if any.
   * @return Optional zero-based atom index in Scene3D::atoms().
   */
  [[nodiscard]] std::optional<std::size_t> pickedAtomIndex() const noexcept {
    const std::scoped_lock lock(mutex_);
    return picked_atom_index_;
  }

  /**
   * @brief Updates X/Y/Z normalized slab clipping bounds and re-renders.
   * @param[in] config Slab clipping configuration.
   */
  void setSlabClipConfig(const SlabClipConfig &config);

  /**
   * @brief Resets slab clipping bounds to [0, 1] and disables clipping.
   */
  void resetSlabClip();

  /**
   * @brief Returns the current slab clipping configuration.
   * @return SlabClipConfig copy.
   */
  [[nodiscard]] SlabClipConfig slabClipConfig() const noexcept {
    const std::scoped_lock lock(mutex_);
    return slab_clip_;
  }

  /**
   * @brief Sets or clears the Dual-View radial-shell brushing distance from g(r) hover.
   * @param[in] r_angstrom Hovered radial distance in Angstroms, or std::nullopt to clear.
   */
  void setRadialBrushDistance(std::optional<double> r_angstrom);

private:
  struct CameraAnimState {
    bool active{false};
    std::chrono::steady_clock::time_point start_time{};
    float duration_sec{0.22F};
    math::Vector3<float> start_target{0.0F, 0.0F, 0.0F};
    math::Vector3<float> end_target{0.0F, 0.0F, 0.0F};
    float start_distance{20.0F};
    float end_distance{20.0F};
    float start_azimuth{0.0F};
    float end_azimuth{0.0F};
    float start_elevation{20.0F};
    float end_elevation{20.0F};
  };

  void requestPhaseCalculation(std::size_t frame_index);
  [[nodiscard]] bool canAdvancePlaybackFrameLocked(bool forward) const;
  void startCameraAnimationLocked(const math::Vector3<float> &target, float distance,
                                  float azimuth_deg, float elevation_deg);
  void cancelCameraAnimationLocked();
  void updateScenePhaseVisuals();
  void updateScenePhaseVisualsLocked();
  void applyAtomInspectionHighlightLocked();
  void applyMeasurementHighlightLocked();
  void applyRadialBrushHighlightLocked(std::size_t &matched_atoms);
  void applyElementVisibilityLocked();
  void applySlabClippingLocked();
  void syncInspectedAtomToWindow();
  void syncMeasurementToWindow();
  void syncElementLegendToWindow();

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
  [[maybe_unused]] const ProgramOptions &options_;

  std::shared_ptr<std::atomic<bool>> alive_;
  std::unique_ptr<Structure3DRenderPipeline> render_pipeline_;
  std::unique_ptr<PhaseSearchService> phase_service_;
  std::shared_ptr<render3d::Scene3D> scene_;
  std::shared_ptr<render3d::Scene3D> published_scene_;
  bool scene_dirty_{true};
  render3d::Camera3D camera_;
  render3d::TrajectoryStreamer streamer_;
  render3d::ViewportDimensions viewport_dims_{.width = 800, .height = 600};

  bool is_active_{false};
  bool render_bonds_{true};
  bool render_wireframe_{true};
  bool in_cell_mode_{true};
  bool depth_cueing_{true};
  bool ssao_enabled_{true};
  bool is_interacting_{false};
  std::array<int, 3> supercell_dims_{1, 1, 1};
  CameraAnimState camera_anim_{};
  std::unordered_set<std::string> hidden_elements_;
  double bond_cutoff_scale_{1.15};
  core::Cell current_cell_;
  std::shared_ptr<const core::Trajectory> current_trajectory_;

  std::shared_ptr<const analysis::PerFrameDistributionStore> distribution_store_;
  analysis::DistributionMetric active_dist_metric_{analysis::DistributionMetric::RDF};
  std::array<render3d::MiniPlotConfig, 4> metric_plot_configs_{};
  bool dist_hud_visible_{true};

  std::optional<std::size_t> picked_atom_index_;
  std::vector<std::size_t> measured_atom_indices_;
  SlabClipConfig slab_clip_{};
  std::optional<double> radial_brush_r_;

  mutable std::mutex mutex_;
  std::thread::id main_thread_id_{std::this_thread::get_id()};
};

} // namespace correlation::app
