/**
 * @file PhaseSearchService.hpp
 * @brief Service managing order parameter calculations, descriptor caching, and phase visuals.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/StructureAnalyzer.hpp"
#include "app/viewmodel/PhaseSearchTypes.hpp"
#include "app/viewmodel/SteinhardtLRUCache.hpp"
#include "calculators/order/SteinhardtCalculator.hpp"
#include "core/Cell.hpp"
#include "plotters/ColormapService.hpp"
#include "render3d/Scene3D.hpp"

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <thread>
#include <utility>

namespace correlation::app {

/**
 * @class PhaseSearchService
 * @brief Manages background order parameter calculations, multi-frame caching, and visual updates.
 */
class PhaseSearchService {
public:
  using CellProvider =
      std::function<std::optional<std::pair<core::Cell, analysis::BondCutoffMatrix>>(std::size_t)>;
  using SceneAccessor = std::function<std::shared_ptr<render3d::Scene3D>()>;
  using ActiveFrameProvider = std::function<std::size_t()>;
  using TotalFramesProvider = std::function<std::size_t()>;
  using UpdateCallback = std::function<void()>;

  /**
   * @brief Constructs PhaseSearchService and spawns the background worker.
   * @param cell_provider Callback returning a cell and its bond cutoffs for a frame.
   * @param scene_accessor Callback returning the current 3D scene.
   * @param active_frame_provider Callback returning the active frame index.
   * @param total_frames_provider Callback returning the total number of frames.
   * @param on_phase_updated Callback invoked when asynchronous calculation updates active scene.
   * @param alive Shared cancellation token indicating controller lifecycle.
   */
  PhaseSearchService(CellProvider cell_provider, SceneAccessor scene_accessor,
                     ActiveFrameProvider active_frame_provider,
                     TotalFramesProvider total_frames_provider, UpdateCallback on_phase_updated,
                     std::shared_ptr<std::atomic<bool>> alive);

  /**
   * @brief Destructor stopping the background worker thread.
   */
  ~PhaseSearchService();

  PhaseSearchService(const PhaseSearchService &) = delete;
  PhaseSearchService &operator=(const PhaseSearchService &) = delete;
  PhaseSearchService(PhaseSearchService &&) = delete;
  PhaseSearchService &operator=(PhaseSearchService &&) = delete;

  // Settings & presets
  void setEnabled(bool enabled);
  [[nodiscard]] bool isEnabled() const noexcept;

  void setDescriptor(PhaseDescriptor descriptor);
  [[nodiscard]] PhaseDescriptor descriptor() const noexcept;

  void setCrystalPreset(CrystalPreset preset);
  [[nodiscard]] CrystalPreset crystalPreset() const noexcept;

  void setThresholds(float min_threshold, float max_threshold);
  [[nodiscard]] float minThreshold() const noexcept;
  [[nodiscard]] float maxThreshold() const noexcept;

  void setColormap(plotters::ColormapType colormap);
  [[nodiscard]] plotters::ColormapType colormap() const noexcept;

  void setHaloIntensity(float intensity);
  [[nodiscard]] float haloIntensity() const noexcept;

  void setHaloRadiusRatio(float ratio);
  [[nodiscard]] float haloRadiusRatio() const noexcept;

  void setAtomScale(float scale);
  [[nodiscard]] float atomScale() const noexcept;

  void setDimNonMatching(bool dim);
  [[nodiscard]] bool isDimNonMatching() const noexcept;

  void setHideNonMatching(bool hide);
  [[nodiscard]] bool isHideNonMatching() const noexcept;

  void set2DMode(bool is_2d);
  [[nodiscard]] bool is2DMode() const noexcept;

  void setDescriptorY(PhaseDescriptor desc_y);
  [[nodiscard]] PhaseDescriptor descriptorY() const noexcept;

  void setThresholds2D(float min_x, float max_x, float min_y, float max_y);
  [[nodiscard]] float minThresholdY() const noexcept;
  [[nodiscard]] float maxThresholdY() const noexcept;

  [[nodiscard]] const PhaseSearchSettings &settings() const noexcept;
  [[nodiscard]] const PhaseHistogram &latestHistogram() const noexcept;
  [[nodiscard]] const PhaseHeatmap2D &latestHeatmap() const noexcept;
  [[nodiscard]] const PhaseDescriptorsLRUCache &lruCache() const noexcept;
  [[nodiscard]] PhaseDescriptorsLRUCache &lruCache() noexcept;

  void requestCalculation(std::size_t frame_index);
  void applyVisualsToScene(render3d::Scene3D &scene, std::size_t active_frame);
  void clearCache();

  [[nodiscard]] std::shared_ptr<const FramePhaseDescriptors>
  getOrComputeDescriptors(std::size_t frame_index);
  [[nodiscard]] std::shared_ptr<const calculators::SteinhardtCalculator::SteinhardtParams>
  getOrComputeSteinhardt(std::size_t frame_index);

  [[nodiscard]] PhaseHistogram
  computePhaseHistogram(const FramePhaseDescriptors *descriptors) const;
  [[nodiscard]] PhaseHeatmap2D computePhaseHeatmap(const FramePhaseDescriptors *descriptors) const;
  [[nodiscard]] std::shared_ptr<const FramePhaseDescriptors>
  computeDescriptorsSync(std::size_t frame_index);

private:
  void phaseWorkerLoop();
  [[nodiscard]] int64_t waitForNextPhaseFrame();
  void processTargetPhaseFrame(int64_t target_frame);
  void prewarmAdjacentFrames(std::size_t center_frame);

  void applyPhaseVisualsToSceneInternal(render3d::Scene3D &scene,
                                        const FramePhaseDescriptors *descriptors);

  CellProvider cell_provider_;
  SceneAccessor scene_accessor_;
  ActiveFrameProvider active_frame_provider_;
  TotalFramesProvider total_frames_provider_;
  UpdateCallback on_phase_updated_;
  std::shared_ptr<std::atomic<bool>> alive_;

  PhaseSearchSettings phase_settings_;
  PhaseHistogram latest_histogram_;
  PhaseHeatmap2D latest_heatmap_;
  PhaseDescriptorsLRUCache lru_cache_;

  mutable std::mutex phase_mutex_;
  std::condition_variable phase_cv_;
  std::atomic<bool> stop_phase_worker_{false};
  std::atomic<int64_t> requested_phase_frame_{-1};
  std::jthread phase_worker_;
};

} // namespace correlation::app
