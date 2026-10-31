/**
 * @file Structure3DController.cpp
 * @brief Implementation of Structure3DController facade for interactive Slint 3D visualization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/Structure3DController.hpp"
#include "AppWindow.h"
#include "app/viewmodel/PhaseExportManager.hpp"
#include "plotters/ColormapService.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <utility>

namespace correlation::app {

Structure3DController::Structure3DController(::AppWindow &window, const ProgramOptions &options)
    : window_(window), options_(options), alive_(std::make_shared<std::atomic<bool>>(true)),
      render_pipeline_(std::make_unique<Structure3DRenderPipeline>(
          [this](const slint::Image &img) {
            safeInvokeOnWindow([this, img] { window_.set_structure_3d_image(img); });
          },
          alive_)) {
  camera_.setViewport(static_cast<float>(viewport_dims_.width),
                      static_cast<float>(viewport_dims_.height));
  camera_.fitToBox(math::Vector3<float>{-10.0F, -10.0F, -10.0F},
                   math::Vector3<float>{10.0F, 10.0F, 10.0F});

  auto cell_provider = [this](std::size_t frame_index)
      -> std::optional<std::pair<core::Cell, analysis::BondCutoffMatrix>> {
    const std::scoped_lock lock(mutex_);
    if (current_trajectory_ && streamer_.totalFrames() > 0) {
      if (frame_index >= current_trajectory_->getFrameCount()) {
        return std::nullopt;
      }
      return std::make_pair(current_trajectory_->getFrame(frame_index),
                            current_trajectory_->getBondCutoffs());
    }
    if (frame_index != 0 || current_cell_.isEmpty()) {
      return std::nullopt;
    }
    return std::make_pair(current_cell_, analysis::BondCutoffMatrix{});
  };

  auto scene_accessor = [this] -> std::shared_ptr<render3d::Scene3D> {
    const std::scoped_lock lock(mutex_);
    return scene_;
  };

  auto active_frame_provider = [this] -> std::size_t {
    const std::scoped_lock lock(mutex_);
    return streamer_.hasTrajectory() ? streamer_.currentFrameIndex() : 0UL;
  };

  auto total_frames_provider = [this] -> std::size_t {
    const std::scoped_lock lock(mutex_);
    if (current_trajectory_) {
      return current_trajectory_->getFrameCount();
    }
    return current_cell_.isEmpty() ? 0UL : 1UL;
  };

  auto on_phase_updated = [this] {
    requestRender();
    syncPhaseSearchToWindow();
  };

  phase_service_ = std::make_unique<PhaseSearchService>(
      std::move(cell_provider), std::move(scene_accessor), std::move(active_frame_provider),
      std::move(total_frames_provider), std::move(on_phase_updated), alive_);

  syncPhaseSearchToWindow();
}

Structure3DController::~Structure3DController() {
  if (alive_ != nullptr) {
    alive_->store(false);
  }
  phase_service_.reset();
  render_pipeline_.reset();
}

void Structure3DController::requestRender(bool immediate) {
  if (render_pipeline_ == nullptr) {
    return;
  }
  const std::scoped_lock lock(mutex_);
  const RenderPipelineSnapshot snapshot{
      .camera = camera_,
      .scene = scene_,
      .is_active = is_active_,
      .render_bonds = render_bonds_,
      .render_wireframe = render_wireframe_,
      .is_interacting = is_interacting_,
  };
  render_pipeline_->updateSnapshot(snapshot, immediate);
}

void Structure3DController::handleViewportResized(float width, float height) {
  if (width < 10.0F || height < 10.0F) {
    return;
  }
  const auto new_width = static_cast<std::uint32_t>(std::round(width));
  const auto new_height = static_cast<std::uint32_t>(std::round(height));

  {
    const std::scoped_lock lock(mutex_);
    if (viewport_dims_.width == new_width && viewport_dims_.height == new_height) {
      return;
    }
    viewport_dims_ = {.width = new_width, .height = new_height};
    camera_.setViewport(static_cast<float>(new_width), static_cast<float>(new_height));
  }
  if (render_pipeline_ != nullptr) {
    render_pipeline_->resize(viewport_dims_);
  }
  requestRender();
}

void Structure3DController::handleRotate(float delta_x, float delta_y) {
  {
    const std::scoped_lock lock(mutex_);
    camera_.orbit(delta_x, delta_y);
  }
  requestRender();
}

void Structure3DController::handlePan(float delta_x, float delta_y) {
  {
    const std::scoped_lock lock(mutex_);
    camera_.pan(delta_x, delta_y);
  }
  requestRender();
}

void Structure3DController::handleZoom(float delta) {
  {
    const std::scoped_lock lock(mutex_);
    camera_.zoom(delta);
  }
  requestRender();
}

void Structure3DController::handleResetView() {
  {
    const std::scoped_lock lock(mutex_);
    if (scene_ != nullptr) {
      camera_.fitToBox(scene_->boundsMin(), scene_->boundsMax());
    } else {
      camera_.fitToBox(math::Vector3<float>{-10.0F, -10.0F, -10.0F},
                       math::Vector3<float>{10.0F, 10.0F, 10.0F});
    }
  }
  requestRender();
}

void Structure3DController::handleToggleBonds(bool enabled) {
  {
    const std::scoped_lock lock(mutex_);
    render_bonds_ = enabled;
    if (render_bonds_ && scene_ != nullptr) {
      scene_->recomputeBonds(bond_cutoff_scale_, in_cell_mode_);
    }
  }
  safeInvokeOnWindow([this, enabled] { window_.set_structure_3d_bonds_enabled(enabled); });
  requestRender();
}

void Structure3DController::handleToggleWireframe(bool enabled) {
  {
    const std::scoped_lock lock(mutex_);
    render_wireframe_ = enabled;
  }
  safeInvokeOnWindow([this, enabled] { window_.set_structure_3d_wireframe_enabled(enabled); });
  requestRender();
}

void Structure3DController::handleToggleInCell(bool enabled) {
  {
    const std::scoped_lock lock(mutex_);
    in_cell_mode_ = enabled;
    if (scene_ != nullptr) {
      if (streamer_.hasTrajectory() && streamer_.totalFrames() > 0) {
        static_cast<void>(streamer_.seekFrame(streamer_.currentFrameIndex(), *scene_, render_bonds_,
                                              bond_cutoff_scale_, in_cell_mode_));
      } else {
        scene_ = render3d::Scene3D::fromCell(current_cell_, bond_cutoff_scale_, render_bonds_,
                                             in_cell_mode_);
      }
      updateScenePhaseVisuals();
    }
  }
  safeInvokeOnWindow([this, enabled] { window_.set_structure_3d_in_cell(enabled); });
  requestRender();
}

void Structure3DController::handleSeekFrame(int frame_index) {
  if (frame_index < 0) {
    return;
  }
  int actual_frame = frame_index;
  {
    const std::scoped_lock lock(mutex_);
    if (streamer_.totalFrames() > 0 && scene_ != nullptr) {
      const auto target = std::min(static_cast<size_t>(frame_index), streamer_.totalFrames() - 1);
      static_cast<void>(
          streamer_.seekFrame(target, *scene_, render_bonds_, bond_cutoff_scale_, in_cell_mode_));
      actual_frame = static_cast<int>(streamer_.currentFrameIndex());
      updateScenePhaseVisuals();
    }
  }
  safeInvokeOnWindow(
      [this, actual_frame] { window_.set_structure_3d_current_frame(actual_frame); });
  requestRender();
}

void Structure3DController::handleStepFrame(bool forward) {
  int current_idx = 0;
  {
    const std::scoped_lock lock(mutex_);
    if (streamer_.totalFrames() > 0 && scene_ != nullptr) {
      if (forward) {
        static_cast<void>(
            streamer_.stepNext(*scene_, true, render_bonds_, bond_cutoff_scale_, in_cell_mode_));
      } else {
        static_cast<void>(
            streamer_.stepPrev(*scene_, true, render_bonds_, bond_cutoff_scale_, in_cell_mode_));
      }
      current_idx = static_cast<int>(streamer_.currentFrameIndex());
      updateScenePhaseVisuals();
    }
  }
  safeInvokeOnWindow([this, current_idx] { window_.set_structure_3d_current_frame(current_idx); });
  requestRender();
}

void Structure3DController::updateBondCutoffScale(double cutoff_scale) {
  if (cutoff_scale <= 0.0) {
    return;
  }
  {
    const std::scoped_lock lock(mutex_);
    bond_cutoff_scale_ = cutoff_scale;
    if (scene_ != nullptr && render_bonds_) {
      scene_->recomputeBonds(bond_cutoff_scale_, in_cell_mode_);
    }
  }
  requestRender();
}

void Structure3DController::setCell(const core::Cell &cell) {
  current_cell_ = cell;
  auto new_scene =
      render3d::Scene3D::fromCell(cell, bond_cutoff_scale_, render_bonds_, in_cell_mode_);
  {
    const std::scoped_lock lock(mutex_);
    current_trajectory_.reset();
    if (phase_service_ != nullptr) {
      phase_service_->clearCache();
    }
    scene_ = std::move(new_scene);
    if (scene_ != nullptr) {
      camera_.fitToBox(scene_->boundsMin(), scene_->boundsMax());
      updateScenePhaseVisuals();
    }
  }
  safeInvokeOnWindow([this] {
    window_.set_structure_3d_current_frame(0);
    window_.set_structure_3d_total_frames(1);
  });
  requestRender();
  syncPhaseSearchToWindow();
}

void Structure3DController::setTrajectory(std::shared_ptr<const core::Trajectory> trajectory) {
  if (!trajectory) {
    return;
  }
  const auto total_frames = trajectory->getFrameCount();
  if (total_frames == 0) {
    return;
  }

  current_trajectory_ = trajectory;
  {
    const std::scoped_lock lock(mutex_);
    if (phase_service_ != nullptr) {
      phase_service_->clearCache();
    }
    const auto &first_cell = trajectory->getFrame(0);
    current_cell_ = first_cell;
    scene_ =
        render3d::Scene3D::fromCell(first_cell, bond_cutoff_scale_, render_bonds_, in_cell_mode_);
    if (scene_ != nullptr) {
      camera_.fitToBox(scene_->boundsMin(), scene_->boundsMax());
    }
    streamer_.setTrajectory(std::move(trajectory));
    updateScenePhaseVisuals();
  }
  safeInvokeOnWindow([this, total_frames] {
    window_.set_structure_3d_current_frame(0);
    window_.set_structure_3d_total_frames(static_cast<int>(total_frames));
  });
  requestRender();
  syncPhaseSearchToWindow();
}

void Structure3DController::set3DModeActive(bool active) {
  {
    const std::scoped_lock lock(mutex_);
    if (is_active_ == active) {
      return;
    }
    is_active_ = active;
  }
  if (render_pipeline_ != nullptr) {
    render_pipeline_->setActive(active);
  }
  requestRender();
}

void Structure3DController::setPhaseSearchEnabled(bool enabled) {
  if (phase_service_ != nullptr) {
    phase_service_->setEnabled(enabled);
    updateScenePhaseVisuals();
  }
  requestRender();
  syncPhaseSearchToWindow();
}

bool Structure3DController::isPhaseSearchEnabled() const noexcept {
  return phase_service_ ? phase_service_->isEnabled() : false;
}

void Structure3DController::setPhaseDescriptor(PhaseDescriptor descriptor) {
  if (phase_service_ != nullptr) {
    phase_service_->setDescriptor(descriptor);
    updateScenePhaseVisuals();
  }
  requestRender();
  syncPhaseSearchToWindow();
}

PhaseDescriptor Structure3DController::phaseDescriptor() const noexcept {
  return phase_service_ ? phase_service_->descriptor() : PhaseDescriptor::Q6Bar;
}

void Structure3DController::setCrystalPreset(CrystalPreset preset) {
  if (phase_service_ != nullptr) {
    phase_service_->setCrystalPreset(preset);
    updateScenePhaseVisuals();
  }
  requestRender();
  syncPhaseSearchToWindow();
}

CrystalPreset Structure3DController::crystalPreset() const noexcept {
  return phase_service_ ? phase_service_->crystalPreset() : CrystalPreset::Custom;
}

void Structure3DController::setPhaseThresholds(float min_threshold, float max_threshold) {
  if (phase_service_ != nullptr) {
    phase_service_->setThresholds(min_threshold, max_threshold);
    updateScenePhaseVisuals();
  }
  requestRender();
  syncPhaseSearchToWindow();
}

float Structure3DController::phaseMinThreshold() const noexcept {
  return phase_service_ ? phase_service_->minThreshold() : 0.0F;
}

float Structure3DController::phaseMaxThreshold() const noexcept {
  return phase_service_ ? phase_service_->maxThreshold() : 1.0F;
}

void Structure3DController::setPhaseColormap(plotters::ColormapType colormap) {
  if (phase_service_ != nullptr) {
    phase_service_->setColormap(colormap);
    updateScenePhaseVisuals();
  }
  requestRender();
  syncPhaseSearchToWindow();
}

plotters::ColormapType Structure3DController::phaseColormap() const noexcept {
  return phase_service_ ? phase_service_->colormap() : plotters::ColormapType::Viridis;
}

void Structure3DController::setPhaseHaloIntensity(float intensity) {
  if (phase_service_ != nullptr) {
    phase_service_->setHaloIntensity(intensity);
    updateScenePhaseVisuals();
  }
  requestRender();
  syncPhaseSearchToWindow();
}

float Structure3DController::phaseHaloIntensity() const noexcept {
  return phase_service_ ? phase_service_->haloIntensity() : 0.8F;
}

void Structure3DController::setPhaseHaloRadiusRatio(float ratio) {
  if (phase_service_ != nullptr) {
    phase_service_->setHaloRadiusRatio(ratio);
    updateScenePhaseVisuals();
  }
  requestRender();
  syncPhaseSearchToWindow();
}

float Structure3DController::phaseHaloRadiusRatio() const noexcept {
  return phase_service_ ? phase_service_->haloRadiusRatio() : 1.5F;
}

void Structure3DController::setAtomScale(float scale) {
  if (phase_service_ != nullptr) {
    phase_service_->setAtomScale(scale);
    updateScenePhaseVisuals();
  }
  requestRender();
  syncPhaseSearchToWindow();
}

float Structure3DController::atomScale() const noexcept {
  return phase_service_ ? phase_service_->atomScale() : 0.35F;
}

void Structure3DController::setPhaseDimNonMatching(bool dim) {
  if (phase_service_ != nullptr) {
    phase_service_->setDimNonMatching(dim);
    updateScenePhaseVisuals();
  }
  requestRender();
  syncPhaseSearchToWindow();
}

bool Structure3DController::isPhaseDimNonMatching() const noexcept {
  return phase_service_ ? phase_service_->isDimNonMatching() : true;
}

void Structure3DController::setPhaseHideNonMatching(bool hide) {
  if (phase_service_ != nullptr) {
    phase_service_->setHideNonMatching(hide);
    updateScenePhaseVisuals();
  }
  requestRender();
  syncPhaseSearchToWindow();
}

bool Structure3DController::isPhaseHideNonMatching() const noexcept {
  return phase_service_ ? phase_service_->isHideNonMatching() : false;
}

PhaseSearchSettings Structure3DController::phaseSearchSettings() const noexcept {
  return phase_service_ ? phase_service_->settings() : PhaseSearchSettings{};
}

PhaseHistogram Structure3DController::phaseHistogram() const {
  return phase_service_ ? phase_service_->latestHistogram() : PhaseHistogram{};
}

PhaseHeatmap2D Structure3DController::phaseHeatmap() const {
  return phase_service_ ? phase_service_->latestHeatmap() : PhaseHeatmap2D{};
}

void Structure3DController::setPhase2DMode(bool is_2d) {
  if (phase_service_ != nullptr) {
    phase_service_->set2DMode(is_2d);
    updateScenePhaseVisuals();
  }
  requestRender();
  syncPhaseSearchToWindow();
}

bool Structure3DController::isPhase2DMode() const noexcept {
  return phase_service_ ? phase_service_->is2DMode() : false;
}

void Structure3DController::setPhaseDescriptorY(PhaseDescriptor desc_y) {
  if (phase_service_ != nullptr) {
    phase_service_->setDescriptorY(desc_y);
    updateScenePhaseVisuals();
  }
  requestRender();
  syncPhaseSearchToWindow();
}

PhaseDescriptor Structure3DController::phaseDescriptorY() const noexcept {
  return phase_service_ ? phase_service_->descriptorY() : PhaseDescriptor::LocalEntropy;
}

void Structure3DController::setPhaseThresholds2D(float min_x, float max_x, float min_y,
                                                 float max_y) {
  if (phase_service_ != nullptr) {
    phase_service_->setThresholds2D(min_x, max_x, min_y, max_y);
    updateScenePhaseVisuals();
  }
  requestRender();
  syncPhaseSearchToWindow();
}

float Structure3DController::phaseMinThresholdY() const noexcept {
  return phase_service_ ? phase_service_->minThresholdY() : 0.0F;
}

float Structure3DController::phaseMaxThresholdY() const noexcept {
  return phase_service_ ? phase_service_->maxThresholdY() : 0.0F;
}

bool Structure3DController::exportHeatmapPng(const std::string &filepath) {
  if (phase_service_ == nullptr) {
    return false;
  }
  const auto heatmap = phase_service_->latestHeatmap();
  const auto settings = phase_service_->settings();
  return PhaseExportManager::exportHeatmapPng(filepath, heatmap, settings.colormap);
}

void Structure3DController::syncPhaseSearchToWindow() {
  PhaseSearchSettings settings_copy;
  PhaseHistogram hist_copy;
  PhaseHeatmap2D heatmap_copy;
  if (phase_service_ != nullptr) {
    settings_copy = phase_service_->settings();
    hist_copy = phase_service_->latestHistogram();
    heatmap_copy = phase_service_->latestHeatmap();
  }

  const auto lut = plotters::ColormapService::getLut(settings_copy.colormap);
  std::array<slint::Rgba8Pixel, 256> pixels{};
  std::size_t pixel_idx = 0;
  for (const auto &lut_entry : lut) {
    if (pixel_idx < pixels.size()) {
      pixels.at(pixel_idx++) =
          slint::Rgba8Pixel{.r = lut_entry.r, .g = lut_entry.g, .b = lut_entry.b, .a = 255};
    }
  }
  const slint::SharedPixelBuffer<slint::Rgba8Pixel> pixel_buffer(256, 1, pixels.data());
  slint::Image colorbar_img(pixel_buffer);

  auto counts_model = std::make_shared<slint::VectorModel<float>>();
  if (!hist_copy.counts.empty()) {
    std::uint32_t max_count = 0;
    for (const auto c : hist_copy.counts) {
      max_count = std::max(max_count, c);
    }
    const float norm_factor = (max_count > 0) ? (1.0F / static_cast<float>(max_count)) : 0.0F;
    for (const auto c : hist_copy.counts) {
      counts_model->push_back(static_cast<float>(c) * norm_factor);
    }
  }

  slint::Image heatmap_img =
      PhaseExportManager::renderHeatmapImage(heatmap_copy, settings_copy.colormap);

  auto markers_model = std::make_shared<slint::VectorModel<CrystalMarker>>();
  const auto ref_markers =
      getCrystalReferenceMarkers(settings_copy.descriptor, settings_copy.descriptor_y);
  for (const auto &ref : ref_markers) {
    markers_model->push_back(CrystalMarker{
        .label = slint::SharedString(ref.label),
        .x_ratio = ref.x_ratio,
        .y_ratio = ref.y_ratio,
    });
  }

  auto apply_to_window = [this, settings = settings_copy, hist = hist_copy, heatmap = heatmap_copy,
                          colorbar_img = std::move(colorbar_img),
                          counts_model = std::move(counts_model),
                          heatmap_img = std::move(heatmap_img),
                          markers_model = std::move(markers_model)] mutable {
    window_.set_structure_3d_phase_enabled(settings.enabled);
    window_.set_structure_3d_phase_preset(static_cast<int>(settings.preset));
    window_.set_structure_3d_phase_descriptor(static_cast<int>(settings.descriptor));
    window_.set_structure_3d_phase_colormap(static_cast<int>(settings.colormap));
    window_.set_structure_3d_phase_min_threshold(settings.min_threshold);
    window_.set_structure_3d_phase_max_threshold(settings.max_threshold);
    window_.set_structure_3d_phase_domain_min(hist.min_val);
    window_.set_structure_3d_phase_domain_max(hist.max_val);
    window_.set_structure_3d_phase_halo_intensity(settings.halo_intensity);
    window_.set_structure_3d_phase_halo_radius_ratio(settings.halo_radius_ratio);
    window_.set_structure_3d_atom_scale(settings.atom_scale);
    window_.set_structure_3d_phase_dim_non_matching(settings.dim_non_matching);
    window_.set_structure_3d_phase_hide_non_matching(settings.hide_non_matching);
    window_.set_structure_3d_phase_matching_atoms(static_cast<int>(
        settings.is_2d_mode ? heatmap.matching_atom_count : hist.matching_atom_count));
    window_.set_structure_3d_phase_total_atoms(
        static_cast<int>(settings.is_2d_mode ? heatmap.total_atom_count : hist.total_atom_count));
    window_.set_structure_3d_phase_descriptor_name(
        slint::SharedString(descriptorToString(settings.descriptor)));
    window_.set_structure_3d_phase_colorbar_image(colorbar_img);
    window_.set_structure_3d_phase_histogram_counts(counts_model);

    window_.set_structure_3d_phase_is_2d_mode(settings.is_2d_mode);
    window_.set_structure_3d_phase_descriptor_y(static_cast<int>(settings.descriptor_y));
    window_.set_structure_3d_phase_min_threshold_y(settings.min_threshold_y);
    window_.set_structure_3d_phase_max_threshold_y(settings.max_threshold_y);
    window_.set_structure_3d_phase_domain_min_y(heatmap.min_y);
    window_.set_structure_3d_phase_domain_max_y(heatmap.max_y);
    window_.set_structure_3d_phase_descriptor_name_y(
        slint::SharedString(descriptorToString(settings.descriptor_y)));
    window_.set_structure_3d_phase_heatmap_image(heatmap_img);
    window_.set_structure_3d_phase_crystal_markers(markers_model);
  };

  safeInvokeOnWindow(std::move(apply_to_window));
}

std::shared_ptr<const FramePhaseDescriptors>
Structure3DController::getOrComputeDescriptors(std::size_t frame_index) {
  return phase_service_ ? phase_service_->getOrComputeDescriptors(frame_index) : nullptr;
}

std::shared_ptr<const calculators::SteinhardtCalculator::SteinhardtParams>
Structure3DController::getOrComputeSteinhardt(std::size_t frame_index) {
  return phase_service_ ? phase_service_->getOrComputeSteinhardt(frame_index) : nullptr;
}

void Structure3DController::requestPhaseCalculation(std::size_t frame_index) {
  if (phase_service_ != nullptr) {
    phase_service_->requestCalculation(frame_index);
  }
}

void Structure3DController::updateScenePhaseVisuals() {
  if (scene_ == nullptr || phase_service_ == nullptr) {
    return;
  }
  const auto active_frame = streamer_.hasTrajectory() ? streamer_.currentFrameIndex() : 0UL;
  phase_service_->applyVisualsToScene(*scene_, active_frame);
}

} // namespace correlation::app
