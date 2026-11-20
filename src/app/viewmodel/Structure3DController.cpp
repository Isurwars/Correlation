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
#include "render3d/RayPicker.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <format>
#include <iostream>
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
    updateScenePhaseVisuals();
    requestRender();
    syncPhaseSearchToWindow();
  };

  phase_service_ = std::make_unique<PhaseSearchService>(
      std::move(cell_provider), std::move(active_frame_provider), std::move(total_frames_provider),
      std::move(on_phase_updated), alive_);

  syncPhaseSearchToWindow();
  syncDistributionHUDToWindow();
  safeInvokeOnWindow([this] {
    window_.set_structure_3d_in_cell(in_cell_mode_);
    window_.set_structure_3d_isometric(camera_.isOrthographic());
  });

  const bool is_gpu = render_pipeline_->isHardwareAccelerated();
  const std::string device_name = render_pipeline_->deviceName();
  const std::string backend_name = render_pipeline_->backendName();
  std::cerr << "[Render3D] Backend: " << backend_name << " | Device: " << device_name << " ("
            << (is_gpu ? "GPU" : "CPU") << ")\n";
  window_.set_structure_3d_gpu_accelerated(is_gpu);
  window_.set_structure_3d_backend_label(
      slint::SharedString(std::format("{}: {}", is_gpu ? "GPU" : "CPU", device_name)));
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
  if (scene_ == nullptr) {
    published_scene_.reset();
    scene_dirty_ = false;
  } else if (scene_dirty_ || published_scene_ == nullptr) {
    published_scene_ = std::make_shared<render3d::Scene3D>(*scene_);
    scene_dirty_ = false;
  }
  const bool gpu_halos = false;
  const RenderPipelineSnapshot snapshot{
      .camera = camera_,
      .scene = published_scene_,
      .is_active = is_active_,
      .render_bonds = render_bonds_,
      .render_wireframe = render_wireframe_,
      .is_interacting = is_interacting_,
      .gpu_halos_enabled = gpu_halos,
      .phase_descriptor = phaseDescriptor(),
      .phase_min_threshold = phaseMinThreshold(),
      .phase_max_threshold = phaseMaxThreshold(),
      .phase_halo_ratio = phaseHaloRadiusRatio(),
      .phase_halo_intensity = phaseHaloIntensity(),
      .phase_halo_color = {1.0F, 0.84F, 0.0F, 1.0F},
      .phase_dim_non_matching = isPhaseDimNonMatching(),
      .phase_hide_non_matching = isPhaseHideNonMatching(),
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
  syncCameraHudToWindow();
  requestRender();
}

void Structure3DController::handleRotate(float delta_x, float delta_y) {
  constexpr float K_ORBIT_SENSITIVITY = 0.4F; // degrees per logical pixel
  {
    const std::scoped_lock lock(mutex_);
    camera_.orbit(delta_x * K_ORBIT_SENSITIVITY, delta_y * K_ORBIT_SENSITIVITY);
  }
  syncCameraHudToWindow();
  requestRender();
}

void Structure3DController::handlePan(float delta_x, float delta_y) {
  {
    const std::scoped_lock lock(mutex_);
    camera_.pan(delta_x, delta_y);
  }
  syncCameraHudToWindow();
  requestRender();
}

void Structure3DController::handleZoom(float delta) {
  constexpr float K_ZOOM_SENSITIVITY = 0.18F; // ~16.5% distance change per wheel notch
  const float clamped_delta = std::clamp(delta, -4.0F, 4.0F);
  {
    const std::scoped_lock lock(mutex_);
    camera_.zoom(clamped_delta * K_ZOOM_SENSITIVITY);
  }
  syncCameraHudToWindow();
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
  syncCameraHudToWindow();
  requestRender();
}

void Structure3DController::handleSetCameraPreset(int preset_idx) {
  {
    const std::scoped_lock lock(mutex_);
    switch (preset_idx) {
    case 0: // Top (XY)
      camera_.setOrientation(0.0F, 89.0F);
      break;
    case 1: // Front (XZ)
      camera_.setOrientation(0.0F, 0.0F);
      break;
    case 2: // Side (YZ)
      camera_.setOrientation(90.0F, 0.0F);
      break;
    case 3: // Isometric
    default:
      camera_.setOrientation(45.0F, 35.26F);
      break;
    }
  }
  syncCameraHudToWindow();
  requestRender();
}

void Structure3DController::handleToggleBonds(bool enabled) {
  {
    const std::scoped_lock lock(mutex_);
    render_bonds_ = enabled;
    if (render_bonds_ && scene_ != nullptr) {
      scene_->recomputeBonds(bond_cutoff_scale_, in_cell_mode_);
      scene_dirty_ = true;
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
      updateScenePhaseVisualsLocked();
    }
  }
  safeInvokeOnWindow([this, enabled] { window_.set_structure_3d_in_cell(enabled); });
  requestRender();
}

void Structure3DController::handleToggleProjection(bool isometric) {
  {
    const std::scoped_lock lock(mutex_);
    camera_.setProjectionMode(isometric ? render3d::ProjectionMode::Orthographic
                                        : render3d::ProjectionMode::Perspective);
  }
  safeInvokeOnWindow([this, isometric] { window_.set_structure_3d_isometric(isometric); });
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
      updateScenePhaseVisualsLocked();
    }
  }
  safeInvokeOnWindow(
      [this, actual_frame] { window_.set_structure_3d_current_frame(actual_frame); });
  syncDistributionHUDToWindow();
  if (isPhaseSearchEnabled()) {
    syncPhaseSearchToWindow();
  }
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
      updateScenePhaseVisualsLocked();
    }
  }
  safeInvokeOnWindow([this, current_idx] { window_.set_structure_3d_current_frame(current_idx); });
  syncDistributionHUDToWindow();
  if (isPhaseSearchEnabled()) {
    syncPhaseSearchToWindow();
  }
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
      scene_dirty_ = true;
    }
  }
  requestRender();
}

void Structure3DController::setCell(const core::Cell &cell) {
  {
    const std::scoped_lock lock(mutex_);
    current_cell_ = cell;
    current_trajectory_.reset();
    picked_atom_index_.reset();
    radial_brush_r_.reset();
    if (phase_service_ != nullptr) {
      phase_service_->clearCache();
    }
    scene_ = render3d::Scene3D::fromCell(cell, bond_cutoff_scale_, render_bonds_, in_cell_mode_);
    if (scene_ != nullptr) {
      camera_.fitToBox(scene_->boundsMin(), scene_->boundsMax());
      updateScenePhaseVisualsLocked();
    }
  }
  safeInvokeOnWindow([this] {
    window_.set_structure_3d_current_frame(0);
    window_.set_structure_3d_total_frames(1);
    window_.set_structure_3d_radial_brush_text("");
  });
  syncInspectedAtomToWindow();
  syncCameraHudToWindow();
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

  {
    const std::scoped_lock lock(mutex_);
    current_trajectory_ = trajectory;
    picked_atom_index_.reset();
    radial_brush_r_.reset();
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
    updateScenePhaseVisualsLocked();
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
        .width_ratio = ref.width_ratio,
        .height_ratio = ref.height_ratio,
        .is_area = ref.is_area,
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
  const std::scoped_lock lock(mutex_);
  updateScenePhaseVisualsLocked();
}

void Structure3DController::updateScenePhaseVisualsLocked() {
  if (scene_ == nullptr) {
    return;
  }
  if (phase_service_ != nullptr) {
    const auto active_frame = streamer_.hasTrajectory() ? streamer_.currentFrameIndex() : 0UL;
    phase_service_->applyVisualsToScene(*scene_, active_frame);
  }
  applyAtomInspectionHighlightLocked();
  std::size_t matched_brush = 0;
  applyRadialBrushHighlightLocked(matched_brush);
  applySlabClippingLocked();
  scene_dirty_ = true;
}

void Structure3DController::syncCameraHudToWindow() {
  render3d::CameraHudMetrics metrics{};
  {
    const std::scoped_lock lock(mutex_);
    metrics = render3d::computeCameraHudMetrics(camera_, static_cast<float>(viewport_dims_.height));
  }
  safeInvokeOnWindow([this, metrics] {
    CameraHudState hud{};
    hud.ax_x = metrics.axis_x[0];
    hud.ax_y = metrics.axis_x[1];
    hud.ay_x = metrics.axis_y[0];
    hud.ay_y = metrics.axis_y[1];
    hud.az_x = metrics.axis_z[0];
    hud.az_y = metrics.axis_z[1];
    hud.scale_bar_px = metrics.scale_bar_px;
    hud.scale_bar_label = slint::SharedString(metrics.scale_bar_label);
    window_.set_structure_3d_camera_hud(hud);
  });
}

void Structure3DController::handlePickAtom(float norm_x, float norm_y) {
  {
    const std::scoped_lock lock(mutex_);
    if (scene_ == nullptr) {
      picked_atom_index_.reset();
    } else {
      picked_atom_index_ = render3d::pickAtom(camera_, *scene_, norm_x, norm_y);
    }
    updateScenePhaseVisualsLocked();
  }
  syncInspectedAtomToWindow();
  requestRender(true);
}

void Structure3DController::clearPickedAtom() {
  {
    const std::scoped_lock lock(mutex_);
    if (!picked_atom_index_.has_value()) {
      return;
    }
    picked_atom_index_.reset();
    updateScenePhaseVisualsLocked();
  }
  syncInspectedAtomToWindow();
  requestRender(true);
}

void Structure3DController::setSlabClipConfig(const SlabClipConfig &config) {
  SlabClipConfig sanitized = config;
  sanitized.x_min = std::clamp(sanitized.x_min, 0.0F, 1.0F);
  sanitized.x_max = std::clamp(sanitized.x_max, sanitized.x_min, 1.0F);
  sanitized.y_min = std::clamp(sanitized.y_min, 0.0F, 1.0F);
  sanitized.y_max = std::clamp(sanitized.y_max, sanitized.y_min, 1.0F);
  sanitized.z_min = std::clamp(sanitized.z_min, 0.0F, 1.0F);
  sanitized.z_max = std::clamp(sanitized.z_max, sanitized.z_min, 1.0F);
  {
    const std::scoped_lock lock(mutex_);
    slab_clip_ = sanitized;
    updateScenePhaseVisualsLocked();
  }
  safeInvokeOnWindow([this, sanitized] {
    SlabClipState ui_state{};
    ui_state.enabled = sanitized.enabled;
    ui_state.x_min = sanitized.x_min;
    ui_state.x_max = sanitized.x_max;
    ui_state.y_min = sanitized.y_min;
    ui_state.y_max = sanitized.y_max;
    ui_state.z_min = sanitized.z_min;
    ui_state.z_max = sanitized.z_max;
    window_.set_structure_3d_slab_clip(ui_state);
  });
  requestRender(true);
}

void Structure3DController::resetSlabClip() { setSlabClipConfig(SlabClipConfig{}); }

void Structure3DController::setRadialBrushDistance(std::optional<double> r_angstrom) {
  std::string brush_label;
  bool changed = false;
  {
    const std::scoped_lock lock(mutex_);
    if (radial_brush_r_.has_value() != r_angstrom.has_value() ||
        (r_angstrom.has_value() && std::abs(*radial_brush_r_ - *r_angstrom) > 0.02)) {
      radial_brush_r_ = r_angstrom;
      changed = true;
      if (scene_ != nullptr) {
        if (phase_service_ != nullptr) {
          const auto active_frame = streamer_.hasTrajectory() ? streamer_.currentFrameIndex() : 0UL;
          phase_service_->applyVisualsToScene(*scene_, active_frame);
        }
        applyAtomInspectionHighlightLocked();
        std::size_t matched_atoms = 0;
        applyRadialBrushHighlightLocked(matched_atoms);
        applySlabClippingLocked();
        if (radial_brush_r_.has_value()) {
          brush_label =
              std::format("Shell r = {:.2f} ± 0.15 Å ({} atoms)", *radial_brush_r_, matched_atoms);
        }
      }
    }
  }
  if (!changed) {
    return;
  }
  safeInvokeOnWindow([this, brush_label] {
    window_.set_structure_3d_radial_brush_text(slint::SharedString(brush_label));
  });
  requestRender(false);
}

void Structure3DController::applyAtomInspectionHighlightLocked() {
  if (scene_ == nullptr || !picked_atom_index_.has_value()) {
    return;
  }
  auto &atoms = scene_->atoms();
  const std::size_t sel_idx = *picked_atom_index_;
  if (sel_idx >= atoms.size()) {
    return;
  }
  const auto &center = atoms[sel_idx];
  const auto scale = static_cast<float>(bond_cutoff_scale_);

  for (std::size_t idx = 0; idx < atoms.size(); ++idx) {
    if (idx == sel_idx) {
      atoms[idx].halo_color = {.r = 64, .g = 215, .b = 215, .a = 230};
      atoms[idx].halo_radius_ratio = 2.4F;
      atoms[idx].dim_factor = 1.0F;
      continue;
    }
    const float cutoff = (center.covalent_radius + atoms[idx].covalent_radius) * scale;
    const float dist = math::norm(atoms[idx].position - center.position);
    if (dist > 0.01F && dist <= cutoff) {
      atoms[idx].halo_color = {.r = 255, .g = 196, .b = 0, .a = 200};
      atoms[idx].halo_radius_ratio = 2.0F;
      atoms[idx].dim_factor = 1.0F;
    } else {
      atoms[idx].dim_factor = std::min(atoms[idx].dim_factor, 0.35F);
    }
  }
}

namespace {

void markShellAroundAtom(const std::vector<render3d::AtomInstance> &atoms, std::size_t ref_idx,
                         float r_lo_sq, float r_hi_sq, std::vector<bool> &in_shell) {
  in_shell[ref_idx] = true;
  const auto &ref_pos = atoms[ref_idx].position;
  for (std::size_t idx = 0; idx < atoms.size(); ++idx) {
    if (idx == ref_idx) {
      continue;
    }
    const auto diff = atoms[idx].position - ref_pos;
    const float dist_sq = math::dot(diff, diff);
    if (dist_sq >= r_lo_sq && dist_sq <= r_hi_sq) {
      in_shell[idx] = true;
    }
  }
}

void markAllShellPairs(const std::vector<render3d::AtomInstance> &atoms, float r_lo_sq,
                       float r_hi_sq, std::vector<bool> &in_shell) {
  const std::size_t limit = std::min<std::size_t>(atoms.size(), 2500);
  for (std::size_t outer = 0; outer < limit; ++outer) {
    for (std::size_t inner = outer + 1; inner < limit; ++inner) {
      const auto diff = atoms[inner].position - atoms[outer].position;
      const float dist_sq = math::dot(diff, diff);
      if (dist_sq >= r_lo_sq && dist_sq <= r_hi_sq) {
        in_shell[outer] = true;
        in_shell[inner] = true;
      }
    }
  }
}

struct NeighborStats {
  int count{0};
  float sum_distance{0.0F};
};

[[nodiscard]] NeighborStats computeFirstShellStats(const std::vector<render3d::AtomInstance> &atoms,
                                                   std::size_t center_idx, float cutoff_scale) {
  NeighborStats stats{};
  const auto &center = atoms[center_idx];
  for (std::size_t other = 0; other < atoms.size(); ++other) {
    if (other == center_idx) {
      continue;
    }
    const float cutoff = (center.covalent_radius + atoms[other].covalent_radius) * cutoff_scale;
    const float dist = math::norm(atoms[other].position - center.position);
    if (dist > 0.01F && dist <= cutoff) {
      ++stats.count;
      stats.sum_distance += dist;
    }
  }
  return stats;
}

} // namespace

void Structure3DController::applyRadialBrushHighlightLocked(std::size_t &matched_atoms) {
  matched_atoms = 0;
  if (scene_ == nullptr || !radial_brush_r_.has_value()) {
    return;
  }
  auto &atoms = scene_->atoms();
  if (atoms.empty()) {
    return;
  }
  constexpr float SHELL_HALF_WIDTH = 0.15F;
  const auto target_r = static_cast<float>(*radial_brush_r_);
  const float r_lo = std::max(0.0F, target_r - SHELL_HALF_WIDTH);
  const float r_hi = target_r + SHELL_HALF_WIDTH;
  const float r_lo_sq = r_lo * r_lo;
  const float r_hi_sq = r_hi * r_hi;

  std::vector<bool> in_shell(atoms.size(), false);
  if (picked_atom_index_.has_value() && *picked_atom_index_ < atoms.size()) {
    markShellAroundAtom(atoms, *picked_atom_index_, r_lo_sq, r_hi_sq, in_shell);
  } else {
    markAllShellPairs(atoms, r_lo_sq, r_hi_sq, in_shell);
  }

  for (std::size_t idx = 0; idx < atoms.size(); ++idx) {
    if (in_shell[idx]) {
      ++matched_atoms;
      atoms[idx].halo_color = {.r = 86, .g = 180, .b = 233, .a = 215};
      atoms[idx].halo_radius_ratio = 2.2F;
      atoms[idx].dim_factor = 1.0F;
    } else {
      atoms[idx].dim_factor = std::min(atoms[idx].dim_factor, 0.25F);
    }
  }
}

void Structure3DController::applySlabClippingLocked() {
  if (scene_ == nullptr || !slab_clip_.enabled) {
    return;
  }
  const auto &min_b = scene_->boundsMin();
  const auto &max_b = scene_->boundsMax();
  const float span_x = std::max(max_b.x() - min_b.x(), 1e-4F);
  const float span_y = std::max(max_b.y() - min_b.y(), 1e-4F);
  const float span_z = std::max(max_b.z() - min_b.z(), 1e-4F);

  const float lo_x = min_b.x() + (slab_clip_.x_min * span_x);
  const float hi_x = min_b.x() + (slab_clip_.x_max * span_x);
  const float lo_y = min_b.y() + (slab_clip_.y_min * span_y);
  const float hi_y = min_b.y() + (slab_clip_.y_max * span_y);
  const float lo_z = min_b.z() + (slab_clip_.z_min * span_z);
  const float hi_z = min_b.z() + (slab_clip_.z_max * span_z);

  for (auto &atom : scene_->atoms()) {
    const auto &pos = atom.position;
    if (pos.x() < lo_x || pos.x() > hi_x || pos.y() < lo_y || pos.y() > hi_y || pos.z() < lo_z ||
        pos.z() > hi_z) {
      atom.dim_factor = 0.0F;
      atom.halo_color.a = 0;
    }
  }
}

void Structure3DController::syncInspectedAtomToWindow() {
  InspectedAtomInfo info{};
  info.valid = false;
  {
    const std::scoped_lock lock(mutex_);
    if (scene_ == nullptr || !picked_atom_index_.has_value() ||
        *picked_atom_index_ >= scene_->atoms().size()) {
      safeInvokeOnWindow([this, info] { window_.set_structure_3d_inspected_atom(info); });
      return;
    }

    const std::size_t idx = *picked_atom_index_;
    const auto &atoms = scene_->atoms();
    const auto &atom = atoms[idx];
    info.valid = true;
    info.atom_id = static_cast<int>(atom.id + 1);
    info.element = slint::SharedString(
        idx < current_cell_.atoms().size() ? current_cell_.atoms()[idx].element().symbol : "X");
    info.coords_text = slint::SharedString(std::format("{:.2f}, {:.2f}, {:.2f}", atom.position.x(),
                                                       atom.position.y(), atom.position.z()));

    const auto stats = computeFirstShellStats(atoms, idx, static_cast<float>(bond_cutoff_scale_));
    info.cn_text = slint::SharedString(std::to_string(stats.count));
    if (stats.count > 0) {
      info.avg_bond_text = slint::SharedString(
          std::format("{:.3f} Å", stats.sum_distance / static_cast<float>(stats.count)));
    }

    if (phase_service_ != nullptr) {
      const auto active_frame = streamer_.hasTrajectory() ? streamer_.currentFrameIndex() : 0UL;
      const auto desc = phase_service_->getOrComputeDescriptors(active_frame);
      if (desc != nullptr && desc->steinhardt != nullptr) {
        if (idx < desc->steinhardt->Q6_bar.size()) {
          info.q6_text = slint::SharedString(std::format("{:.3f}", desc->steinhardt->Q6_bar[idx]));
        }
        if (idx < desc->steinhardt->Q4_bar.size()) {
          info.q4_text = slint::SharedString(std::format("{:.3f}", desc->steinhardt->Q4_bar[idx]));
        }
      }
    }
  }
  safeInvokeOnWindow([this, info] { window_.set_structure_3d_inspected_atom(info); });
}

void Structure3DController::setDistributionStore(
    std::shared_ptr<const analysis::PerFrameDistributionStore> store) {
  {
    const std::scoped_lock lock(mutex_);
    distribution_store_ = std::move(store);
  }
  syncDistributionHUDToWindow();
}

void Structure3DController::handleSelectDistributionMetric(int metric_index) {
  if (metric_index < 0 || metric_index > 3) {
    return;
  }
  {
    const std::scoped_lock lock(mutex_);
    active_dist_metric_ = static_cast<analysis::DistributionMetric>(metric_index);
  }
  syncDistributionHUDToWindow();
}

void Structure3DController::handleToggleDistributionHUD(bool visible) {
  {
    const std::scoped_lock lock(mutex_);
    dist_hud_visible_ = visible;
  }
  syncDistributionHUDToWindow();
}

void Structure3DController::setMetricPlotLimits(analysis::DistributionMetric metric,
                                                const render3d::MiniPlotConfig &config) {
  const auto idx = static_cast<std::size_t>(metric);
  if (idx < metric_plot_configs_.size()) {
    {
      const std::scoped_lock lock(mutex_);
      metric_plot_configs_.at(idx) = config;
    }
    syncDistributionHUDToWindow();
  }
}

void Structure3DController::syncDistributionHUDToWindow() {
  safeInvokeOnWindow([this] {
    std::shared_ptr<const analysis::PerFrameDistributionStore> store;
    analysis::DistributionMetric metric = analysis::DistributionMetric::RDF;
    render3d::MiniPlotConfig plot_cfg{};
    bool visible = false;
    std::size_t active_frame = 0;

    {
      const std::scoped_lock lock(mutex_);
      store = distribution_store_;
      metric = active_dist_metric_;
      visible = dist_hud_visible_;
      const auto idx = static_cast<std::size_t>(metric);
      if (idx < metric_plot_configs_.size()) {
        plot_cfg = metric_plot_configs_.at(idx);
      }
      active_frame = streamer_.hasTrajectory() ? streamer_.currentFrameIndex() : 0UL;
    }

    window_.set_structure_3d_dist_hud_visible(visible);

    if (!store || !visible || !store->hasData(metric)) {
      window_.set_structure_3d_dist_has_data(false);
      return;
    }

    const auto &metric_data = store->getMetricData(metric);
    const auto buf =
        render3d::MiniPlotRenderer::render(metric_data, active_frame, plot_cfg, 260, 126);

    const slint::SharedPixelBuffer<slint::Rgba8Pixel> pixel_buffer(
        buf.width, buf.height, reinterpret_cast<const slint::Rgba8Pixel *>(buf.pixels.data()));
    const slint::Image chart_img(pixel_buffer);

    window_.set_structure_3d_dist_image(chart_img);
    window_.set_structure_3d_dist_metric(static_cast<int>(metric));
    window_.set_structure_3d_dist_hud_visible(visible);
    window_.set_structure_3d_dist_has_data(true);

    const auto peak = render3d::MiniPlotRenderer::findPeak(
        store->getBins(metric), store->getFrameCurve(metric, active_frame));
    if (peak.valid) {
      window_.set_structure_3d_dist_peak_text(
          slint::SharedString(std::format("Peak: {:.2f} {}", peak.x, metric_data.x_unit)));
    } else {
      window_.set_structure_3d_dist_peak_text(slint::SharedString(""));
    }
  });
}

} // namespace correlation::app
