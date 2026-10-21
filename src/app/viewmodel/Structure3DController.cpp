/**
 * @file Structure3DController.cpp
 * @brief Implementation of Structure3DController for interactive Slint 3D visualization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/Structure3DController.hpp"
#include "AppWindow.h"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/order/SteinhardtCalculator.hpp"
#include "physics/PhysicalData.hpp"
#include "plotters/ColormapService.hpp"
#include "render3d/RenderEngineFactory.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <utility>

namespace correlation::app {

namespace {

real_t getElementCovalentRadius(const std::string &symbol) {
  try {
    return physics::getCovalentRadius(symbol);
  } catch (const std::out_of_range &) {
    return static_cast<real_t>(1.5);
  }
}

analysis::BondCutoffMatrix buildCellBondCutoffs(const core::Cell &cell) {
  const auto &elements = cell.elements();
  const std::size_t num_elements = elements.size();
  analysis::BondCutoffMatrix matrix(num_elements);

  for (std::size_t i = 0; i < num_elements; ++i) {
    const real_t radius_a = getElementCovalentRadius(elements[i].symbol);
    for (std::size_t j = i; j < num_elements; ++j) {
      const real_t radius_b = getElementCovalentRadius(elements[j].symbol);
      const real_t sum_radii = radius_a + radius_b;
      const real_t min_bond_dist = sum_radii * static_cast<real_t>(0.6);
      const real_t max_bond_dist = sum_radii * static_cast<real_t>(1.3);
      matrix[i][j] = analysis::BondCutoffRange{
          .min_sq = min_bond_dist * min_bond_dist,
          .max_sq = max_bond_dist * max_bond_dist,
      };
      matrix[j][i] = matrix[i][j];
    }
  }
  return matrix;
}

real_t computeMaxBondDistance(const analysis::BondCutoffMatrix &matrix, std::size_t num_elements) {
  real_t max_dist = 0.0;
  for (std::size_t i = 0; i < num_elements; ++i) {
    for (std::size_t j = i; j < num_elements; ++j) {
      max_dist = std::max(max_dist, static_cast<real_t>(std::sqrt(matrix[i][j].max_sq)));
    }
  }
  return max_dist;
}

std::span<const real_t>
extractDescriptorValues(const calculators::SteinhardtCalculator::SteinhardtParams &params,
                        PhaseDescriptor descriptor) noexcept {
  switch (descriptor) {
  case PhaseDescriptor::Q6Bar:
    return params.Q6_bar;
  case PhaseDescriptor::Q4Bar:
    return params.Q4_bar;
  case PhaseDescriptor::Q6:
    return params.Q6;
  case PhaseDescriptor::Q4:
    return params.Q4;
  case PhaseDescriptor::W6Hat:
    return params.W6_hat;
  }
  return {};
}

} // namespace

Structure3DController::Structure3DController(::AppWindow &window, const ProgramOptions &options)
    : window_(window), options_(options),
      engine_(render3d::RenderEngineFactory::createEngine(render3d::BackendPreference::Auto)),
      alive_(std::make_shared<std::atomic<bool>>(true)) {
  if (engine_ != nullptr) {
    static_cast<void>(engine_->initialize(viewport_dims_));
  }
  camera_.setViewport(static_cast<float>(viewport_dims_.width),
                      static_cast<float>(viewport_dims_.height));
  camera_.fitToBox(math::Vector3<float>{-10.0F, -10.0F, -10.0F},
                   math::Vector3<float>{10.0F, 10.0F, 10.0F});

  syncPhaseSearchToWindow();

  phase_worker_ = std::jthread([this] { phaseWorkerLoop(); });
  render_worker_ = std::jthread([this] { renderLoop(); });
}

Structure3DController::~Structure3DController() {
  if (alive_ != nullptr) {
    alive_->store(false);
  }
  stop_phase_worker_.store(true);
  phase_cv_.notify_all();
  if (phase_worker_.joinable()) {
    phase_worker_.join();
  }

  {
    const std::scoped_lock lock(mutex_);
    stop_worker_.store(true);
  }
  cv_.notify_all();
  if (render_worker_.joinable()) {
    render_worker_.join();
  }
}

void Structure3DController::requestRender(bool /*immediate*/) {
  {
    const std::scoped_lock lock(mutex_);
    render_pending_ = true;
  }
  cv_.notify_one();
}

void Structure3DController::handleViewportResized(float width, float height) {
  if (width < 10.0F || height < 10.0F) {
    return;
  }
  const auto new_width = static_cast<uint32_t>(std::round(width));
  const auto new_height = static_cast<uint32_t>(std::round(height));

  {
    const std::scoped_lock lock(mutex_);
    if (viewport_dims_.width == new_width && viewport_dims_.height == new_height) {
      return;
    }
    viewport_dims_ = {.width = new_width, .height = new_height};
    camera_.setViewport(static_cast<float>(new_width), static_cast<float>(new_height));
    if (engine_ != nullptr) {
      engine_->resize(viewport_dims_);
    }
    render_pending_ = true;
  }
  cv_.notify_one();
}

void Structure3DController::handleRotate(float delta_x, float delta_y) {
  {
    const std::scoped_lock lock(mutex_);
    camera_.orbit(delta_x, delta_y);
    render_pending_ = true;
  }
  cv_.notify_one();
}

void Structure3DController::handlePan(float delta_x, float delta_y) {
  {
    const std::scoped_lock lock(mutex_);
    camera_.pan(delta_x, delta_y);
    render_pending_ = true;
  }
  cv_.notify_one();
}

void Structure3DController::handleZoom(float delta) {
  {
    const std::scoped_lock lock(mutex_);
    camera_.zoom(delta);
    render_pending_ = true;
  }
  cv_.notify_one();
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
    render_pending_ = true;
  }
  cv_.notify_one();
}

void Structure3DController::handleToggleBonds(bool enabled) {
  {
    const std::scoped_lock lock(mutex_);
    render_bonds_ = enabled;
    if (render_bonds_ && scene_ != nullptr) {
      scene_->recomputeBonds(bond_cutoff_scale_, in_cell_mode_);
    }
    render_pending_ = true;
  }
  window_.set_structure_3d_bonds_enabled(enabled);
  cv_.notify_one();
}

void Structure3DController::handleToggleWireframe(bool enabled) {
  {
    const std::scoped_lock lock(mutex_);
    render_wireframe_ = enabled;
    render_pending_ = true;
  }
  window_.set_structure_3d_wireframe_enabled(enabled);
  cv_.notify_one();
}

void Structure3DController::handleToggleInCell(bool enabled) {
  {
    const std::scoped_lock lock(phase_mutex_, mutex_);
    in_cell_mode_ = enabled;
    if (scene_ != nullptr) {
      if (streamer_.hasTrajectory() && streamer_.totalFrames() > 0) {
        static_cast<void>(streamer_.seekFrame(streamer_.currentFrameIndex(), *scene_, render_bonds_,
                                              bond_cutoff_scale_, in_cell_mode_));
      } else {
        scene_ = render3d::Scene3D::fromCell(current_cell_, bond_cutoff_scale_, render_bonds_,
                                             in_cell_mode_);
        if (engine_ != nullptr) {
          engine_->setScene(scene_);
        }
      }
      updateScenePhaseVisuals();
      render_pending_ = true;
    }
  }
  window_.set_structure_3d_in_cell(enabled);
  cv_.notify_one();
}

void Structure3DController::handleSeekFrame(int frame_index) {
  if (frame_index < 0) {
    return;
  }
  int actual_frame = frame_index;
  {
    const std::scoped_lock lock(phase_mutex_, mutex_);
    if (streamer_.totalFrames() > 0 && scene_ != nullptr) {
      const auto target = std::min(static_cast<size_t>(frame_index), streamer_.totalFrames() - 1);
      static_cast<void>(
          streamer_.seekFrame(target, *scene_, render_bonds_, bond_cutoff_scale_, in_cell_mode_));
      actual_frame = static_cast<int>(streamer_.currentFrameIndex());
      updateScenePhaseVisuals();
      render_pending_ = true;
    }
  }
  window_.set_structure_3d_current_frame(actual_frame);
  cv_.notify_one();
}

void Structure3DController::handleStepFrame(bool forward) {
  int current_idx = 0;
  {
    const std::scoped_lock lock(phase_mutex_, mutex_);
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
      render_pending_ = true;
    }
  }
  window_.set_structure_3d_current_frame(current_idx);
  cv_.notify_one();
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
      render_pending_ = true;
    }
  }
  cv_.notify_one();
}

void Structure3DController::setCell(const core::Cell &cell) {
  current_cell_ = cell;
  auto new_scene =
      render3d::Scene3D::fromCell(cell, bond_cutoff_scale_, render_bonds_, in_cell_mode_);
  {
    const std::scoped_lock lock(phase_mutex_, mutex_);
    current_trajectory_.reset();
    lru_cache_.clear();
    scene_ = std::move(new_scene);
    if (scene_ != nullptr) {
      camera_.fitToBox(scene_->boundsMin(), scene_->boundsMax());
      if (engine_ != nullptr) {
        engine_->setScene(scene_);
      }
      updateScenePhaseVisuals();
    }
    render_pending_ = true;
  }
  window_.set_structure_3d_current_frame(0);
  window_.set_structure_3d_total_frames(1);
  cv_.notify_one();
  syncPhaseSearchToWindow();
}

void Structure3DController::setTrajectory(std::shared_ptr<const core::Trajectory> trajectory) {
  if (!trajectory) {
    return;
  }
  const auto total_frames = static_cast<int>(trajectory->getFrameCount());
  {
    const std::scoped_lock lock(phase_mutex_, mutex_);
    current_trajectory_ = trajectory;
    lru_cache_.clear();
    if (total_frames > 0) {
      const auto &first_cell = trajectory->firstFrame();
      current_cell_ = first_cell;
      scene_ =
          render3d::Scene3D::fromCell(first_cell, bond_cutoff_scale_, render_bonds_, in_cell_mode_);
      if (scene_ != nullptr) {
        camera_.fitToBox(scene_->boundsMin(), scene_->boundsMax());
        if (engine_ != nullptr) {
          engine_->setScene(scene_);
        }
      }
    }
    streamer_.setTrajectory(std::move(trajectory));
    updateScenePhaseVisuals();
    render_pending_ = true;
  }
  window_.set_structure_3d_current_frame(0);
  window_.set_structure_3d_total_frames(total_frames);
  cv_.notify_one();
  syncPhaseSearchToWindow();
}

void Structure3DController::set3DModeActive(bool active) {
  {
    const std::scoped_lock lock(mutex_);
    if (is_active_ == active) {
      return;
    }
    is_active_ = active;
    if (is_active_) {
      render_pending_ = true;
    }
  }
  if (active) {
    cv_.notify_one();
  }
}

void Structure3DController::setPhaseSearchEnabled(bool enabled) {
  {
    const std::scoped_lock lock(phase_mutex_, mutex_);
    phase_settings_.enabled = enabled;
    updateScenePhaseVisuals();
  }
  cv_.notify_one();
  syncPhaseSearchToWindow();
}

bool Structure3DController::isPhaseSearchEnabled() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.enabled;
}

void Structure3DController::setPhaseDescriptor(PhaseDescriptor descriptor) {
  {
    const std::scoped_lock lock(phase_mutex_, mutex_);
    phase_settings_.descriptor = descriptor;
    phase_settings_.preset = CrystalPreset::Custom;
    updateScenePhaseVisuals();
  }
  cv_.notify_one();
  syncPhaseSearchToWindow();
}

PhaseDescriptor Structure3DController::phaseDescriptor() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.descriptor;
}

void Structure3DController::setCrystalPreset(CrystalPreset preset) {
  {
    const std::scoped_lock lock(phase_mutex_, mutex_);
    phase_settings_.preset = preset;
    if (preset != CrystalPreset::Custom) {
      const auto config = getPresetConfig(preset);
      phase_settings_.descriptor = config.descriptor;
      phase_settings_.min_threshold = config.min_threshold;
      phase_settings_.max_threshold = config.max_threshold;
    }
    updateScenePhaseVisuals();
  }
  cv_.notify_one();
  syncPhaseSearchToWindow();
}

CrystalPreset Structure3DController::crystalPreset() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.preset;
}

void Structure3DController::setPhaseThresholds(float min_threshold, float max_threshold) {
  {
    const std::scoped_lock lock(phase_mutex_, mutex_);
    phase_settings_.preset = CrystalPreset::Custom;
    phase_settings_.min_threshold = std::min(min_threshold, max_threshold);
    phase_settings_.max_threshold = std::max(min_threshold, max_threshold);
    updateScenePhaseVisuals();
  }
  cv_.notify_one();
  syncPhaseSearchToWindow();
}

float Structure3DController::phaseMinThreshold() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.min_threshold;
}

float Structure3DController::phaseMaxThreshold() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.max_threshold;
}

void Structure3DController::setPhaseColormap(plotters::ColormapType colormap) {
  {
    const std::scoped_lock lock(phase_mutex_, mutex_);
    phase_settings_.colormap = colormap;
    updateScenePhaseVisuals();
  }
  cv_.notify_one();
  syncPhaseSearchToWindow();
}

plotters::ColormapType Structure3DController::phaseColormap() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.colormap;
}

void Structure3DController::setPhaseHaloIntensity(float intensity) {
  {
    const std::scoped_lock lock(phase_mutex_, mutex_);
    phase_settings_.halo_intensity = std::clamp(intensity, 0.1F, 1.0F);
    updateScenePhaseVisuals();
  }
  cv_.notify_one();
  syncPhaseSearchToWindow();
}

float Structure3DController::phaseHaloIntensity() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.halo_intensity;
}

void Structure3DController::setPhaseHaloRadiusRatio(float ratio) {
  {
    const std::scoped_lock lock(phase_mutex_, mutex_);
    phase_settings_.halo_radius_ratio = std::clamp(ratio, 1.0F, 3.0F);
    updateScenePhaseVisuals();
  }
  cv_.notify_one();
  syncPhaseSearchToWindow();
}

float Structure3DController::phaseHaloRadiusRatio() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.halo_radius_ratio;
}

void Structure3DController::setPhaseDimNonMatching(bool dim) {
  {
    const std::scoped_lock lock(phase_mutex_, mutex_);
    phase_settings_.dim_non_matching = dim;
    updateScenePhaseVisuals();
  }
  cv_.notify_one();
  syncPhaseSearchToWindow();
}

bool Structure3DController::isPhaseDimNonMatching() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.dim_non_matching;
}

PhaseSearchSettings Structure3DController::phaseSearchSettings() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_;
}

PhaseHistogram Structure3DController::phaseHistogram() const {
  const std::scoped_lock lock(phase_mutex_);
  return latest_histogram_;
}

void Structure3DController::syncPhaseSearchToWindow() {
  PhaseSearchSettings settings_copy;
  PhaseHistogram hist_copy;
  {
    const std::scoped_lock lock(phase_mutex_);
    settings_copy = phase_settings_;
    hist_copy = latest_histogram_;
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
    uint32_t max_count = 0;
    for (const auto c : hist_copy.counts) {
      max_count = std::max(max_count, c);
    }
    const float norm_factor = (max_count > 0) ? (1.0F / static_cast<float>(max_count)) : 0.0F;
    for (const auto c : hist_copy.counts) {
      counts_model->push_back(static_cast<float>(c) * norm_factor);
    }
  }

  auto apply_to_window = [this, settings = settings_copy, hist = hist_copy,
                          colorbar_img = std::move(colorbar_img),
                          counts_model = std::move(counts_model)] mutable {
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
    window_.set_structure_3d_phase_dim_non_matching(settings.dim_non_matching);
    window_.set_structure_3d_phase_matching_atoms(static_cast<int>(hist.matching_atom_count));
    window_.set_structure_3d_phase_total_atoms(static_cast<int>(hist.total_atom_count));
    window_.set_structure_3d_phase_descriptor_name(
        slint::SharedString(descriptorToString(settings.descriptor)));
    window_.set_structure_3d_phase_colorbar_image(colorbar_img);
    window_.set_structure_3d_phase_histogram_counts(counts_model);
  };

  if (std::this_thread::get_id() == main_thread_id_) {
    apply_to_window();
  } else {
    slint::invoke_from_event_loop(
        [alive = alive_, sync_action = std::move(apply_to_window)] mutable {
          if (!alive->load()) {
            return;
          }
          sync_action();
        });
  }
}

std::shared_ptr<const calculators::SteinhardtCalculator::SteinhardtParams>
Structure3DController::getOrComputeSteinhardt(std::size_t frame_index) {
  auto cached = lru_cache_.get(frame_index);
  if (cached) {
    return cached;
  }
  auto computed = computeSteinhardtSync(frame_index);
  if (computed) {
    lru_cache_.put(frame_index, computed);
  }
  return computed;
}

void Structure3DController::requestPhaseCalculation(std::size_t frame_index) {
  requested_phase_frame_.store(static_cast<int64_t>(frame_index));
  phase_cv_.notify_one();
}

void Structure3DController::updateScenePhaseVisuals() {
  if (scene_ == nullptr) {
    return;
  }
  const auto active_frame = streamer_.hasTrajectory() ? streamer_.currentFrameIndex() : 0UL;
  const auto params = lru_cache_.peek(active_frame);
  if (phase_settings_.enabled && !params) {
    requestPhaseCalculation(active_frame);
  }
  applyPhaseVisualsToScene(*scene_, params.get());
  latest_histogram_ = computePhaseHistogram(params.get());
  render_pending_ = true;
}

void Structure3DController::applyPhaseVisualsToScene(
    render3d::Scene3D &target_scene,
    const calculators::SteinhardtCalculator::SteinhardtParams *params) const {
  auto &atoms = target_scene.atoms();
  if (!phase_settings_.enabled || params == nullptr) {
    for (auto &atom : atoms) {
      atom.halo_color = render3d::RgbaColor{.r = 0, .g = 0, .b = 0, .a = 0};
      atom.dim_factor = 1.0F;
    }
    return;
  }

  const auto values = extractDescriptorValues(*params, phase_settings_.descriptor);
  if (values.size() != atoms.size()) {
    return;
  }

  const float min_thresh = phase_settings_.min_threshold;
  const float max_thresh = phase_settings_.max_threshold;
  const float diff = max_thresh - min_thresh;
  const float halo_alpha = std::clamp(phase_settings_.halo_intensity, 0.0F, 1.0F);
  const auto alpha_byte = static_cast<uint8_t>(std::round(halo_alpha * 255.0F));

  for (std::size_t i = 0; i < atoms.size(); ++i) {
    const auto val = static_cast<float>(values[i]);
    if (val >= min_thresh && val <= max_thresh) {
      const double norm_t =
          (diff > 1e-6F) ? std::clamp(static_cast<double>((val - min_thresh) / diff), 0.0, 1.0)
                         : 0.5;
      const auto sampled = plotters::ColormapService::sample(norm_t, phase_settings_.colormap);
      atoms[i].halo_color = render3d::RgbaColor{
          .r = sampled.r,
          .g = sampled.g,
          .b = sampled.b,
          .a = alpha_byte,
      };
      atoms[i].halo_radius_ratio = phase_settings_.halo_radius_ratio;
      atoms[i].dim_factor = 1.0F;
    } else {
      atoms[i].halo_color = render3d::RgbaColor{.r = 0, .g = 0, .b = 0, .a = 0};
      atoms[i].dim_factor = phase_settings_.dim_non_matching ? phase_settings_.dim_factor : 1.0F;
    }
  }
}

PhaseHistogram Structure3DController::computePhaseHistogram(
    const calculators::SteinhardtCalculator::SteinhardtParams *params) const {
  PhaseHistogram hist;
  constexpr std::size_t K_BIN_COUNT = 50;
  hist.counts.assign(K_BIN_COUNT, 0);
  hist.bin_centers.resize(K_BIN_COUNT);

  const auto [domain_min, domain_max] = getDescriptorDomain(phase_settings_.descriptor);
  hist.min_val = domain_min;
  hist.max_val = domain_max;
  hist.bin_width = (domain_max - domain_min) / static_cast<float>(K_BIN_COUNT);

  for (std::size_t i = 0; i < K_BIN_COUNT; ++i) {
    hist.bin_centers[i] = domain_min + (static_cast<float>(i) + 0.5F) * hist.bin_width;
  }

  if (params == nullptr) {
    return hist;
  }

  const auto values = extractDescriptorValues(*params, phase_settings_.descriptor);
  hist.total_atom_count = static_cast<uint32_t>(values.size());

  const float min_thresh = phase_settings_.min_threshold;
  const float max_thresh = phase_settings_.max_threshold;

  for (const auto val_raw : values) {
    const auto val = static_cast<float>(val_raw);
    if (val >= min_thresh && val <= max_thresh) {
      ++hist.matching_atom_count;
    }
    const auto bin_idx =
        static_cast<std::size_t>(std::clamp(std::floor((val - domain_min) / hist.bin_width), 0.0F,
                                            static_cast<float>(K_BIN_COUNT - 1)));
    ++hist.counts[bin_idx];
  }

  return hist;
}

std::shared_ptr<const calculators::SteinhardtCalculator::SteinhardtParams>
Structure3DController::computeSteinhardtSync(std::size_t frame_index) {
  core::Cell target_cell;
  analysis::BondCutoffMatrix cutoffs;
  {
    const std::scoped_lock lock(mutex_);
    if (current_trajectory_ && streamer_.totalFrames() > 0) {
      if (frame_index >= current_trajectory_->getFrameCount()) {
        return nullptr;
      }
      target_cell = current_trajectory_->getFrame(frame_index);
      if (!current_trajectory_->getBondCutoffs().empty()) {
        cutoffs = current_trajectory_->getBondCutoffs();
      }
    } else {
      if (frame_index != 0 || current_cell_.isEmpty()) {
        return nullptr;
      }
      target_cell = current_cell_;
    }
  }

  if (target_cell.isEmpty() || target_cell.atomCount() == 0) {
    return nullptr;
  }

  if (cutoffs.empty() || cutoffs.size() != target_cell.elements().size()) {
    cutoffs = buildCellBondCutoffs(target_cell);
  }

  const real_t max_dist = computeMaxBondDistance(cutoffs, target_cell.elements().size());
  try {
    const analysis::StructureAnalyzer analyzer(target_cell, max_dist, cutoffs, false);
    auto params = calculators::SteinhardtCalculator::calculatePerAtom(target_cell, &analyzer);
    return std::make_shared<const calculators::SteinhardtCalculator::SteinhardtParams>(
        std::move(params));
  } catch (const std::exception &) {
    return nullptr;
  }
}

int64_t Structure3DController::waitForNextPhaseFrame() {
  std::unique_lock lock(phase_mutex_);
  phase_cv_.wait_for(lock, std::chrono::milliseconds(50), [this] {
    return stop_phase_worker_.load(std::memory_order_relaxed) ||
           (requested_phase_frame_.load(std::memory_order_relaxed) >= 0);
  });

  if (stop_phase_worker_.load(std::memory_order_relaxed)) {
    return -1;
  }

  return requested_phase_frame_.exchange(-1);
}

void Structure3DController::processTargetPhaseFrame(int64_t target_frame) {
  const auto frame_idx = static_cast<std::size_t>(target_frame);
  auto params = lru_cache_.get(frame_idx);
  if (!params) {
    params = computeSteinhardtSync(frame_idx);
    if (params) {
      lru_cache_.put(frame_idx, params);
    }
  }

  if (!params) {
    return;
  }

  bool need_render = false;
  {
    const std::scoped_lock lock(phase_mutex_, mutex_);
    const auto active_frame = streamer_.hasTrajectory() ? streamer_.currentFrameIndex() : 0UL;
    if (active_frame == frame_idx && scene_ != nullptr) {
      applyPhaseVisualsToScene(*scene_, params.get());
      latest_histogram_ = computePhaseHistogram(params.get());
      render_pending_ = true;
      need_render = true;
    }
  }
  if (need_render) {
    cv_.notify_one();
    syncPhaseSearchToWindow();
  }
}

void Structure3DController::prewarmAdjacentFrames(int64_t target_frame) {
  const std::size_t total_frames = streamer_.hasTrajectory() ? streamer_.totalFrames() : 1UL;
  if (total_frames <= 1) {
    return;
  }

  constexpr std::array<int64_t, 4> K_NEIGHBOR_OFFSETS = {-1, 1, -2, 2};
  for (const auto offset : K_NEIGHBOR_OFFSETS) {
    if (stop_phase_worker_.load(std::memory_order_relaxed) ||
        requested_phase_frame_.load(std::memory_order_relaxed) >= 0) {
      break;
    }
    const int64_t candidate = target_frame + offset;
    if (std::cmp_greater_equal(candidate, 0) && std::cmp_less(candidate, total_frames)) {
      const auto cand_idx = static_cast<std::size_t>(candidate);
      if (!lru_cache_.contains(cand_idx)) {
        auto pre_params = computeSteinhardtSync(cand_idx);
        if (pre_params) {
          lru_cache_.put(cand_idx, std::move(pre_params));
        }
      }
    }
  }
}

void Structure3DController::phaseWorkerLoop() {
  while (!stop_phase_worker_.load(std::memory_order_relaxed)) {
    const int64_t target_frame = waitForNextPhaseFrame();
    if (target_frame < 0) {
      continue;
    }

    processTargetPhaseFrame(target_frame);
    prewarmAdjacentFrames(target_frame);
  }
}

void Structure3DController::dispatchFrame(const uint8_t *pixel_data, uint32_t width,
                                          uint32_t height) {
  if (pixel_data == nullptr || width == 0 || height == 0) {
    return;
  }

  const auto *rgba_pixels = reinterpret_cast<const slint::Rgba8Pixel *>(pixel_data);
  const slint::SharedPixelBuffer<slint::Rgba8Pixel> buffer(width, height, rgba_pixels);
  slint::Image img(buffer);

  slint::invoke_from_event_loop([this, alive = alive_, img = std::move(img)] mutable {
    if (!alive->load()) {
      return;
    }
    window_.set_structure_3d_image(img);
  });
}

void Structure3DController::renderLoop() {
  while (!stop_worker_.load(std::memory_order_relaxed)) {
    render3d::Camera3D cam_copy;
    std::shared_ptr<render3d::Scene3D> scene_copy;
    bool active = false;
    bool bonds_vis = true;
    bool wireframe_vis = true;
    bool interacting = false;

    {
      std::unique_lock lock(mutex_);
      cv_.wait_for(lock, std::chrono::milliseconds(50), [this] {
        return stop_worker_.load(std::memory_order_relaxed) || (render_pending_ && is_active_);
      });

      if (stop_worker_.load(std::memory_order_relaxed)) {
        break;
      }

      render_pending_ = false;
      active = is_active_;
      cam_copy = camera_;
      scene_copy = scene_;
      bonds_vis = render_bonds_;
      wireframe_vis = render_wireframe_;
      interacting = is_interacting_;
    }

    if (!active || scene_copy == nullptr || engine_ == nullptr) {
      continue;
    }

    const auto lod_policy = lod_controller_.evaluatePolicy(scene_copy->atoms().size(), interacting);
    const bool show_bonds = bonds_vis && lod_policy.render_bonds;

    engine_->setFeatureVisibility(show_bonds, wireframe_vis);
    engine_->setCamera(cam_copy);
    engine_->setScene(scene_copy);

    const auto pixel_span = engine_->renderFrame();
    if (!pixel_span.empty()) {
      const auto dims = engine_->dimensions();
      dispatchFrame(pixel_span.data(), dims.width, dims.height);
    }
  }
}

} // namespace correlation::app
