/**
 * @file PhaseSearchService.cpp
 * @brief Implementation of PhaseSearchService for order parameter calculations and phase visuals.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/PhaseSearchService.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/order/LocalEntropyCalculator.hpp"
#include "calculators/order/SteinhardtCalculator.hpp"
#include "calculators/order/TenWoldeFrenkelCalculator.hpp"
#include "calculators/order/TetrahedralCalculator.hpp"
#include "physics/PhysicalData.hpp"
#include "plotters/ColormapService.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
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

std::span<const real_t> extractDescriptorValues(const FramePhaseDescriptors &descriptors,
                                                PhaseDescriptor descriptor) noexcept {
  switch (descriptor) {
  case PhaseDescriptor::Q6Bar:
    return descriptors.steinhardt ? descriptors.steinhardt->Q6_bar : std::span<const real_t>{};
  case PhaseDescriptor::Q4Bar:
    return descriptors.steinhardt ? descriptors.steinhardt->Q4_bar : std::span<const real_t>{};
  case PhaseDescriptor::Q6:
    return descriptors.steinhardt ? descriptors.steinhardt->Q6 : std::span<const real_t>{};
  case PhaseDescriptor::Q4:
    return descriptors.steinhardt ? descriptors.steinhardt->Q4 : std::span<const real_t>{};
  case PhaseDescriptor::W6Hat:
    return descriptors.steinhardt ? descriptors.steinhardt->W6_hat : std::span<const real_t>{};
  case PhaseDescriptor::W4Hat:
    return descriptors.steinhardt ? descriptors.steinhardt->W4_hat : std::span<const real_t>{};
  case PhaseDescriptor::LocalEntropy:
    return descriptors.local_entropy;
  case PhaseDescriptor::TwfNconn:
    return descriptors.twf_nconn;
  case PhaseDescriptor::QTet:
    return descriptors.q_tet;
  }
  return {};
}

void resetPhaseVisuals(std::span<render3d::AtomInstance> atoms) noexcept {
  for (auto &atom : atoms) {
    atom.halo_color = render3d::RgbaColor{.r = 0, .g = 0, .b = 0, .a = 0};
    atom.dim_factor = 1.0F;
  }
}

struct PhaseVisualContext {
  float min_x{0.0F};
  float max_x{0.0F};
  float min_y{0.0F};
  float max_y{0.0F};
  float diff_x{0.0F};
  std::uint8_t alpha_byte{0};
  float halo_radius_ratio{1.0F};
  float non_matching_dim{1.0F};
  plotters::ColormapType colormap{plotters::ColormapType::Viridis};
  bool is_2d{false};
};

[[nodiscard]] float computeNonMatchingDimFactor(const PhaseSearchSettings &settings) noexcept {
  if (settings.hide_non_matching) {
    return 0.0F;
  }
  return settings.dim_non_matching ? settings.dim_factor : 1.0F;
}

[[nodiscard]] float computeLargeCellHaloBoost(const render3d::Scene3D &scene) noexcept {
  const float box_diag = math::norm(scene.boundsMax() - scene.boundsMin());
  if (box_diag <= 40.0F) {
    return 1.0F;
  }
  return std::clamp(1.0F + 0.35F * std::log2(box_diag / 40.0F), 1.0F, 2.0F);
}

[[nodiscard]] bool isPhaseMatch(float val_x, const PhaseVisualContext &ctx,
                                std::span<const real_t> values_y, std::size_t index) noexcept {
  const bool match_x = (val_x >= ctx.min_x && val_x <= ctx.max_x);
  if (!match_x) {
    return false;
  }
  if (ctx.is_2d && !values_y.empty()) {
    const auto val_y = static_cast<float>(values_y[index]);
    return (val_y >= ctx.min_y && val_y <= ctx.max_y);
  }
  return true;
}

void applyAtomVisual(render3d::AtomInstance &atom, float val_x, bool matches,
                     const PhaseVisualContext &ctx) {
  if (matches) {
    const double norm_t =
        (ctx.diff_x > 1e-6F)
            ? std::clamp(static_cast<double>((val_x - ctx.min_x) / ctx.diff_x), 0.0, 1.0)
            : 0.5;
    const auto sampled = plotters::ColormapService::sample(norm_t, ctx.colormap);
    atom.halo_color = render3d::RgbaColor{
        .r = sampled.r,
        .g = sampled.g,
        .b = sampled.b,
        .a = ctx.alpha_byte,
    };
    atom.halo_radius_ratio = ctx.halo_radius_ratio;
    atom.dim_factor = 1.0F;
  } else {
    atom.halo_color = render3d::RgbaColor{.r = 0, .g = 0, .b = 0, .a = 0};
    atom.dim_factor = ctx.non_matching_dim;
  }
}

} // namespace

PhaseSearchService::PhaseSearchService(CellProvider cell_provider, SceneAccessor scene_accessor,
                                       ActiveFrameProvider active_frame_provider,
                                       TotalFramesProvider total_frames_provider,
                                       UpdateCallback on_phase_updated,
                                       std::shared_ptr<std::atomic<bool>> alive)
    : cell_provider_(std::move(cell_provider)), scene_accessor_(std::move(scene_accessor)),
      active_frame_provider_(std::move(active_frame_provider)),
      total_frames_provider_(std::move(total_frames_provider)),
      on_phase_updated_(std::move(on_phase_updated)), alive_(std::move(alive)) {
  phase_worker_ = std::jthread([this] { phaseWorkerLoop(); });
}

PhaseSearchService::~PhaseSearchService() {
  stop_phase_worker_.store(true);
  phase_cv_.notify_all();
  if (phase_worker_.joinable()) {
    phase_worker_.join();
  }
}

void PhaseSearchService::setEnabled(bool enabled) {
  const std::scoped_lock lock(phase_mutex_);
  phase_settings_.enabled = enabled;
}

bool PhaseSearchService::isEnabled() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.enabled;
}

void PhaseSearchService::setDescriptor(PhaseDescriptor descriptor) {
  const std::scoped_lock lock(phase_mutex_);
  phase_settings_.preset = CrystalPreset::Custom;
  phase_settings_.descriptor = descriptor;
  const auto [domain_min, domain_max] = getDescriptorDomain(descriptor);
  phase_settings_.min_threshold = domain_min;
  phase_settings_.max_threshold = domain_max;
}

PhaseDescriptor PhaseSearchService::descriptor() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.descriptor;
}

void PhaseSearchService::setCrystalPreset(CrystalPreset preset) {
  const std::scoped_lock lock(phase_mutex_);
  phase_settings_.preset = preset;
  if (preset != CrystalPreset::Custom) {
    const auto config = getPresetConfig(preset);
    if (phase_settings_.is_2d_mode && config.is_2d) {
      phase_settings_.descriptor = config.descriptor_2d_x;
      phase_settings_.min_threshold = config.min_threshold_2d_x;
      phase_settings_.max_threshold = config.max_threshold_2d_x;
      phase_settings_.descriptor_y = config.descriptor_2d_y;
      phase_settings_.min_threshold_y = config.min_threshold_2d_y;
      phase_settings_.max_threshold_y = config.max_threshold_2d_y;
    } else {
      phase_settings_.descriptor = config.descriptor;
      phase_settings_.min_threshold = config.min_threshold;
      phase_settings_.max_threshold = config.max_threshold;
      if (config.is_2d) {
        phase_settings_.descriptor_y = config.descriptor_2d_y;
        phase_settings_.min_threshold_y = config.min_threshold_2d_y;
        phase_settings_.max_threshold_y = config.max_threshold_2d_y;
      }
    }
  }
}

CrystalPreset PhaseSearchService::crystalPreset() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.preset;
}

void PhaseSearchService::setThresholds(float min_threshold, float max_threshold) {
  const std::scoped_lock lock(phase_mutex_);
  phase_settings_.preset = CrystalPreset::Custom;
  phase_settings_.min_threshold = std::min(min_threshold, max_threshold);
  phase_settings_.max_threshold = std::max(min_threshold, max_threshold);
}

float PhaseSearchService::minThreshold() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.min_threshold;
}

float PhaseSearchService::maxThreshold() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.max_threshold;
}

void PhaseSearchService::setColormap(plotters::ColormapType colormap) {
  const std::scoped_lock lock(phase_mutex_);
  phase_settings_.colormap = colormap;
}

plotters::ColormapType PhaseSearchService::colormap() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.colormap;
}

void PhaseSearchService::setHaloIntensity(float intensity) {
  const std::scoped_lock lock(phase_mutex_);
  phase_settings_.halo_intensity = std::clamp(intensity, 0.1F, 1.0F);
}

float PhaseSearchService::haloIntensity() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.halo_intensity;
}

void PhaseSearchService::setHaloRadiusRatio(float ratio) {
  const std::scoped_lock lock(phase_mutex_);
  phase_settings_.halo_radius_ratio = std::clamp(ratio, 1.0F, 5.0F);
}

float PhaseSearchService::haloRadiusRatio() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.halo_radius_ratio;
}

void PhaseSearchService::setAtomScale(float scale) {
  const std::scoped_lock lock(phase_mutex_);
  phase_settings_.atom_scale = std::clamp(scale, 0.05F, 2.0F);
}

float PhaseSearchService::atomScale() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.atom_scale;
}

void PhaseSearchService::setDimNonMatching(bool dim) {
  const std::scoped_lock lock(phase_mutex_);
  phase_settings_.dim_non_matching = dim;
}

bool PhaseSearchService::isDimNonMatching() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.dim_non_matching;
}

void PhaseSearchService::setHideNonMatching(bool hide) {
  const std::scoped_lock lock(phase_mutex_);
  phase_settings_.hide_non_matching = hide;
}

bool PhaseSearchService::isHideNonMatching() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.hide_non_matching;
}

void PhaseSearchService::set2DMode(bool is_2d) {
  const std::scoped_lock lock(phase_mutex_);
  phase_settings_.is_2d_mode = is_2d;
}

bool PhaseSearchService::is2DMode() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.is_2d_mode;
}

void PhaseSearchService::setDescriptorY(PhaseDescriptor desc_y) {
  const std::scoped_lock lock(phase_mutex_);
  phase_settings_.preset = CrystalPreset::Custom;
  phase_settings_.descriptor_y = desc_y;
  const auto [dmin, dmax] = getDescriptorDomain(desc_y);
  phase_settings_.min_threshold_y = dmin;
  phase_settings_.max_threshold_y = dmax;
}

PhaseDescriptor PhaseSearchService::descriptorY() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.descriptor_y;
}

void PhaseSearchService::setThresholds2D(float min_x, float max_x, float min_y, float max_y) {
  const std::scoped_lock lock(phase_mutex_);
  phase_settings_.preset = CrystalPreset::Custom;
  phase_settings_.min_threshold = std::min(min_x, max_x);
  phase_settings_.max_threshold = std::max(min_x, max_x);
  phase_settings_.min_threshold_y = std::min(min_y, max_y);
  phase_settings_.max_threshold_y = std::max(min_y, max_y);
}

float PhaseSearchService::minThresholdY() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.min_threshold_y;
}

float PhaseSearchService::maxThresholdY() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.max_threshold_y;
}

const PhaseSearchSettings &PhaseSearchService::settings() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_;
}

const PhaseHistogram &PhaseSearchService::latestHistogram() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return latest_histogram_;
}

const PhaseHeatmap2D &PhaseSearchService::latestHeatmap() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return latest_heatmap_;
}

const PhaseDescriptorsLRUCache &PhaseSearchService::lruCache() const noexcept { return lru_cache_; }

PhaseDescriptorsLRUCache &PhaseSearchService::lruCache() noexcept { return lru_cache_; }

void PhaseSearchService::clearCache() { lru_cache_.clear(); }

void PhaseSearchService::requestCalculation(std::size_t frame_index) {
  requested_phase_frame_.store(static_cast<int64_t>(frame_index));
  phase_cv_.notify_one();
}

void PhaseSearchService::applyVisualsToScene(render3d::Scene3D &scene, std::size_t active_frame) {
  const std::scoped_lock lock(phase_mutex_);
  const auto descriptors = lru_cache_.peek(active_frame);
  if (phase_settings_.enabled && !descriptors) {
    requestCalculation(active_frame);
  }
  applyPhaseVisualsToSceneInternal(scene, descriptors.get());
  latest_histogram_ = computePhaseHistogram(descriptors.get());
  latest_heatmap_ = computePhaseHeatmap(descriptors.get());
}

void PhaseSearchService::applyPhaseVisualsToSceneInternal(
    render3d::Scene3D &scene, const FramePhaseDescriptors *descriptors) {
  scene.setAtomScale(phase_settings_.atom_scale);
  auto &atoms = scene.atoms();
  if (!phase_settings_.enabled || descriptors == nullptr) {
    resetPhaseVisuals(atoms);
    return;
  }

  const auto values_x = extractDescriptorValues(*descriptors, phase_settings_.descriptor);
  if (values_x.size() != atoms.size()) {
    return;
  }

  const bool is_2d = phase_settings_.is_2d_mode;
  const auto values_y = is_2d ? extractDescriptorValues(*descriptors, phase_settings_.descriptor_y)
                              : std::span<const real_t>{};

  const float halo_alpha = std::clamp(phase_settings_.halo_intensity, 0.0F, 1.0F);
  const auto alpha_byte = static_cast<std::uint8_t>(std::round(halo_alpha * 255.0F));
  const float effective_halo_ratio =
      phase_settings_.halo_radius_ratio * computeLargeCellHaloBoost(scene);

  const PhaseVisualContext ctx{
      .min_x = phase_settings_.min_threshold,
      .max_x = phase_settings_.max_threshold,
      .min_y = phase_settings_.min_threshold_y,
      .max_y = phase_settings_.max_threshold_y,
      .diff_x = phase_settings_.max_threshold - phase_settings_.min_threshold,
      .alpha_byte = alpha_byte,
      .halo_radius_ratio = effective_halo_ratio,
      .non_matching_dim = computeNonMatchingDimFactor(phase_settings_),
      .colormap = phase_settings_.colormap,
      .is_2d = is_2d,
  };

  for (std::size_t i = 0; i < atoms.size(); ++i) {
    const auto val_x = static_cast<float>(values_x[i]);
    const bool matches = isPhaseMatch(val_x, ctx, values_y, i);
    applyAtomVisual(atoms[i], val_x, matches, ctx);
  }
}

std::shared_ptr<const FramePhaseDescriptors>
PhaseSearchService::computeDescriptorsSync(std::size_t frame_index) {
  if (!cell_provider_) {
    return nullptr;
  }

  auto maybe_cell_and_cutoffs = cell_provider_(frame_index);
  if (!maybe_cell_and_cutoffs) {
    return nullptr;
  }

  auto [target_cell, cutoffs] = std::move(*maybe_cell_and_cutoffs);
  if (target_cell.isEmpty() || target_cell.atomCount() == 0) {
    return nullptr;
  }

  if (cutoffs.empty() || cutoffs.size() != target_cell.elements().size()) {
    cutoffs = buildCellBondCutoffs(target_cell);
  }

  const real_t max_dist = computeMaxBondDistance(cutoffs, target_cell.elements().size());
  try {
    const analysis::StructureAnalyzer analyzer(target_cell, max_dist, cutoffs, false);
    const calculators::SteinhardtCalculator::SteinhardtOptions steinhardt_opts{
        .keep_q6m = true,
    };
    auto steinhardt_params = calculators::SteinhardtCalculator::calculatePerAtom(
        target_cell, &analyzer, steinhardt_opts);
    auto lef_values = calculators::LocalEntropyCalculator::calculatePerAtom(target_cell, &analyzer);
    const calculators::TenWoldeFrenkelParams twf_params{};
    auto twf_values = calculators::TenWoldeFrenkelCalculator::calculatePerAtom(
        target_cell, &analyzer, twf_params, &steinhardt_params);
    auto q_tet_values =
        calculators::TetrahedralCalculator::calculatePerAtom(target_cell, &analyzer);

    auto descriptors = std::make_shared<FramePhaseDescriptors>();
    descriptors->steinhardt =
        std::make_shared<const calculators::SteinhardtCalculator::SteinhardtParams>(
            std::move(steinhardt_params));
    descriptors->local_entropy = std::move(lef_values);
    descriptors->twf_nconn = std::move(twf_values);
    descriptors->q_tet = std::move(q_tet_values);
    return descriptors;
  } catch (const std::exception &) {
    return nullptr;
  }
}

std::shared_ptr<const FramePhaseDescriptors>
PhaseSearchService::getOrComputeDescriptors(std::size_t frame_index) {
  auto cached = lru_cache_.get(frame_index);
  if (cached) {
    return cached;
  }
  auto computed = computeDescriptorsSync(frame_index);
  if (computed) {
    lru_cache_.put(frame_index, computed);
  }
  return computed;
}

std::shared_ptr<const calculators::SteinhardtCalculator::SteinhardtParams>
PhaseSearchService::getOrComputeSteinhardt(std::size_t frame_index) {
  const auto descriptors = getOrComputeDescriptors(frame_index);
  return descriptors ? descriptors->steinhardt : nullptr;
}

PhaseHistogram
PhaseSearchService::computePhaseHistogram(const FramePhaseDescriptors *descriptors) const {
  PhaseHistogram hist;
  constexpr std::size_t K_BIN_COUNT = 50;
  hist.counts.assign(K_BIN_COUNT, 0);
  hist.bin_centers.resize(K_BIN_COUNT);

  const auto [domain_min, domain_max] = getDescriptorDomain(phase_settings_.descriptor);
  hist.min_val = domain_min;
  hist.max_val = domain_max;
  const float domain_span = std::max(domain_max - domain_min, 1e-6F);
  hist.bin_width = domain_span / static_cast<float>(K_BIN_COUNT);

  for (std::size_t i = 0; i < K_BIN_COUNT; ++i) {
    hist.bin_centers[i] = domain_min + (static_cast<float>(i) + 0.5F) * hist.bin_width;
  }

  if (descriptors == nullptr) {
    return hist;
  }

  const auto values = extractDescriptorValues(*descriptors, phase_settings_.descriptor);
  hist.total_atom_count = static_cast<std::uint32_t>(values.size());

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

PhaseHeatmap2D
PhaseSearchService::computePhaseHeatmap(const FramePhaseDescriptors *descriptors) const {
  constexpr std::size_t K_GRID = PhaseHeatmap2D::K_GRID_SIZE;
  PhaseHeatmap2D heatmap{};
  heatmap.counts.assign(K_GRID * K_GRID, 0);

  const auto [domain_min_x, domain_max_x] = getDescriptorDomain(phase_settings_.descriptor);
  const auto [domain_min_y, domain_max_y] = getDescriptorDomain(phase_settings_.descriptor_y);

  heatmap.min_x = domain_min_x;
  heatmap.max_x = domain_max_x;
  heatmap.min_y = domain_min_y;
  heatmap.max_y = domain_max_y;
  const float span_x = std::max(domain_max_x - domain_min_x, 1e-6F);
  const float span_y = std::max(domain_max_y - domain_min_y, 1e-6F);
  heatmap.bin_width_x = span_x / static_cast<float>(K_GRID);
  heatmap.bin_width_y = span_y / static_cast<float>(K_GRID);

  if (descriptors == nullptr) {
    return heatmap;
  }

  const auto values_x = extractDescriptorValues(*descriptors, phase_settings_.descriptor);
  const auto values_y = extractDescriptorValues(*descriptors, phase_settings_.descriptor_y);
  if (values_x.size() != values_y.size()) {
    return heatmap;
  }

  heatmap.total_atom_count = static_cast<std::uint32_t>(values_x.size());
  const float min_x_thresh = phase_settings_.min_threshold;
  const float max_x_thresh = phase_settings_.max_threshold;
  const float min_y_thresh = phase_settings_.min_threshold_y;
  const float max_y_thresh = phase_settings_.max_threshold_y;

  for (std::size_t i = 0; i < values_x.size(); ++i) {
    const auto v_x = static_cast<float>(values_x[i]);
    const auto v_y = static_cast<float>(values_y[i]);

    if (v_x >= min_x_thresh && v_x <= max_x_thresh && v_y >= min_y_thresh && v_y <= max_y_thresh) {
      ++heatmap.matching_atom_count;
    }

    const auto bin_x =
        static_cast<std::size_t>(std::clamp(std::floor((v_x - domain_min_x) / heatmap.bin_width_x),
                                            0.0F, static_cast<float>(K_GRID - 1)));
    const auto bin_y =
        static_cast<std::size_t>(std::clamp(std::floor((v_y - domain_min_y) / heatmap.bin_width_y),
                                            0.0F, static_cast<float>(K_GRID - 1)));

    const std::size_t cell_idx = bin_y * K_GRID + bin_x;
    if (cell_idx < heatmap.counts.size()) {
      ++heatmap.counts[cell_idx];
      heatmap.max_count = std::max(heatmap.max_count, heatmap.counts[cell_idx]);
    }
  }

  return heatmap;
}

int64_t PhaseSearchService::waitForNextPhaseFrame() {
  std::unique_lock lock(phase_mutex_);
  phase_cv_.wait(
      lock, [this] { return stop_phase_worker_.load() || requested_phase_frame_.load() >= 0; });
  if (stop_phase_worker_.load()) {
    return -1;
  }
  return requested_phase_frame_.exchange(-1);
}

void PhaseSearchService::processTargetPhaseFrame(int64_t target_frame) {
  if (target_frame < 0) {
    return;
  }
  const auto frame_idx = static_cast<std::size_t>(target_frame);
  auto descriptors = lru_cache_.get(frame_idx);
  if (!descriptors) {
    descriptors = computeDescriptorsSync(frame_idx);
    if (descriptors) {
      lru_cache_.put(frame_idx, descriptors);
    }
  }

  bool need_update = false;
  {
    const std::scoped_lock lock(phase_mutex_);
    const auto active_frame = active_frame_provider_ ? active_frame_provider_() : 0UL;
    if (active_frame == frame_idx && scene_accessor_) {
      const auto scene = scene_accessor_();
      if (scene != nullptr) {
        applyPhaseVisualsToSceneInternal(*scene, descriptors.get());
        latest_histogram_ = computePhaseHistogram(descriptors.get());
        latest_heatmap_ = computePhaseHeatmap(descriptors.get());
        need_update = true;
      }
    }
  }

  if (need_update && on_phase_updated_) {
    on_phase_updated_();
  }

  prewarmAdjacentFrames(frame_idx);
}

void PhaseSearchService::prewarmAdjacentFrames(std::size_t center_frame) {
  const std::size_t total = total_frames_provider_ ? total_frames_provider_() : 0UL;
  if (total <= 1) {
    return;
  }
  constexpr std::array<int64_t, 4> K_OFFSETS{-1, 1, -2, 2};
  for (const auto offset : K_OFFSETS) {
    if (stop_phase_worker_.load() || requested_phase_frame_.load() >= 0) {
      break;
    }
    const auto cand_idx = static_cast<int64_t>(center_frame) + offset;
    if (cand_idx >= 0 && std::cmp_less(cand_idx, total)) {
      const auto frame_idx = static_cast<std::size_t>(cand_idx);
      if (!lru_cache_.contains(frame_idx)) {
        auto pre_descriptors = computeDescriptorsSync(frame_idx);
        if (pre_descriptors) {
          lru_cache_.put(frame_idx, std::move(pre_descriptors));
        }
      }
    }
  }
}

void PhaseSearchService::phaseWorkerLoop() {
  while (!stop_phase_worker_.load()) {
    const auto target_frame = waitForNextPhaseFrame();
    if (target_frame < 0 || stop_phase_worker_.load()) {
      break;
    }
    processTargetPhaseFrame(target_frame);
  }
}

} // namespace correlation::app
