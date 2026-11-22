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
#include <cmath>
#include <cstdint>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
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
  case PhaseDescriptor::Displacement:
    return descriptors.displacement;
  case PhaseDescriptor::D2Min:
    return descriptors.d2_min;
  }
  return {};
}

constexpr std::size_t K_LUT_SIZE = 256;
constexpr std::size_t K_PARALLEL_ATOM_THRESHOLD = 2048;
constexpr std::size_t K_CACHE_MEMORY_BUDGET_BYTES = 256UL * 1024UL * 1024UL;
constexpr std::size_t K_BYTES_PER_ATOM_ESTIMATE = 72UL;
constexpr std::size_t K_MIN_CACHE_FRAMES = 32UL;
constexpr std::size_t K_MAX_CACHE_FRAMES = 512UL;

[[nodiscard]] std::size_t computeAdaptiveCacheCapacity(std::size_t atom_count) noexcept {
  const std::size_t per_frame_bytes =
      std::max<std::size_t>(1UL, atom_count * K_BYTES_PER_ATOM_ESTIMATE);
  return std::clamp(K_CACHE_MEMORY_BUDGET_BYTES / per_frame_bytes, K_MIN_CACHE_FRAMES,
                    K_MAX_CACHE_FRAMES);
}

[[nodiscard]] bool isSteinhardtDescriptor(PhaseDescriptor descriptor) noexcept {
  return descriptor == PhaseDescriptor::Q6Bar || descriptor == PhaseDescriptor::Q4Bar ||
         descriptor == PhaseDescriptor::Q6 || descriptor == PhaseDescriptor::Q4 ||
         descriptor == PhaseDescriptor::W6Hat || descriptor == PhaseDescriptor::W4Hat ||
         descriptor == PhaseDescriptor::TwfNconn;
}

[[nodiscard]] bool hasSingleDescriptor(const FramePhaseDescriptors &desc,
                                       PhaseDescriptor descriptor) noexcept {
  if (isSteinhardtDescriptor(descriptor)) {
    if (descriptor == PhaseDescriptor::TwfNconn) {
      return !desc.twf_nconn.empty();
    }
    return desc.steinhardt != nullptr;
  }
  if (descriptor == PhaseDescriptor::LocalEntropy) {
    return !desc.local_entropy.empty();
  }
  if (descriptor == PhaseDescriptor::QTet) {
    return !desc.q_tet.empty();
  }
  if (descriptor == PhaseDescriptor::Displacement) {
    return !desc.displacement.empty();
  }
  if (descriptor == PhaseDescriptor::D2Min) {
    return !desc.d2_min.empty();
  }
  return false;
}

[[nodiscard]] bool hasRequiredDescriptors(const FramePhaseDescriptors &desc,
                                          const PhaseSearchSettings &settings) noexcept {
  if (!hasSingleDescriptor(desc, settings.descriptor)) {
    return false;
  }
  if (settings.is_2d_mode && !hasSingleDescriptor(desc, settings.descriptor_y)) {
    return false;
  }
  return true;
}

void resetPhaseVisuals(std::span<render3d::AtomInstance> atoms) noexcept {
  if (atoms.size() >= K_PARALLEL_ATOM_THRESHOLD) {
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, atoms.size()),
                      [&atoms](const tbb::blocked_range<std::size_t> &range) {
                        for (std::size_t i = range.begin(); i < range.end(); ++i) {
                          atoms[i].halo_color = render3d::RgbaColor{.r = 0, .g = 0, .b = 0, .a = 0};
                          atoms[i].dim_factor = 1.0F;
                        }
                      });
    return;
  }
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
  float inv_diff_x{0.0F};
  float halo_radius_ratio{1.0F};
  float non_matching_dim{1.0F};
  std::array<render3d::RgbaColor, K_LUT_SIZE> color_lut{};
  bool is_2d{false};
};

[[nodiscard]] std::array<render3d::RgbaColor, K_LUT_SIZE>
buildColormapLut(plotters::ColormapType colormap, std::uint8_t alpha_byte) noexcept {
  std::array<render3d::RgbaColor, K_LUT_SIZE> lut{};
  constexpr double K_INV_MAX = 1.0 / static_cast<double>(K_LUT_SIZE - 1);
  for (std::size_t i = 0; i < K_LUT_SIZE; ++i) {
    const auto sampled =
        plotters::ColormapService::sample(static_cast<double>(i) * K_INV_MAX, colormap);
    lut.at(i) = render3d::RgbaColor{
        .r = sampled.r,
        .g = sampled.g,
        .b = sampled.b,
        .a = alpha_byte,
    };
  }
  return lut;
}

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
                     const PhaseVisualContext &ctx) noexcept {
  if (matches) {
    const float norm_t = (ctx.inv_diff_x > 0.0F)
                             ? std::clamp((val_x - ctx.min_x) * ctx.inv_diff_x, 0.0F, 1.0F)
                             : 0.5F;
    const auto lut_idx = static_cast<std::size_t>(norm_t * static_cast<float>(K_LUT_SIZE - 1));
    atom.halo_color = ctx.color_lut.at(lut_idx);
    atom.halo_radius_ratio = ctx.halo_radius_ratio;
    atom.dim_factor = 1.0F;
  } else {
    atom.halo_color = render3d::RgbaColor{.r = 0, .g = 0, .b = 0, .a = 0};
    atom.dim_factor = ctx.non_matching_dim;
  }
}

} // namespace

PhaseSearchService::PhaseSearchService(CellProvider cell_provider,
                                       ActiveFrameProvider active_frame_provider,
                                       TotalFramesProvider total_frames_provider,
                                       UpdateCallback on_phase_updated,
                                       std::shared_ptr<std::atomic<bool>> alive)
    : cell_provider_(std::move(cell_provider)),
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
  {
    const std::scoped_lock lock(phase_mutex_);
    phase_settings_.enabled = enabled;
    if (enabled) {
      startPrebakeLocked();
    } else {
      cancelPrebakeLocked();
    }
  }
  phase_cv_.notify_one();
}

bool PhaseSearchService::isEnabled() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.enabled;
}

void PhaseSearchService::setDescriptor(PhaseDescriptor descriptor) {
  {
    const std::scoped_lock lock(phase_mutex_);
    phase_settings_.preset = CrystalPreset::Custom;
    phase_settings_.descriptor = descriptor;
    const auto [domain_min, domain_max] = getDescriptorDomain(descriptor);
    phase_settings_.min_threshold = domain_min;
    phase_settings_.max_threshold = domain_max;
    if (phase_settings_.enabled) {
      startPrebakeLocked();
    }
  }
  phase_cv_.notify_one();
}

PhaseDescriptor PhaseSearchService::descriptor() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.descriptor;
}

void PhaseSearchService::setCrystalPreset(CrystalPreset preset) {
  {
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
    if (phase_settings_.enabled) {
      startPrebakeLocked();
    }
  }
  phase_cv_.notify_one();
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
  {
    const std::scoped_lock lock(phase_mutex_);
    phase_settings_.is_2d_mode = is_2d;
    if (phase_settings_.enabled) {
      startPrebakeLocked();
    }
  }
  phase_cv_.notify_one();
}

bool PhaseSearchService::is2DMode() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_.is_2d_mode;
}

void PhaseSearchService::setDescriptorY(PhaseDescriptor desc_y) {
  {
    const std::scoped_lock lock(phase_mutex_);
    phase_settings_.preset = CrystalPreset::Custom;
    phase_settings_.descriptor_y = desc_y;
    const auto [dmin, dmax] = getDescriptorDomain(desc_y);
    phase_settings_.min_threshold_y = dmin;
    phase_settings_.max_threshold_y = dmax;
    if (phase_settings_.enabled) {
      startPrebakeLocked();
    }
  }
  phase_cv_.notify_one();
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

PhaseSearchSettings PhaseSearchService::settings() const noexcept {
  const std::scoped_lock lock(phase_mutex_);
  return phase_settings_;
}

PhaseHistogram PhaseSearchService::latestHistogram() const {
  const std::scoped_lock lock(phase_mutex_);
  return latest_histogram_;
}

PhaseHeatmap2D PhaseSearchService::latestHeatmap() const {
  const std::scoped_lock lock(phase_mutex_);
  return latest_heatmap_;
}

const PhaseDescriptorsLRUCache &PhaseSearchService::lruCache() const noexcept { return lru_cache_; }

PhaseDescriptorsLRUCache &PhaseSearchService::lruCache() noexcept { return lru_cache_; }

void PhaseSearchService::setGpuDescriptorEvaluator(GpuDescriptorEvaluator evaluator) {
  const std::scoped_lock lock(phase_mutex_);
  gpu_evaluator_ = std::move(evaluator);
}

void PhaseSearchService::clearCache() {
  cancelPrebake();
  {
    const std::scoped_lock lock(phase_mutex_);
    ref_frame0_analyzer_.reset();
  }
  lru_cache_.clear();
}

void PhaseSearchService::startPrebakeLocked() { prebake_active_.store(true); }

void PhaseSearchService::cancelPrebakeLocked() {
  prebake_active_.store(false);
  is_computing_.store(false);
}

void PhaseSearchService::startPrebake() {
  {
    const std::scoped_lock lock(phase_mutex_);
    if (!phase_settings_.enabled) {
      return;
    }
    startPrebakeLocked();
  }
  phase_cv_.notify_one();
}

void PhaseSearchService::cancelPrebake() {
  const std::scoped_lock lock(phase_mutex_);
  cancelPrebakeLocked();
}

bool PhaseSearchService::isFrameReady(std::size_t frame_index) const {
  PhaseSearchSettings active_settings;
  {
    const std::scoped_lock lock(phase_mutex_);
    active_settings = phase_settings_;
  }
  if (!active_settings.enabled) {
    return true;
  }
  const auto cached = lru_cache_.peek(frame_index);
  return cached != nullptr && hasRequiredDescriptors(*cached, active_settings);
}

PhasePrebakeStatus PhaseSearchService::prebakeStatus() const {
  PhaseSearchSettings active_settings;
  {
    const std::scoped_lock lock(phase_mutex_);
    active_settings = phase_settings_;
  }
  const std::size_t total = total_frames_provider_ ? total_frames_provider_() : 0UL;
  if (!active_settings.enabled || total == 0) {
    return PhasePrebakeStatus{};
  }
  const std::size_t target_total = std::min(total, lru_cache_.capacity());
  const std::size_t cached =
      lru_cache_.countMatching([&active_settings](const FramePhaseDescriptors &desc) {
        return hasRequiredDescriptors(desc, active_settings);
      });
  const float progress =
      (target_total > 0)
          ? std::clamp(static_cast<float>(cached) / static_cast<float>(target_total), 0.0F, 1.0F)
          : 1.0F;
  return PhasePrebakeStatus{
      .progress = progress,
      .cached_frames = cached,
      .total_frames = target_total,
      .is_computing = is_computing_.load() || (prebake_active_.load() && cached < target_total),
  };
}

void PhaseSearchService::requestCalculation(std::size_t frame_index) {
  requested_phase_frame_.store(static_cast<int64_t>(frame_index));
  phase_cv_.notify_one();
}

void PhaseSearchService::applyVisualsToScene(render3d::Scene3D &scene, std::size_t active_frame) {
  const std::scoped_lock lock(phase_mutex_);
  const auto descriptors = lru_cache_.peek(active_frame);
  if (phase_settings_.enabled &&
      (!descriptors || !hasRequiredDescriptors(*descriptors, phase_settings_))) {
    requestCalculation(active_frame);
  }
  applyPhaseVisualsToSceneInternal(scene, descriptors.get());
  if (!phase_settings_.enabled || descriptors != nullptr) {
    if (phase_settings_.is_2d_mode) {
      latest_heatmap_ = computePhaseHeatmap(descriptors.get());
    } else {
      latest_histogram_ = computePhaseHistogram(descriptors.get());
    }
  }
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
  const float diff_x = phase_settings_.max_threshold - phase_settings_.min_threshold;

  const PhaseVisualContext ctx{
      .min_x = phase_settings_.min_threshold,
      .max_x = phase_settings_.max_threshold,
      .min_y = phase_settings_.min_threshold_y,
      .max_y = phase_settings_.max_threshold_y,
      .inv_diff_x = (diff_x > 1e-6F) ? (1.0F / diff_x) : 0.0F,
      .halo_radius_ratio = effective_halo_ratio,
      .non_matching_dim = computeNonMatchingDimFactor(phase_settings_),
      .color_lut = buildColormapLut(phase_settings_.colormap, alpha_byte),
      .is_2d = is_2d,
  };

  if (atoms.size() >= K_PARALLEL_ATOM_THRESHOLD) {
    tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, atoms.size()),
        [&atoms, &values_x, &values_y, &ctx](const tbb::blocked_range<std::size_t> &range) {
          for (std::size_t i = range.begin(); i < range.end(); ++i) {
            const auto val_x = static_cast<float>(values_x[i]);
            const bool matches = isPhaseMatch(val_x, ctx, values_y, i);
            applyAtomVisual(atoms[i], val_x, matches, ctx);
          }
        });
    return;
  }

  for (std::size_t i = 0; i < atoms.size(); ++i) {
    const auto val_x = static_cast<float>(values_x[i]);
    const bool matches = isPhaseMatch(val_x, ctx, values_y, i);
    applyAtomVisual(atoms[i], val_x, matches, ctx);
  }
}

namespace {

[[nodiscard]] math::Vector3<double> minimumImageDiffDouble(const math::Vector3<real_t> &origin,
                                                           const math::Vector3<real_t> &target,
                                                           const core::Cell &cell) noexcept {
  const math::Vector3<double> diff(
      static_cast<double>(target.x()) - static_cast<double>(origin.x()),
      static_cast<double>(target.y()) - static_cast<double>(origin.y()),
      static_cast<double>(target.z()) - static_cast<double>(origin.z()));
  if (cell.volume() <= static_cast<real_t>(1e-6)) {
    return diff;
  }
  const auto &inv = cell.inverseLatticeVectors();
  const auto &lat = cell.latticeVectors();
  const double fx = static_cast<double>(inv[0][0]) * diff.x() +
                    static_cast<double>(inv[0][1]) * diff.y() +
                    static_cast<double>(inv[0][2]) * diff.z();
  const double fy = static_cast<double>(inv[1][0]) * diff.x() +
                    static_cast<double>(inv[1][1]) * diff.y() +
                    static_cast<double>(inv[1][2]) * diff.z();
  const double fz = static_cast<double>(inv[2][0]) * diff.x() +
                    static_cast<double>(inv[2][1]) * diff.y() +
                    static_cast<double>(inv[2][2]) * diff.z();

  const double wx = fx - std::round(fx);
  const double wy = fy - std::round(fy);
  const double wz = fz - std::round(fz);

  return math::Vector3<double>(
      static_cast<double>(lat[0][0]) * wx + static_cast<double>(lat[1][0]) * wy +
          static_cast<double>(lat[2][0]) * wz,
      static_cast<double>(lat[0][1]) * wx + static_cast<double>(lat[1][1]) * wy +
          static_cast<double>(lat[2][1]) * wz,
      static_cast<double>(lat[0][2]) * wx + static_cast<double>(lat[1][2]) * wy +
          static_cast<double>(lat[2][2]) * wz);
}

[[nodiscard]] std::vector<real_t> computePerAtomDisplacement(const core::Cell &target_cell,
                                                             const core::Cell &ref_cell) {
  const auto &target_atoms = target_cell.atoms();
  const auto &ref_atoms = ref_cell.atoms();
  const std::size_t count = target_atoms.size();
  std::vector<real_t> disp(count, static_cast<real_t>(0.0));
  if (ref_atoms.size() != count || count == 0) {
    return disp;
  }
  const bool has_pbc = target_cell.volume() > static_cast<real_t>(1e-6) &&
                       ref_cell.volume() > static_cast<real_t>(1e-6);
  const auto &inv_t = target_cell.inverseLatticeVectors();
  const auto &inv_0 = ref_cell.inverseLatticeVectors();
  const auto &lat_t = target_cell.latticeVectors();

  auto compute_one = [&](std::size_t idx) {
    const auto &pos_t = target_atoms[idx].position();
    const auto &pos_0 = ref_atoms[idx].position();
    if (!has_pbc) {
      const double dx = static_cast<double>(pos_t.x()) - static_cast<double>(pos_0.x());
      const double dy = static_cast<double>(pos_t.y()) - static_cast<double>(pos_0.y());
      const double dz = static_cast<double>(pos_t.z()) - static_cast<double>(pos_0.z());
      disp[idx] = static_cast<real_t>(std::sqrt(dx * dx + dy * dy + dz * dz));
      return;
    }
    const auto frac_t = inv_t * pos_t;
    const auto frac_0 = inv_0 * pos_0;
    const double wx = static_cast<double>(frac_t.x()) - static_cast<double>(frac_0.x()) -
                      std::round(static_cast<double>(frac_t.x()) - static_cast<double>(frac_0.x()));
    const double wy = static_cast<double>(frac_t.y()) - static_cast<double>(frac_0.y()) -
                      std::round(static_cast<double>(frac_t.y()) - static_cast<double>(frac_0.y()));
    const double wz = static_cast<double>(frac_t.z()) - static_cast<double>(frac_0.z()) -
                      std::round(static_cast<double>(frac_t.z()) - static_cast<double>(frac_0.z()));
    const double dx = static_cast<double>(lat_t[0][0]) * wx +
                      static_cast<double>(lat_t[1][0]) * wy + static_cast<double>(lat_t[2][0]) * wz;
    const double dy = static_cast<double>(lat_t[0][1]) * wx +
                      static_cast<double>(lat_t[1][1]) * wy + static_cast<double>(lat_t[2][1]) * wz;
    const double dz = static_cast<double>(lat_t[0][2]) * wx +
                      static_cast<double>(lat_t[1][2]) * wy + static_cast<double>(lat_t[2][2]) * wz;
    disp[idx] = static_cast<real_t>(std::sqrt(dx * dx + dy * dy + dz * dz));
  };

  if (count >= K_PARALLEL_ATOM_THRESHOLD) {
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, count),
                      [&compute_one](const tbb::blocked_range<std::size_t> &range) {
                        for (std::size_t idx = range.begin(); idx < range.end(); ++idx) {
                          compute_one(idx);
                        }
                      });
  } else {
    for (std::size_t idx = 0; idx < count; ++idx) {
      compute_one(idx);
    }
  }
  return disp;
}

[[nodiscard]] bool invert3x3Double(const std::array<std::array<double, 3>, 3> &mat,
                                   std::array<std::array<double, 3>, 3> &inv) noexcept {
  const double c00 = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
  const double c01 = mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2];
  const double c02 = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];
  const double det = mat[0][0] * c00 + mat[0][1] * c01 + mat[0][2] * c02;
  if (std::abs(det) < 1e-10) {
    return false;
  }
  const double inv_det = 1.0 / det;
  inv[0][0] = c00 * inv_det;
  inv[1][0] = c01 * inv_det;
  inv[2][0] = c02 * inv_det;
  inv[0][1] = (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]) * inv_det;
  inv[1][1] = (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) * inv_det;
  inv[2][1] = (mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1]) * inv_det;
  inv[0][2] = (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) * inv_det;
  inv[1][2] = (mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]) * inv_det;
  inv[2][2] = (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) * inv_det;
  return true;
}

[[nodiscard]] double computeSingleAtomD2Min(std::size_t atom_idx, const core::Cell &target_cell,
                                            const core::NeighborGraph &ref_graph) noexcept {
  const auto neighbors = ref_graph.getNeighbors(atom_idx);
  if (neighbors.size() < 4) {
    return 0.0;
  }
  const auto &target_atoms = target_cell.atoms();
  const auto &origin_pos = target_atoms[atom_idx].position();

  std::array<std::array<double, 3>, 3> x_mat{};
  std::array<std::array<double, 3>, 3> y_mat{};
  for (const auto &nbr : neighbors) {
    if (nbr.index >= target_atoms.size()) {
      continue;
    }
    const std::array<double, 3> ref_bond{static_cast<double>(nbr.r_ij.x()),
                                         static_cast<double>(nbr.r_ij.y()),
                                         static_cast<double>(nbr.r_ij.z())};
    const auto cur_vec =
        minimumImageDiffDouble(origin_pos, target_atoms[nbr.index].position(), target_cell);
    const std::array<double, 3> cur_bond{cur_vec.x(), cur_vec.y(), cur_vec.z()};
    for (std::size_t row = 0; row < 3; ++row) {
      for (std::size_t col = 0; col < 3; ++col) {
        x_mat.at(row).at(col) += cur_bond.at(row) * ref_bond.at(col);
        y_mat.at(row).at(col) += ref_bond.at(row) * ref_bond.at(col);
      }
    }
  }

  std::array<std::array<double, 3>, 3> y_inv{};
  std::array<std::array<double, 3>, 3> j_mat{{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
  if (invert3x3Double(y_mat, y_inv)) {
    j_mat = {};
    for (std::size_t row = 0; row < 3; ++row) {
      for (std::size_t col = 0; col < 3; ++col) {
        for (std::size_t mid = 0; mid < 3; ++mid) {
          j_mat.at(row).at(col) += x_mat.at(row).at(mid) * y_inv.at(mid).at(col);
        }
      }
    }
  }

  double d2_sum = 0.0;
  for (const auto &nbr : neighbors) {
    if (nbr.index >= target_atoms.size()) {
      continue;
    }
    const std::array<double, 3> ref_bond{static_cast<double>(nbr.r_ij.x()),
                                         static_cast<double>(nbr.r_ij.y()),
                                         static_cast<double>(nbr.r_ij.z())};
    const auto cur_vec =
        minimumImageDiffDouble(origin_pos, target_atoms[nbr.index].position(), target_cell);
    const double res_x = cur_vec.x() - (j_mat[0][0] * ref_bond[0] + j_mat[0][1] * ref_bond[1] +
                                        j_mat[0][2] * ref_bond[2]);
    const double res_y = cur_vec.y() - (j_mat[1][0] * ref_bond[0] + j_mat[1][1] * ref_bond[1] +
                                        j_mat[1][2] * ref_bond[2]);
    const double res_z = cur_vec.z() - (j_mat[2][0] * ref_bond[0] + j_mat[2][1] * ref_bond[1] +
                                        j_mat[2][2] * ref_bond[2]);
    d2_sum += res_x * res_x + res_y * res_y + res_z * res_z;
  }
  return d2_sum;
}

[[nodiscard]] std::vector<real_t> computePerAtomD2Min(const core::Cell &target_cell,
                                                      const core::NeighborGraph &ref_graph) {
  const std::size_t count = target_cell.atoms().size();
  std::vector<real_t> d2_vals(count, static_cast<real_t>(0.0));
  if (count == 0) {
    return d2_vals;
  }
  if (count >= K_PARALLEL_ATOM_THRESHOLD) {
    tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, count),
        [&d2_vals, &target_cell, &ref_graph](const tbb::blocked_range<std::size_t> &range) {
          for (std::size_t idx = range.begin(); idx < range.end(); ++idx) {
            d2_vals[idx] = static_cast<real_t>(computeSingleAtomD2Min(idx, target_cell, ref_graph));
          }
        });
  } else {
    for (std::size_t idx = 0; idx < count; ++idx) {
      d2_vals[idx] = static_cast<real_t>(computeSingleAtomD2Min(idx, target_cell, ref_graph));
    }
  }
  return d2_vals;
}

void populateMissingDescriptors(FramePhaseDescriptors &descriptors, const core::Cell &target_cell,
                                const PhaseSearchSettings &settings) {
  const auto &analyzer = *descriptors.analyzer;
  const auto active_desc = settings.descriptor;
  const bool is_2d = settings.is_2d_mode;
  const auto active_desc_y = settings.descriptor_y;

  const bool need_twf = (active_desc == PhaseDescriptor::TwfNconn) ||
                        (is_2d && active_desc_y == PhaseDescriptor::TwfNconn);
  const bool need_steinhardt =
      isSteinhardtDescriptor(active_desc) || (is_2d && isSteinhardtDescriptor(active_desc_y));

  if (need_steinhardt && (!descriptors.steinhardt ||
                          (need_twf && !descriptors.steinhardt->q6m_normalized.has_value()))) {
    const calculators::SteinhardtCalculator::SteinhardtOptions steinhardt_opts{
        .keep_q6m = need_twf,
    };
    auto computed = calculators::SteinhardtCalculator::calculatePerAtom(target_cell, &analyzer,
                                                                        steinhardt_opts);
    descriptors.steinhardt =
        std::make_shared<const calculators::SteinhardtCalculator::SteinhardtParams>(
            std::move(computed));
  }

  const bool need_entropy = (active_desc == PhaseDescriptor::LocalEntropy) ||
                            (is_2d && active_desc_y == PhaseDescriptor::LocalEntropy);
  if (need_entropy && descriptors.local_entropy.empty()) {
    descriptors.local_entropy =
        calculators::LocalEntropyCalculator::calculatePerAtom(target_cell, &analyzer);
  }

  if (need_twf && descriptors.twf_nconn.empty()) {
    const calculators::TenWoldeFrenkelParams twf_params{};
    descriptors.twf_nconn = calculators::TenWoldeFrenkelCalculator::calculatePerAtom(
        target_cell, &analyzer, twf_params, descriptors.steinhardt.get());
  }

  const bool need_qtet =
      (active_desc == PhaseDescriptor::QTet) || (is_2d && active_desc_y == PhaseDescriptor::QTet);
  if (need_qtet && descriptors.q_tet.empty()) {
    descriptors.q_tet =
        calculators::TetrahedralCalculator::calculatePerAtom(target_cell, &analyzer);
  }
}

} // namespace

std::shared_ptr<const analysis::StructureAnalyzer> PhaseSearchService::getOrBuildFrame0Analyzer() {
  {
    const std::scoped_lock lock(phase_mutex_);
    if (ref_frame0_analyzer_) {
      return ref_frame0_analyzer_;
    }
  }
  auto ref_opt = cell_provider_(0);
  if (!ref_opt.has_value()) {
    return nullptr;
  }
  auto [ref_cell, ref_cutoffs] = std::move(*ref_opt);
  if (ref_cell.isEmpty() || ref_cell.atomCount() == 0) {
    return nullptr;
  }
  if (ref_cutoffs.empty() || ref_cutoffs.size() != ref_cell.elements().size()) {
    ref_cutoffs = buildCellBondCutoffs(ref_cell);
  }
  const real_t ref_max_dist = computeMaxBondDistance(ref_cutoffs, ref_cell.elements().size());
  auto built = std::make_shared<const analysis::StructureAnalyzer>(
      ref_cell, ref_max_dist, ref_cutoffs, false, static_cast<real_t>(0.0));
  const std::scoped_lock lock(phase_mutex_);
  ref_frame0_analyzer_ = built;
  return built;
}

void PhaseSearchService::populateDeformationDescriptors(FramePhaseDescriptors &descriptors,
                                                        std::size_t frame_index,
                                                        const core::Cell &target_cell,
                                                        const PhaseSearchSettings &settings) {
  const bool need_disp =
      (settings.descriptor == PhaseDescriptor::Displacement) ||
      (settings.is_2d_mode && settings.descriptor_y == PhaseDescriptor::Displacement);
  const bool need_d2min = (settings.descriptor == PhaseDescriptor::D2Min) ||
                          (settings.is_2d_mode && settings.descriptor_y == PhaseDescriptor::D2Min);

  if (need_disp && descriptors.displacement.empty()) {
    if (frame_index == 0) {
      descriptors.displacement.assign(target_cell.atomCount(), static_cast<real_t>(0.0));
    } else if (auto ref_opt = cell_provider_(0); ref_opt.has_value()) {
      descriptors.displacement = computePerAtomDisplacement(target_cell, ref_opt->first);
    }
  }

  if (need_d2min && descriptors.d2_min.empty()) {
    if (frame_index == 0) {
      descriptors.d2_min.assign(target_cell.atomCount(), static_cast<real_t>(0.0));
    } else if (const auto ref_analyzer = getOrBuildFrame0Analyzer()) {
      descriptors.d2_min = computePerAtomD2Min(target_cell, ref_analyzer->neighborGraph());
    }
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

  PhaseSearchSettings active_settings;
  {
    const std::scoped_lock lock(phase_mutex_);
    active_settings = phase_settings_;
  }

  const auto existing = lru_cache_.peek(frame_index);
  auto descriptors = existing ? std::make_shared<FramePhaseDescriptors>(*existing)
                              : std::make_shared<FramePhaseDescriptors>();

  try {
    GpuDescriptorEvaluator gpu_eval;
    {
      const std::scoped_lock lock(phase_mutex_);
      gpu_eval = gpu_evaluator_;
    }

    if (cutoffs.empty() || cutoffs.size() != target_cell.elements().size()) {
      cutoffs = buildCellBondCutoffs(target_cell);
    }
    const real_t max_dist = computeMaxBondDistance(cutoffs, target_cell.elements().size());

    // GPU-accelerated pre-bake for large cells (>= 4k atoms) on scalar GPU-supported descriptors
    constexpr std::size_t K_GPU_PREBAKE_THRESHOLD = 4096;
    if (gpu_eval && target_cell.atomCount() >= K_GPU_PREBAKE_THRESHOLD) {
      auto try_gpu = [&](PhaseDescriptor desc) {
        std::vector<real_t> vals;
        if (desc == PhaseDescriptor::LocalEntropy && descriptors->local_entropy.empty() &&
            gpu_eval(target_cell, desc, max_dist, vals) && vals.size() == target_cell.atomCount()) {
          descriptors->local_entropy = std::move(vals);
        } else if (desc == PhaseDescriptor::QTet && descriptors->q_tet.empty() &&
                   gpu_eval(target_cell, desc, max_dist, vals) &&
                   vals.size() == target_cell.atomCount()) {
          descriptors->q_tet = std::move(vals);
        } else if (desc == PhaseDescriptor::TwfNconn && descriptors->twf_nconn.empty() &&
                   gpu_eval(target_cell, desc, max_dist, vals) &&
                   vals.size() == target_cell.atomCount()) {
          descriptors->twf_nconn = std::move(vals);
        }
      };
      try_gpu(active_settings.descriptor);
      if (active_settings.is_2d_mode) {
        try_gpu(active_settings.descriptor_y);
      }
    }

    if (!hasRequiredDescriptors(*descriptors, active_settings)) {
      if (!descriptors->analyzer) {
        // Pass r_bin_width = 0.0 to skip allocating/filling unused RDF raw_histograms_ tensor
        descriptors->analyzer = std::make_shared<const analysis::StructureAnalyzer>(
            target_cell, max_dist, cutoffs, false, static_cast<real_t>(0.0));
      }
      populateMissingDescriptors(*descriptors, target_cell, active_settings);
      populateDeformationDescriptors(*descriptors, frame_index, target_cell, active_settings);
    }

    // Release heavy StructureAnalyzer neighbor graphs so scalar descriptors fit 32..512 frames
    descriptors->analyzer.reset();
    lru_cache_.setCapacity(computeAdaptiveCacheCapacity(target_cell.atomCount()));
    return descriptors;
  } catch (const std::exception &) {
    return nullptr;
  }
}

std::shared_ptr<const FramePhaseDescriptors>
PhaseSearchService::getOrComputeDescriptors(std::size_t frame_index) {
  auto cached = lru_cache_.get(frame_index);
  PhaseSearchSettings active_settings;
  {
    const std::scoped_lock lock(phase_mutex_);
    active_settings = phase_settings_;
  }
  if (cached && hasRequiredDescriptors(*cached, active_settings)) {
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

int64_t
PhaseSearchService::findNextUncachedPrebakeFrame(const PhaseSearchSettings &settings) const {
  const std::size_t total = total_frames_provider_ ? total_frames_provider_() : 0UL;
  if (total == 0) {
    return -1;
  }
  const std::size_t cap = std::min(total, lru_cache_.capacity());
  const std::size_t start = active_frame_provider_ ? (active_frame_provider_() % total) : 0UL;
  for (std::size_t offset = 0; offset < cap; ++offset) {
    const std::size_t cand = (start + offset) % total;
    const auto cached = lru_cache_.peek(cand);
    if (!cached || !hasRequiredDescriptors(*cached, settings)) {
      return static_cast<int64_t>(cand);
    }
  }
  return -1;
}

int64_t PhaseSearchService::waitForNextPhaseFrame() {
  while (!stop_phase_worker_.load()) {
    PhaseSearchSettings active_settings;
    {
      std::unique_lock lock(phase_mutex_);
      phase_cv_.wait(lock, [this] {
        return stop_phase_worker_.load() || requested_phase_frame_.load() >= 0 ||
               prebake_active_.load();
      });
      if (stop_phase_worker_.load()) {
        return -1;
      }
      if (const int64_t req = requested_phase_frame_.exchange(-1); req >= 0) {
        return req;
      }
      active_settings = phase_settings_;
    }

    if (!prebake_active_.load() || !active_settings.enabled) {
      prebake_active_.store(false);
      is_computing_.store(false);
      continue;
    }

    const int64_t next_frame = findNextUncachedPrebakeFrame(active_settings);
    if (next_frame >= 0) {
      return next_frame;
    }

    prebake_active_.store(false);
    is_computing_.store(false);
    if (on_phase_updated_) {
      on_phase_updated_();
    }
  }
  return -1;
}

void PhaseSearchService::processTargetPhaseFrame(int64_t target_frame) {
  if (target_frame < 0) {
    return;
  }
  is_computing_.store(true);
  const auto frame_idx = static_cast<std::size_t>(target_frame);
  PhaseSearchSettings active_settings;
  {
    const std::scoped_lock lock(phase_mutex_);
    active_settings = phase_settings_;
  }
  auto descriptors = lru_cache_.get(frame_idx);
  if (!descriptors || !hasRequiredDescriptors(*descriptors, active_settings)) {
    descriptors = computeDescriptorsSync(frame_idx);
    if (descriptors) {
      lru_cache_.put(frame_idx, descriptors);
    }
  }

  const std::size_t active_frame = active_frame_provider_ ? active_frame_provider_() : 0UL;
  if (descriptors && frame_idx == active_frame) {
    const std::scoped_lock lock(phase_mutex_);
    if (phase_settings_.is_2d_mode) {
      latest_heatmap_ = computePhaseHeatmap(descriptors.get());
    } else {
      latest_histogram_ = computePhaseHistogram(descriptors.get());
    }
  }

  is_computing_.store(prebake_active_.load() || requested_phase_frame_.load() >= 0);
  if (on_phase_updated_) {
    on_phase_updated_();
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
