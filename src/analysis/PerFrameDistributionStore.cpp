// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "analysis/PerFrameDistributionStore.hpp"

#include <algorithm>
#include <mutex>
#include <utility>

namespace correlation::analysis {

const std::vector<real_t> PerFrameDistributionStore::EMPTY_CURVE{};

PerFrameDistributionStore::PerFrameDistributionStore(std::size_t total_frames)
    : total_frames_(total_frames) {
  // Initialize descriptors lazily: frame_curves vectors are allocated only when a metric is written
  metrics_.at(toIndex(DistributionMetric::RDF)) = {
      .bins = {},
      .mean_curve = {},
      .frame_curves = {},
      .name = "g(r)",
      .title = "Radial Distribution Function",
      .x_unit = "Å",
      .y_unit = "",
  };

  metrics_.at(toIndex(DistributionMetric::PAD)) = {
      .bins = {},
      .mean_curve = {},
      .frame_curves = {},
      .name = "PAD",
      .title = "PAD(θ) Normalized",
      .x_unit = "deg",
      .y_unit = "deg⁻¹",
  };

  metrics_.at(toIndex(DistributionMetric::SQ)) = {
      .bins = {},
      .mean_curve = {},
      .frame_curves = {},
      .name = "S(Q)",
      .title = "Structure Factor",
      .x_unit = "Å⁻¹",
      .y_unit = "",
  };

  metrics_.at(toIndex(DistributionMetric::XRD)) = {
      .bins = {},
      .mean_curve = {},
      .frame_curves = {},
      .name = "XRD",
      .title = "X-Ray Diffraction Pattern",
      .x_unit = "deg",
      .y_unit = "a.u.",
  };
}

std::size_t PerFrameDistributionStore::toIndex(DistributionMetric metric) noexcept {
  return static_cast<std::size_t>(metric);
}

void PerFrameDistributionStore::setBins(DistributionMetric metric, std::vector<real_t> bins) {
  const std::unique_lock lock(mutex_);
  metrics_.at(toIndex(metric)).bins = std::move(bins);
}

void PerFrameDistributionStore::setMeanCurve(DistributionMetric metric,
                                             std::vector<real_t> mean_curve) {
  const std::unique_lock lock(mutex_);
  metrics_.at(toIndex(metric)).mean_curve = std::move(mean_curve);
}

void PerFrameDistributionStore::setFrameCurve(DistributionMetric metric, std::size_t frame_idx,
                                              std::vector<real_t> curve) {
  if (frame_idx >= total_frames_) {
    return;
  }
  const std::unique_lock lock(mutex_);
  auto &frame_curves = metrics_.at(toIndex(metric)).frame_curves;
  if (frame_curves.empty()) {
    frame_curves.resize(total_frames_);
  }
  frame_curves[frame_idx] = std::move(curve);
}

std::size_t PerFrameDistributionStore::getTotalFrames() const noexcept { return total_frames_; }

bool PerFrameDistributionStore::hasData(DistributionMetric metric) const noexcept {
  const std::shared_lock lock(mutex_);
  const auto &data = metrics_.at(toIndex(metric));
  return !data.bins.empty() && (!data.mean_curve.empty() || !data.frame_curves.empty());
}

bool PerFrameDistributionStore::hasFrame(DistributionMetric metric,
                                         std::size_t frame_idx) const noexcept {
  const std::shared_lock lock(mutex_);
  const auto &frame_curves = metrics_.at(toIndex(metric)).frame_curves;
  if (frame_idx >= frame_curves.size()) {
    return false;
  }
  return !frame_curves[frame_idx].empty();
}

const std::vector<real_t> &PerFrameDistributionStore::getBins(DistributionMetric metric) const {
  const std::shared_lock lock(mutex_);
  return metrics_.at(toIndex(metric)).bins;
}

const std::vector<real_t> &
PerFrameDistributionStore::getMeanCurve(DistributionMetric metric) const {
  const std::shared_lock lock(mutex_);
  return metrics_.at(toIndex(metric)).mean_curve;
}

const std::vector<real_t> &PerFrameDistributionStore::getFrameCurve(DistributionMetric metric,
                                                                    std::size_t frame_idx) const {
  const std::shared_lock lock(mutex_);
  const auto &frame_curves = metrics_.at(toIndex(metric)).frame_curves;
  if (frame_idx >= frame_curves.size() || frame_curves[frame_idx].empty()) {
    return EMPTY_CURVE;
  }
  return frame_curves[frame_idx];
}

const MetricCurveData &PerFrameDistributionStore::getMetricData(DistributionMetric metric) const {
  const std::shared_lock lock(mutex_);
  return metrics_.at(toIndex(metric));
}

MetricCurveData &PerFrameDistributionStore::getMetricData(DistributionMetric metric) {
  const std::unique_lock lock(mutex_);
  return metrics_.at(toIndex(metric));
}

void PerFrameDistributionStore::configureMeanFrames(MeanFrameWindow window) {
  const std::unique_lock lock(mutex_);
  has_explicit_mean_mask_ = true;
  in_mean_mask_.assign(total_frames_, 0U);
  const std::size_t effective_stride = (window.stride > 0U) ? window.stride : 1U;
  for (std::size_t sample_idx = 0; sample_idx < window.sample_count; ++sample_idx) {
    const std::size_t frame_idx = window.start_frame + (sample_idx * effective_stride);
    if (frame_idx < total_frames_) {
      in_mean_mask_[frame_idx] = 1U;
    }
  }
}

void PerFrameDistributionStore::markFrameInMean(std::size_t frame_idx, bool in_mean) {
  if (frame_idx >= total_frames_) {
    return;
  }
  const std::unique_lock lock(mutex_);
  if (!has_explicit_mean_mask_) {
    has_explicit_mean_mask_ = true;
    in_mean_mask_.assign(total_frames_, 0U);
  }
  in_mean_mask_[frame_idx] = in_mean ? 1U : 0U;
}

bool PerFrameDistributionStore::isFrameInMean(std::size_t frame_idx) const noexcept {
  if (frame_idx >= total_frames_) {
    return false;
  }
  const std::shared_lock lock(mutex_);
  if (!has_explicit_mean_mask_) {
    return true;
  }
  return frame_idx < in_mean_mask_.size() && in_mean_mask_[frame_idx] != 0U;
}

void PerFrameDistributionStore::setFrameGenerator(FrameCurveGenerator generator) {
  const std::unique_lock lock(mutex_);
  frame_generator_ = std::move(generator);
}

void PerFrameDistributionStore::clearFrameGenerator() noexcept {
  const std::unique_lock lock(mutex_);
  frame_generator_ = nullptr;
}

bool PerFrameDistributionStore::needsFrameGenerationUnlocked(std::size_t frame_idx) const noexcept {
  if (frame_idx >= total_frames_) {
    return false;
  }
  return std::ranges::any_of(metrics_, [frame_idx](const MetricCurveData &data) {
    return !data.mean_curve.empty() &&
           (frame_idx >= data.frame_curves.size() || data.frame_curves[frame_idx].empty());
  });
}

bool PerFrameDistributionStore::needsFrameGeneration(std::size_t frame_idx) const noexcept {
  const std::shared_lock lock(mutex_);
  return needsFrameGenerationUnlocked(frame_idx);
}

void PerFrameDistributionStore::ensureFrame(std::size_t frame_idx) const {
  {
    const std::shared_lock lock(mutex_);
    if (!frame_generator_ || !needsFrameGenerationUnlocked(frame_idx)) {
      return;
    }
  }

  const std::scoped_lock gen_lock(generation_mutex_);
  FrameCurveGenerator generator;
  {
    const std::shared_lock lock(mutex_);
    if (!frame_generator_ || !needsFrameGenerationUnlocked(frame_idx)) {
      return;
    }
    generator = frame_generator_;
  }

  auto generated = generator(frame_idx);

  const std::unique_lock lock(mutex_);
  for (std::size_t idx = 0; idx < METRIC_COUNT; ++idx) {
    if (generated.at(idx).empty()) {
      continue;
    }
    auto &frame_curves = metrics_.at(idx).frame_curves;
    if (frame_curves.empty()) {
      frame_curves.resize(total_frames_);
    }
    if (frame_idx < frame_curves.size() && frame_curves[frame_idx].empty()) {
      frame_curves[frame_idx] = std::move(generated.at(idx));
    }
  }
}

} // namespace correlation::analysis
